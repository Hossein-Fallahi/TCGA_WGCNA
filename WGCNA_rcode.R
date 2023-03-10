# This tutorial is modified from the https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html

# 1: Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")

#this line was added to address a warning:>> Please consider using the new package "gprofiler2". At the moment you are using a deprecated package relying on outdated data.
install.packages("gprofiler2")
library("gprofiler2")

# 2: TCGA data

GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

TCGAbiolinks:::getProjectSummary("TCGA-LIHC")
query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts")
lihc_res = getResults(query_TCGA) # make results as table

colnames(lihc_res) # columns present in the table
head(lihc_res$sample_type) # first 6 types of tissue.
summary(factor(lihc_res$sample_type)) # summary of distinct tissues types present in this study
query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))
GDCdownload(
  query = query_TCGA
)

tcga_data = GDCprepare(query_TCGA)


dim(tcga_data)
colnames(colData(tcga_data))


table(tcga_data@colData$vital_status)
table(tcga_data@colData$ajcc_pathologic_stage)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)

dim(assay(tcga_data))     # gene expression matrices.
head(assay(tcga_data)[,1:10]) 
head(rowData(tcga_data))     

saveRDS(object = tcga_data,
        file = "tcga_data.RDS",
        compress = FALSE)


tcga_data = readRDS(file = "tcga_data.RDS")


# 3: RNAseq data analysis
limma_pipeline = function(
    tcga_data,
    condition_variable,
    reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}

limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)


gender_limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="gender",
  reference_group="female"
)

plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}

res_pca = plot_PCA(limma_res$voomObj, "definition")
colnames(colData(tcga_data))


# 4: Classification of data


d_mat = as.matrix(t(limma_res$voomObj$E))
d_resp = as.factor(limma_res$voomObj$targets$definition)

set.seed(42)
train_ids = createDataPartition(d_resp, p=0.75, list=FALSE)

x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

y_train = d_resp[train_ids]
y_test  = d_resp[-train_ids]


res = cv.glmnet(
  x = x_train,
  y = y_train,
  alpha = 0.5,
  family = "binomial"
)

y_pred = predict(res, newx=x_test, type="class", s="lambda.min")

confusion_matrix = table(y_pred, y_test)

print(confusion_matrix)

print(paste0("Sensitivity: ",sensitivity(confusion_matrix)))
print(paste0("Specificity: ",specificity(confusion_matrix)))
print(paste0("Precision: ",precision(confusion_matrix)))


res_coef = coef(res, s="lambda.min") 
dim(res_coef)
head(res_coef) # in a sparse matrix the "." represents the value of zero
res_coef = res_coef[res_coef[,1] != 0,]
head(res_coef)
res_coef = res_coef[-1]

relevant_genes = names(res_coef) # get names of the (non-zero) variables.
length(relevant_genes) # number of selected genes
head(relevant_genes) # few select genes
head(limma_res$voomObj$genes)

relevant_gene_names = limma_res$voomObj$genes[relevant_genes,"gene_name"]
head(relevant_gene_names) # few select genes (with readable names now)

print(intersect(limma_res$topGenes$gene_id, relevant_genes))

# define the color palette for the plot
hmcol = colorRampPalette(rev(brewer.pal(9, "RdBu")))(256)


# perform complete linkage clustering
clust = function(x) hclust(x, method="complete")
# use the inverse of correlation as distance.
dist = function(x) as.dist((1-cor(t(x)))/2)

# Show green color for genes that also show up in DE analysis
colorLimmaGenes = ifelse(
  # Given a vector of boolean values
  (relevant_genes %in% limma_res$topGenes$ensembl_gene_id),
  "green", # if true, return green for that value
  "white" # if false, return white for that value
)




# As you've seen a good looking heatmap involves a lot of parameters
gene_heatmap = heatmap.2(
  t(d_mat[,relevant_genes]),
  scale="row",          # scale the values for each gene (row)
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  col=hmcol,            # define the color map
  labRow=relevant_gene_names, # use gene names instead of ensembl annotation
  RowSideColors=colorLimmaGenes,
  labCol=FALSE,         # Not showing column labels
  ColSideColors=as.character(as.numeric(d_resp)), # Show colors for each response class
  dendrogram="both",    # Show dendrograms for both axis
  hclust = clust,       # Define hierarchical clustering method
  distfun = dist,       # Using correlation coefficient for distance function
  cexRow=.6,            # Resize row labels
  margins=c(1,5)        # Define margin spaces
)

# Using the same method as in Day-2, get the dendrogram from the heatmap
# and cut it to get the 2 classes of genes

# Extract the hierarchical cluster from heatmap to class "hclust"
hc = as.hclust(gene_heatmap$rowDendrogram)

# Cut the tree into 2 groups, up-regulated in tumor and up-regulated in control
clusters = cutree(hc, k=2)
table(clusters)


# selecting just a few columns so that its easier to visualize the table
gprofiler_cols = c("significant","p.value","overlap.size","term.id","term.name")

# make sure the URL uses https
set_base_url("https://biit.cs.ut.ee/gprofiler")


gprofiler(names(clusters[clusters %in% 1]))[, gprofiler_cols]
gprofiler(names(clusters[clusters %in% 2]))[, gprofiler_cols]


# retain only a small subset of the genes (see documentation for ?varFilter)
d_mat = varFilter(limma_res$voomObj$E, var.func=IQR, var.cutoff=0.95)

# transpose the matrix, so that it has the same shape as the d_mat we used at the beginning
d_mat = t(d_mat)

print(dim(d_mat))




# size before
print(dim(x_train))
print(dim(x_test))

x_train = d_mat[train_ids, ]
x_test  = d_mat[-train_ids, ]

# size after
print(dim(x_train))
print(dim(x_test))


# 5: Survival Analysis

# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")

# NB: make sure you set the working directory of RStudio correctly


tcga_data = readRDS(file = "tcga_data.RDS")
limma_res = readRDS(file = "limma_res.RDS")

# extract clinical data
clinical = tcga_data@colData

dim(clinical)

# we are only interested in the "Primary solid Tumor" cases for survival
clin_df = clinical[clinical$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "ajcc_pathologic_stage")]

# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"
# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive


clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)


Surv(clin_df$overall_survival, clin_df$deceased)


Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender

# fit a survival model
fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)
print(fit)

# we produce a Kaplan Meier plot

ggsurvplot(fit, data=clin_df)


ggsurvplot(fit, data=clin_df, pval=T)

ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, risk.table.col="strata")

#Another question could be: how does tumor stage affect survival?
clin_df$ajcc_pathologic_stage = gsub("[abc]$", "", clin_df$ajcc_pathologic_stage)

# we remove those with stage "not reported", since they are unknown
clin_df[which(clin_df$ajcc_pathologic_stage == "not reported"), "ajcc_pathologic_stage"] = NA

# finally, we also remove those with tumor stage 4, since they are too few
clin_df[which(clin_df$ajcc_pathologic_stage == "stage iv"), "ajcc_pathologic_stage"] = NA

table(clin_df$ajcc_pathologic_stage)

fit = survfit(Surv(overall_survival, deceased) ~ ajcc_pathologic_stage, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)

# we produce a Kaplan-Meier plot from the fitted model
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T)


# 6: Gene expression and survival


# let's extract the table of differential expression we got earlier
expr_df = limma_res$topGenes

# print the first row, to see the gene name, the logFC value and the p-value
print(expr_df[1, ])

# get the ensembl gene id of the first row
gene_id = expr_df[1, "gene_id"]

# also get the common gene name of the first row
gene_name = expr_df[1, "gene_name"]

# visualize the gene expression distribution on the diseased samples (in black)
# versus the healthy samples (in red)
expr_diseased = d_mat[rownames(clin_df), gene_id]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")

# get the expression values for the selected gene
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]

# find the median value of the gene and print it
median_value = median(clin_df$gene_value)
print(median_value)

# divide patients in two groups, up and down regulated.
# if the patient expression is greater or equal to them median we put it
# among the "up-regulated", otherwise among the "down-regulated"

clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased) ~ gene, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)

# and finally, we produce a Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, pval=T, risk.table=T, title=paste(gene_name))

