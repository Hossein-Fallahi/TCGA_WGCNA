# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("clusterProfiler")
# BiocManager::install("EDASeq")
# BiocManager::install("edgeR")
# BiocManager::install("DO.db")
# BiocManager::install("GO.db")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("pathview")
# install.packages("here")

library(TCGAbiolinks)


barcodes1 <- c("TCGA-5L-AAT1", "TCGA-UU-A93S","TCGA-A8-A07W", "TCGA-AN-A0FJ", "TCGA-AC-A62V", "TCGA-BH-A18J", 
               "TCGA-B6-A0I9", "TCGA-AO-A0J5", "TCGA-B6-A0IB", "TCGA-A8-A08T", "TCGA-A2-A0SW", "TCGA-LL-A73Z", 
               "TCGA-A2-A0CS", "TCGA-PL-A8LX", "TCGA-A2-A0T2", "TCGA-A8-A08J", "TCGA-A2-A0SV", "TCGA-A8-A08O",
               "TCGA-BH-A1FH", "TCGA-B6-A3ZX")

query.exp.BRCAiv <- GDCquery(project = "TCGA-BRCA",
                             data.category = "Gene expression",
                             data.type = "Gene expression quantification",
                             platform = "Illumina HiSeq", 
                             file.type  = "results", 
                             sample.type = c("Primary Tumor"),
                             legacy = TRUE,
                             barcode=barcodes1)


GDCdownload(query.exp.BRCAiv)

barcodes2 <- c ("TCGA-BH-A0B6",
                "TCGA-GM-A2D9",
                "TCGA-BH-A0HA",
                "TCGA-BH-A18P",
                "TCGA-AR-A2LE",
                "TCGA-BH-A0W7",
                "TCGA-E2-A1IN",
                "TCGA-E2-A1LH",
                "TCGA-GM-A2DH",
                "TCGA-GM-A2DO",
                "TCGA-A8-A095",
                "TCGA-BH-A0BP",
                "TCGA-E2-A152",
                "TCGA-BH-A0H5",
                "TCGA-BH-A209",
                "TCGA-BH-A0WA",
                "TCGA-BH-A0BW",
                "TCGA-Z7-A8R6",
                "TCGA-BH-A0B8",
                "TCGA-E2-A14Z",
                "TCGA-A8-A06O",
                "TCGA-BH-A0AV",
                "TCGA-E2-A1IH",
                "TCGA-AR-A255",
                "TCGA-BH-A0B0",
                "TCGA-BH-A0HY",
                "TCGA-OL-A66J",
                "TCGA-AR-A2LR",
                "TCGA-AR-A1AK",
                "TCGA-AR-A1AY",
                "TCGA-GM-A2DI",
                "TCGA-GM-A2DD",
                "TCGA-B6-A40B",
                "TCGA-BH-A0BR",
                "TCGA-BH-A0DX",
                "TCGA-BH-A0BL",
                "TCGA-BH-A18K",
                "TCGA-AO-A03V",
                "TCGA-E2-A154",
                "TCGA-A2-A3XZ",
                "TCGA-AR-A24S",
                "TCGA-GM-A2DL",
                "TCGA-AR-A1AJ",
                "TCGA-EW-A1J6",
                "TCGA-BH-A1FD",
                "TCGA-A2-A0YF",
                "TCGA-E2-A15C",
                "TCGA-EW-A1IY",
                "TCGA-A7-A0CD",
                "TCGA-AR-A1AP",
                "TCGA-AR-A24N",
                "TCGA-A8-A0AD",
                "TCGA-E2-A15O",
                "TCGA-BH-A0BG",
                "TCGA-S3-AA14",
                "TCGA-E2-A1II",
                "TCGA-GM-A2DK",
                "TCGA-BH-A0DO",
                "TCGA-AR-A1AX",
                "TCGA-BH-A1FG",
                "TCGA-A1-A0SE",
                "TCGA-BH-A18S",
                "TCGA-BH-A0H6",
                "TCGA-BH-A0BO",
                "TCGA-AR-A24P",
                "TCGA-BH-A1ET",
                "TCGA-BH-A1EU",
                "TCGA-E2-A156",
                "TCGA-B6-A1KI",
                "TCGA-B6-A0X0",
                "TCGA-E2-A15J",
                "TCGA-BH-A0BQ",
                "TCGA-A2-A259",
                "TCGA-E2-A1IF",
                "TCGA-A2-A0EP",
                "TCGA-BH-A0C3",
                "TCGA-E2-A1IJ",
                "TCGA-BH-A0H3",
                "TCGA-E2-A1IO",
                "TCGA-A2-A0YI",
                "TCGA-AR-A252",
                "TCGA-A8-A08A",
                "TCGA-B6-A402",
                "TCGA-AO-A03M",
                "TCGA-E2-A14U",
                "TCGA-A7-A3IY",
                "TCGA-AO-A03U",
                "TCGA-E2-A15F",
                "TCGA-E2-A14S",
                "TCGA-A1-A0SB"
)

query.exp.BRCAi <- GDCquery(project = "TCGA-BRCA",
                            data.category = "Gene expression",
                            data.type = "Gene expression quantification",
                            platform = "Illumina HiSeq", 
                            file.type  = "results", 
                            sample.type = c("Primary Tumor"),
                            legacy = TRUE,
                            barcode=barcodes2)


GDCdownload(query.exp.BRCAi)


BRCAi.exp <- GDCprepare(query.exp.BRCAi, 
                        save = TRUE, 
                        summarizedExperiment = TRUE, 
                        save.filename = "BRCAillumina_HiSeq.rda")

BRCAiv.exp <- GDCprepare(query.exp.BRCAiv, 
                         save = TRUE, 
                         summarizedExperiment = TRUE, 
                         save.filename = "BRCAivIllumina_HiSeq.rda")


dataPrep_BRCAi <- TCGAanalyze_Preprocessing(object = BRCAi.exp,
                                            cor.cut = 0.6,    
                                            datatype = "raw_count",
                                            filename = "BRCAi_IlluminaHiSeq_RNASeqV2.png")

dataPrep_BRCAiv <- TCGAanalyze_Preprocessing(object = BRCAiv.exp,
                                             cor.cut = 0.6, 
                                             datatype = "raw_count",
                                             filename = "BRCAiv_IlluminaHiSeq_RNASeqV2.png")

dataNorm <- TCGAanalyze_Normalization(tabDF = cbind(dataPrep_BRCAi, dataPrep_BRCAiv),
                                      geneInfo = TCGAbiolinks::geneInfo,
                                      method = "gcContent") 

knitr::include_graphics(here::here("BRCAiv_IlluminaHiSeq_RNASeqV2.png"))


dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.05)  


save(dataFilt, file = paste0("BRCAi_BRCAiv_Norm_IlluminaHiSeq.rda"))


dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "filter2")  


save(dataFilt, file = paste0("BRCAi_BRCAiv_Norm_IlluminaHiSeq.rda"))

dataFiltBRCAi <- subset(dataFilt, select = substr(colnames(dataFilt),1,12) %in% barcodes2)
dataFiltBRCAiv <- subset(dataFilt, select = substr(colnames(dataFilt),1,12) %in% barcodes1)


dataDEGs <- TCGAanalyze_DEA(mat1 = dataFiltBRCAi,
                            mat2 = dataFiltBRCAiv,
                            Cond1type = "BRCAi",
                            Cond2type = "BRCAiv",
                            fdr.cut = 0.1,
                            method = "exactTest")

nrow(dataDEGs)


ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes BRCAi Vs BRCAiv", RegulonList = rownames(dataDEGs))


TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP, 
                        nRGTab = rownames(dataDEGs),
                        nBar = 20,
                        filename="BP.pdf")

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResCC),
                        GOCCTab = ansEA$ResCC,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20,
                        filename="CC.pdf")

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResPat),
                        PathTab = ansEA$ResPat,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20,
                        filename="Pat.pdf")

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResMF),
                        GOMFTab = ansEA$ResMF, 
                        nRGTab = rownames(dataDEGs),
                        nBar = 20,
                        filename="MF.pdf")

library(SummarizedExperiment)

library(clusterProfiler)


dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"BRCAi","BRCAiv",
                                          dataFilt[,colnames(dataFiltBRCAi)],
                                          dataFilt[,colnames(dataFiltBRCAiv)])

dataDEGsFiltLevel$GeneID <- 0


eg = as.data.frame(bitr(dataDEGsFiltLevel$mRNA,
                        fromType="SYMBOL",
                        toType="ENTREZID",
                        OrgDb="org.Hs.eg.db"))


eg <- eg[!duplicated(eg$SYMBOL),]

dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]

dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]


all(eg$SYMBOL == dataDEGsFiltLevel$mRNA)

dataDEGsFiltLevel$GeneID <- eg$ENTREZID

dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID



output <- names(genelistDEGs)
write.csv(output, "gene_list.csv")

data("bods", package = "pathview")
hsa03050 <- pathview::pathview(gene.data  = genelistDEGs,
                               pathway.id = "hsa03050",
                               species    = "hsa",
                               limit = list(gene=as.integer(max(abs(genelistDEGs)))))

hsa05204 <- pathview::pathview(gene.data  = genelistDEGs,
                               pathway.id = "hsa05204",
                               species    = "hsa",
                               limit = list(gene=as.integer(max(abs(genelistDEGs)))))

knitr::include_graphics(here::here("hsa03050.pathview.png"))
