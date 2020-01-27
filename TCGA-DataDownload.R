library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

##Download RNA-Seq read counts from TCGA projects. Below is a project of Pancreatic Cancer. Check TCGAbiolinks for more information.
query.exp <- GDCquery(project = "TCGA-PAAD", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")

GDCdownload(query.exp)
gene.count=GDCprepare(query=query.exp, save=T,save.filename="geneExpPAAD.rda")

##Extract gene raw read counts
gs=assay(gene.count)

##Normalize read counts using function vst provided by DESeq2
gs_normalized=vst(gs,blind=F)
