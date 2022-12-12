#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
source("packages.R")

# Gene expression
gdc.ge<- GDCquery(project = "TCGA-BRCA",data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", barcode = c("TCGA-AR-A252-01A-11R-A169-07", "TCGA-D8-A1JP-01A-11R-A13Q-07","TCGA-A2-A0ER-01A-21R-A034-07","TCGA-AR-A2LM-01A-11R-A180-07","TCGA-AN-A0AK-01A-21R-A00Z-07","TCGA-AR-A0U0-01A-11R-A109-07","TCGA-AR-A2LH-01A-31R-A18M-07","TCGA-E2-A1L7-01A-11R-A144-0","TCGA-GI-A2C8-01A-11R-A16F-07","TCGA-E2-A1LE-01A-12R-A19W-07","TCGA-A2-A0SV-01A-11R-A084-07","TCGA-A8-A099-01A-11R-A00Z-07","TCGA-AR-A0TZ-01A-12R-A084-07"))
GDCdownload(gdc.ge, files.per.chunk = 200)



# Protein expression
gdc.pe<- GDCquery(project = "TCGA-BRCA",data.category = "Proteome Profiling",data.type = "Protein Expression Quantification", barcode = c("TCGA-AR-A252-01A", "TCGA-D8-A1JP-01A","TCGA-A2-A0ER-01A","TCGA-AR-A2LM-01A","TCGA-AN-A0AK-01A","TCGA-AR-A0U0-01A","TCGA-AR-A2LH-01A","TCGA-E2-A1L7-01A","TCGA-GI-A2C8-01A","TCGA-E2-A1LE-01A","TCGA-A2-A0SV-01A","TCGA-A8-A099-01A","TCGA-AR-A0TZ-01A"))
GDCdownload(gdc.pe, files.per.chunk = 200)



# Collecting Sample Sheet Information
SampleSheet.information=c("id","file_name","data_category","data_type", "project","cases.submitter_id","sample.submitter_id","sample_type")
SampleSheet=NULL
id=NULL
CollectSamplesheetInfo=getResults(gdc.ge)
L=ncol(CollectSamplesheetInfo)
for(i in 1:L)    id=c(id,which(colnames(CollectSamplesheetInfo)==SampleSheet.information[i]))
SampleSheet=CollectSamplesheetInfo[,id]
id=NULL
CollectSamplesheetInfo=getResults(gdc.pe)
L=ncol(CollectSamplesheetInfo)
for(i in 1:L)    id=c(id,which(colnames(CollectSamplesheetInfo)==SampleSheet.information[i]))
SampleSheet=rbind(SampleSheet,CollectSamplesheetInfo[,id])
colnames(SampleSheet)<-paste(c("File ID","File Name","Data Category","Data Type","Project ID","Case ID","Sample ID","Sample Type"))
write.table(SampleSheet,"samplesheet_file.txt",quote=F,col.names=T,row.names=F,sep="\t")


# Collecting clinical file
gdc.clinical <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", barcode = c("TCGA-AR-A252", "TCGA-D8-A1JP","TCGA-A2-A0ER","TCGA-AR-A2LM","TCGA-AN-A0AK","TCGA-AR-A0U0","TCGA-AR-A2LH-01A","TCGA-E2-A1L7","TCGA-GI-A2C8","TCGA-E2-A1LE","TCGA-A2-A0SV","TCGA-A8-A099","TCGA-AR-A0TZ"))
GDCdownload(gdc.clinical, files.per.chunk = 200)
clinical_InfoFile <- GDCprepare_clinic(gdc.clinical,"patient")
write.table(clinical_InfoFile,"clinical_InfoFile.txt",quote=F,col.names=T,row.names=F,sep="\t")


# Download Antibody map
TCGA_antibodies_descriptions.gencode.v36 <- read.csv(url("https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300"),sep="\t")

antibody_map="https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300"


samplesheet="./samplesheet_file.txt"
sampledata_folder="./GDCdata"
clinical_samplesheet="./clinical_InfoFile.txt"
antibody_map="https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300"


 A=read.csv(samplesheet,sep="\t")
 B=read.csv(clinical_samplesheet,sep="\t")

 samplesheet.data.category.column=which(colnames(A)=="Data.Category")
 samplesheet.sample.type.column=which(colnames(A)=="Sample.Type")
 samplesheet.case.id.column=which(colnames(A)=="Case.ID")
 samplesheet.data.type.column=which(colnames(A)=="Data.Type")
 clin.samplesheet.case.id.column=which(colnames(B)=="bcr_patient_barcode") #case_submitter_id
 clin.samplesheet.stage.column=which(colnames(B)=="stage_event_pathologic_stage") #ajcc_pathologic_stage
 clin.samplesheet.age.column=which(colnames(B)=="days_to_birth") #age_at_index


create_Omic1_inputfiles=Omic1_rmNoOmic2(samplesheet,data.type1="Protein Expression Quantification", data.type2="Transcriptome Profiling",sample.type="Primary Tumor",sampledata_folder,clinical_samplesheet,pattern="RPPA",samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.case.id.column,samplesheet.data.type.column,clin.samplesheet.case.id.column,clin.samplesheet.stage.column,clin.samplesheet.age.column,omic1.expr.column=6)


pheno="./status.txt"
subject.id.0="./subid.all.0.txt"
subject.id.1="./subid.all.1.txt"

create_Omic2_inputfiles=omic2(samplesheet,sampledata_folder,pheno,subject.id.0,subject.id.1,data.type1="Protein Expression Quantification",data.type2="Transcriptome Profiling",sample.type="Primary Tumor",gene.type="protein_coding",pattern = "star_gene_counts",gene.type.column=3,omic2.expr.column=8,samplesheet.data.category.column, samplesheet.sample.type.column, samplesheet.case.id.column,samplesheet.data.type.column)


pathway.genes="./pathway_genes.csv"
pheno="./status.txt"
omic1.expr_Group1="./omic1_Group1.txt"
omic1.expr_Group2="./omic1_Group2.txt"
get.pathway.omic1=pathwayOmic1_expr(pathway.genes, antibody_map,omic1.expr_Group1,omic1.expr_Group2)

omic2.expr_Group1="./omic2_Group1.txt"
omic2.expr_Group2="./omic2_Group2.txt"
get.pathway.omic2=pathwayOmic2_expr(pathway.genes,omic2.expr_Group1,omic2.expr_Group2)


 pathway_omic1.exprfile = "./pathway_omic1.expr.csv"
 pathway_omic2.exprfile = "./pathway_omic2.expr.csv"
 Result <- iosearchR(pathway.genes, pheno, pathway_omic1.exprfile, pathway_omic2.exprfile, K=1, ntopOmics2=5, ntopOmics1=3)
