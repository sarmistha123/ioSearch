#' Finds significant pathways, omic1 and omic2
#'
#' @param pathway.genes Read the information obtained from find_pathwayGenes().
#' @param pheno Reads the phenotype file generated from Omic1_rmNA().
#' @param pathway_omic1_expr Reads matrix with pathway information and omic1 for all subjects generated from pathwayOmic1_expr().
#' @param pathway_omic2_expr Reads matrix with pathway information and omic2 for all subjects generated from pathwayOmic2_expr().
#' @param K Default is 1. Tuning parameter specific to elasticnet package.
#' @param ntopOmics2 Number of top omics1 variables to be selected from each pathway (user-defined).
#' @param ntopOmics1 Number of top omics1 variables to be selected from each pathway (user-defined).
#'
#' @return Matrices for significant pathways, omic1 and omic2.
#' @export iosearchR
#' @importFrom BiocManager install
#' @importFrom TCGAbiolinks GDCquery GDCdownload getResults GDCprepare_clinic
#' @examples
#'
#' ### Download a small breast cancer data from TCGA (GDC portal)
#'
#' ### Run the following codes to download ~26Mb data from TCGA to find ioSearch outcome:
#'
#' # BiocManager::install("TCGAbiolinks",force = TRUE)
#' # library(TCGAbiolinks)
#' # cat("Downloading small dataset containing gene and protein expression and clinical
#' # information from TCGA...","\n")
#'
#' ### First, downloading Gene Expression from TCGA
#'
#' # gdc.ge<- TCGAbiolinks::GDCquery(project = "TCGA-BRCA",data.category = "Transcriptome Profiling",
#' # data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts",
#' # barcode = c("TCGA-AR-A252-01A-11R-A169-07", "TCGA-D8-A1JP-01A-11R-A13Q-07",
#' # "TCGA-A2-A0ER-01A-21R-A034-07","TCGA-GI-A2C8-01A-11R-A16F-07",
#' # "TCGA-E2-A1LE-01A-12R-A19W-07", "TCGA-A2-A0SV-01A-11R-A084-07"))
#'
#' # TCGAbiolinks::GDCdownload(gdc.ge, files.per.chunk = 200)
#'
#' ### Second, downloading Protein Expression from TCGA
#'
#' # gdc.pe<- TCGAbiolinks::GDCquery(project = "TCGA-BRCA",data.category = "Proteome Profiling",
#' # data.type = "Protein Expression Quantification",
#' # barcode = c("TCGA-AR-A252-01A", "TCGA-D8-A1JP-01A","TCGA-A2-A0ER-01A",
#' # "TCGA-GI-A2C8-01A","TCGA-E2-A1LE-01A","TCGA-A2-A0SV-01A"))
#'
#' # TCGAbiolinks::GDCdownload(gdc.pe, files.per.chunk = 200)
#'
#' ### Now, collate SampleSheet Information using TCGAbiolinks or Download from TCGA directly
#'
#' # SampleSheet.information=c("id","file_name","data_category","data_type", "project",
#' # "cases.submitter_id","sample.submitter_id","sample_type")
#' # SampleSheet=id=NULL
#' # CollectSamplesheetInfo=TCGAbiolinks::getResults(gdc.ge)
#' # L=ncol(CollectSamplesheetInfo)
#' # for(i in 1:L) id=c(id,which(colnames(CollectSamplesheetInfo)==SampleSheet.information[i]))
#' # SampleSheet=CollectSamplesheetInfo[,id]
#'
#' # id=NULL
#' # CollectSamplesheetInfo=TCGAbiolinks::getResults(gdc.pe)
#' # L=ncol(CollectSamplesheetInfo)
#' # for(i in 1:L) id=c(id,which(colnames(CollectSamplesheetInfo)==SampleSheet.information[i]))
#' # SampleSheet=rbind(SampleSheet,CollectSamplesheetInfo[,id])
#' # colnames(SampleSheet)<-paste(c("File ID","File Name","Data Category","Data Type","Project ID",
#' # "Case ID","Sample ID","Sample Type"))
#' # write.table(SampleSheet,"samplesheet_file.txt",quote=FALSE,col.names=TRUE,
#' # row.names=FALSE,sep="\t")
#'
#' ### Next, collecting clinical information form TCGA
#' # gdc.clinical <- TCGAbiolinks::GDCquery(project = "TCGA-BRCA", data.category = "Clinical",
#' # barcode = c("TCGA-AR-A252", "TCGA-D8-A1JP","TCGA-A2-A0ER","TCGA-GI-A2C8",
#' # "TCGA-E2-A1LE","TCGA-A2-A0SV"))
#' # TCGAbiolinks::GDCdownload(gdc.clinical, files.per.chunk = 200)
#' # clinical_InfoFile <- TCGAbiolinks::GDCprepare_clinic(gdc.clinical,"patient")
#' # write.table(clinical_InfoFile,"clinical_InfoFile.txt",quote=FALSE,
#' # col.names=TRUE,row.names=FALSE,sep="\t")
#'
#' # cat("Downloading data is complete","\n")
#' # cat("Starting to collect required information..","\n")
#'
#' # samplesheet="./samplesheet_file.txt"
#' # sampledata_folder="./GDCdata"
#' # clinical_samplesheet="./clinical_InfoFile.txt"
#' # antibody_map="https://api.gdc.cancer.gov/v0/data/62647302-b4d3-4a81-a7c0-d141f5dbd300"
#'
#' # cat("Reading the collect required information from the downloads..","\n")
#' # A=read.csv(samplesheet,sep="\t")
#' # B=read.csv(clinical_samplesheet,sep="\t")
#'
#' # samplesheet.data.category.column=which(colnames(A)=="Data.Category")
#' # samplesheet.sample.type.column=which(colnames(A)=="Sample.Type")
#' # samplesheet.case.id.column=which(colnames(A)=="Case.ID")
#' # samplesheet.data.type.column=which(colnames(A)=="Data.Type")
#' # clin.samplesheet.case.id.column=which(colnames(B)=="bcr_patient_barcode") #case_submitter_id
#' # ## ajcc_pathologic_stage
#' # clin.samplesheet.stage.column=which(colnames(B)=="stage_event_pathologic_stage")
#' # clin.samplesheet.age.column=which(colnames(B)=="days_to_birth") #age_at_index
#'
#' ### Run Omic1_rmNoOmic2() function to create omic1 data matrix for both disease groups by collecting
#' ### information from all subject-specific omic1 files.
#'
#' # create_Omic1_inputfiles=Omic1_rmNoOmic2(samplesheet,data.type1="Protein Expression
#' # Quantification", data.type2="Transcriptome Profiling",sample.type="Primary Tumor",
#' # sampledata_folder, clinical_samplesheet,pattern="RPPA",samplesheet.sample.type.column,
#' # samplesheet.data.category.column,samplesheet.case.id.column,
#' # samplesheet.data.type.column, clin.samplesheet.case.id.column,clin.samplesheet.stage.column,
#' # clin.samplesheet.age.column,omic1.expr.column=6)
#'
#' ### Create the phenotype file
#'
#' # pheno=dstatus(samplesheet, clinical_samplesheet,data.type1="Protein Expression Quantification",
#' # samplesheet.data.type.column,samplesheet.case.id.column,
#' # clin.samplesheet.case.id.column, clin.samplesheet.stage.column,clin.samplesheet.age.column)
#'
### Create the Subject IDs file for both groups
#'
#' # subject.id.0=create_Omic1_inputfiles$subid.all.0
#' # subject.id.1=create_Omic1_inputfiles$subid.all.1
#'
#' ### Run omic2() function to create omic2 data matrix for both disease groups by collecting
#' ### information from all subject-specific omic2 files.
#'
#' # create_Omic2_inputfiles=omic2(samplesheet,sampledata_folder,pheno,subject.id.0,subject.id.1,
#' # data.type1="Protein Expression Quantification",data.type2="Transcriptome Profiling",
#' # sample.type="Primary Tumor",gene.type="protein_coding",pattern = "star_gene_counts",
#' # gene.type.column=3,omic2.expr.column=8,samplesheet.data.category.column,
#' # samplesheet.sample.type.column, samplesheet.case.id.column,samplesheet.data.type.column)
#'
#' ### Extracting pathway information for each gene from KEGG database
#'
#' # library(AnnotationDbi)
#' # library(KEGGREST)
#' # library(org.Hs.eg.db)
#' # find_pathwayGenes()
#'
#' # pathway.genes=find_pathwayGenes()
#'
#' ### Collate data from KEGG and omics-expression
#'
#' # omic1.expr_Group1=create_Omic1_inputfiles$omic1_Group1
#' # omic1.expr_Group2=create_Omic1_inputfiles$omic1_Group2
#'
#' # omic2.expr_Group1=create_Omic2_inputfiles$omic2_Group1
#' # omic2.expr_Group2=create_Omic2_inputfiles$omic2_Group2
#'
#' # get.pathway.omic1=pathwayOmic1_expr(pathway.genes, antibody_map,omic1.expr_Group1,
#' # omic1.expr_Group2)
#' # get.pathway.omic2=pathwayOmic2_expr(pathway.genes,omic2.expr_Group1,omic2.expr_Group2)
#'
#' # cat("Collecting required information complete","\n")
#'
#' # ioSearch_Result=iosearchR(pathway.genes, pheno, get.pathway.omic1,
#' # get.pathway.omic2, K=1, ntopOmics2=5, ntopOmics1=3)
#'
#' # DO NOT RUN
#' # if(.Platform$OS.type == "windows") withAutoprint({
#' ## gc()
#' # memory.size()
#' # memory.size(TRUE)
#' # memory.limit(35000)
#' # })
iosearchR=function(pathway.genes, pheno, pathway_omic1_expr, pathway_omic2_expr, K=1, ntopOmics2, ntopOmics1)
{
    cat("Reading all data files...","\n")
    #rdata_omic2=utils::read.csv(pathway_omic2_expr)
    #rdata_omic1=utils::read.csv(pathway_omic1_expr)
    #pheno0=utils::read.csv(pheno,sep="\t")
    #pathway_genes=utils::read.csv(pathway.genes)
    
    rdata_omic2= pathway_omic2_expr
    rdata_omic1= pathway_omic1_expr
    pheno0=pheno
    pathway_genes=pathway.genes
    cat("Finding the common sets/pathways for Omic1 and Omic2...","\n")

    pathwayID.omic2=unique(unlist(strsplit(rdata_omic2[,2],split="path:hsa"))[seq(2,2*nrow(rdata_omic2),2)])
    pathwayID.omic1=unique(unlist(strsplit(rdata_omic1[,2],split="path:hsa"))[seq(2,2*nrow(rdata_omic1),2)])
    # Omic1 and Omic2 in the same pathway
    pathway_common=pathwayID.omic1[as.numeric(stats::na.omit(match(pathwayID.omic2,pathwayID.omic1)))]
    pathway_common=paste0("path:hsa",pathway_common,sep="")
    u_path=pathway_common
    # make the order of subjects same in the three files: Omic2, Omic1, phenotype
    id0=NULL
    id0=which(pheno0[,4]==0)
    pheno1=rbind(pheno0[id0,],pheno0[-id0,])
    cn.omic2=colnames(rdata_omic2)[-(1:5)]
    cn.omic1=colnames(rdata_omic1)[-(1:5)]
    id.g=NULL;for (i in 1:nrow(pheno1)) id.g=c(id.g,which(cn.omic2==pheno1[i,1]))
    id.p=NULL;for (i in 1:nrow(pheno1)) id.p=c(id.p,which(cn.omic1==pheno1[i,1]))
    rdata_omic2.11=rdata_omic1.11=NULL
    rdata_omic2.11=cbind(rdata_omic2[,1:5],rdata_omic2[,(id.g+5)])
    rdata_omic1.11=cbind(rdata_omic1[,1:5],rdata_omic1[,(id.p+5)])
    subjID=intersect(pheno1[,1],colnames(rdata_omic2.11)[-(1:5)])
    pheno2=NULL
    for(i in 1:length(subjID))
    pheno2=rbind(pheno2,pheno1[which(pheno1[,1]==subjID[i]),])
    
    
    # ioSearch test statistic calculation
    cat("Starting calulation of ioSearch test statistic for set...","\n")
    test.statistic=matrix(0,nrow=length(u_path),ncol=4)
    for (i in 1:length(u_path))
    {
        tmp.f=NULL;tmp.f=as.vector(func_iosearch(i,rdata_omic1.11,rdata_omic2.11,phen=pheno2,sets=u_path,K, ntopOmics2, ntopOmics1))
        test.statistic[i,]=tmp.f
        rm(tmp.f)
    }
    cat("Test statistic calculation completed","\n")
    n=dim(pheno2)[1]
    pval.1=2*stats::pt(abs(test.statistic[,1]), df=n-1,lower.tail = F)
    pval.2=stats::pchisq(test.statistic[,2], df=test.statistic[,3],lower.tail = F)
    pval.3=2*stats::pt(abs(test.statistic[,4]), df=n-1,lower.tail = F)
    p.matrix=cbind(pval.1,pval.2,pval.3)
    H01.rejected=length(which(p.matrix[,1]<0.05))/nrow(p.matrix)
    H02.rejected=length(which(p.matrix[,2]<0.05))/nrow(p.matrix)
    H03.rejected=length(which(p.matrix[,3]<0.05))/nrow(p.matrix)
    pmax=apply(p.matrix,1,max)
    pmax.rejected=length(which(pmax<0.05))/length(pmax)
    PvalMat <- as.matrix(p.matrix)
    if(length(c(which(PvalMat[,1]==0),which(PvalMat[,2]==0),which(PvalMat[,3]==0),which(PvalMat[,1]==1),which(PvalMat[,2]==1),which(PvalMat[,3]==1)))>0)  PvalMat<-PvalMat[-c(which(PvalMat[,1]==0),which(PvalMat[,2]==0),which(PvalMat[,3]==0),which(PvalMat[,1]==1),which(PvalMat[,2]==1),which(PvalMat[,3]==1)),]

     Q <- ncol(PvalMat)
     AtLeast <- Q
     Tmp <- qch::GetHinfo(Q,AtLeast)
     Hconfig <- Tmp$Hconfig
     Hconfig.H1 <- Tmp$Hconfig.H1
     rm(Tmp)
     ResMMP <- qch::qch.fit(PvalMat,Hconfig)
     Alpha=0.05
     List <- qch::qch.test(ResMMP$posterior,Hconfig.H1,Alpha)
     I=which(List$Rejection==1)
     sig.path=u_path[I]
     Sig.path=NULL
     for(i in 1:length(sig.path))
     Sig.path=rbind(Sig.path,pathway_genes[which(pathway_genes[,1]==sig.path[i])[1],])
     Sig.path=Sig.path[,c(1,5)]
     #utils::write.csv(Sig.path,"signif.iosearch_pathways.csv",row.names=F)
     cat("Finding list of significant Omics...","\n")
     omic1.list=omic2.list=NULL
     for (i in I)
     {
         temp=temp.1=temp.2=NULL
         temp=f(i,mat1=rdata_omic1.11,mat2=rdata_omic2.11,phen=pheno2,sets=u_path,K, ntopOmics2, ntopOmics1)
         temp.2=temp$A[,1]
         temp.1=temp$B[,1]
         omic2.list=c(omic2.list,temp.2)
         omic1.list=c(omic1.list,temp.1)
     }
     result=list("Sig.path"=Sig.path,"signif.omic1"=unique(omic1.list),"signif.omic2"=unique(omic2.list))
     #files=dir()
     #unlink(files, recursive=TRUE,force=TRUE)
     return(result)
     #utils::write.csv(unique(omic2.list),"signif.iosearch_omic2.csv",row.names=FALSE,col.names=FALSE,quote=FALSE)
     #utils::write.csv(unique(omic1.list),"signif.iosearch_omic1.csv",row.names=FALSE,col.names=FALSE,quote=FALSE)
     cat("Done","\n")
}

