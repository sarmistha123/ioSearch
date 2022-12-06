#' Create omic1 data matrix for both disease groups by collecting information from all subject-specific omic1 files.
#'
#' @param samplesheet Samplesheet matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param data.type1 Omic1 data type name (as in TCGA).
#' @param data.type2 Omic2 data type name (as in TCGA).
#' @param sample.type Sample type name (Tumor type such as "Primary Tumor", "Solid Tumor" etc.).
#' @param sampledata_folder Omics data folder downloaded from TCGA or using TCGAbiolinks (as shown in the example).
#' @param clinical_samplesheet Clinical data matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param pattern String to extract omic1 data files.
#' @param omic1.expr.column Column number containing the omic1 expression values in the omic-specific datafile
#' @param samplesheet.sample.type.column Column number containing the sample type in the samplesheet file.
#' @param samplesheet.data.category.column Column number containing the data category in the samplesheet file.
#' @param samplesheet.case.id.column Column number containing the case IDs in the samplesheet file.
#' @param samplesheet.data.type.column Column number containing the data type in the samplesheet file.
#' @param clin.samplesheet.case.id.column Column number containing the case IDs in the clinical file.
#' @param clin.samplesheet.stage.column Column number containing the tumor stage in the clinical file.
#' @param clin.samplesheet.age.column Column number containing the age values in the clinical file.
#'
#' @return Individual matrices for phenotype, subject IDs in each disease group, omic1 expression data for each disease group.
#' @export Omic1_rmNoOmic2
#'
Omic1_rmNoOmic2=function(samplesheet,data.type1, data.type2,sample.type,sampledata_folder,clinical_samplesheet,pattern, omic1.expr.column,samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.case.id.column,samplesheet.data.type.column,clin.samplesheet.case.id.column,clin.samplesheet.stage.column,clin.samplesheet.age.column)
{
    B=nnormal(samplesheet,sampledata_folder,clinical_samplesheet,data.type1,sample.type,pattern, omic1.expr.column,samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.case.id.column,samplesheet.data.type.column,clin.samplesheet.case.id.column,clin.samplesheet.stage.column,clin.samplesheet.age.column)
    B.00=B$B.00
    B.11=B$B.11
    
    a=utils::read.csv(samplesheet,sep="\t")
    a1=a[which(a[,samplesheet.data.category.column]==data.type2),]
    a1=a1[which(a1[,samplesheet.sample.type.column]==sample.type),]
    a1=a1[!duplicated(a1[,samplesheet.case.id.column]),]
    subids=unique(a1[,samplesheet.case.id.column])
    subids=gsub("-", ".", subids)

    subid.0=colnames(B.00)
    subid.1=colnames(B.11)
    subid.all.0=intersect(subids,subid.0)
    subid.all.1=intersect(subids,subid.1)


    id.0=id.1=NULL
    for (i in 1:length(subid.all.0)) id.0=c(id.0,which(subid.0==subid.all.0[i]))
    for (i in 1:length(subid.all.1)) id.1=c(id.1,which(subid.1==subid.all.1[i]))

    if (length(id.0)>0) {B.000=B.00[,id.0]} else {B.000=B.00}
    if (length(id.1)>0) {B.111=B.11[,id.1]} else {B.111=B.11}
    
    cat("Creating data files to be used as ioSearch function inputs...","\n")
    #utils::write.table(B.000,"omic1_Group1.txt",quote=F,row.names=T,col.names=T,sep="\t")
    #utils::write.table(B.111,"omic1_Group2.txt",quote=F,row.names=T,col.names=T,sep="\t")
    #utils::write.table(subid.all.0,"subid.all.0.txt",quote=F,row.names=F,col.names=F,sep="\t")
    #utils::write.table(subid.all.1,"subid.all.1.txt",quote=F,row.names=F,col.names=F,sep="\t")
    
    return(list("omic1_Group1"=B.000, "omic1_Group2"=B.111, "subid.all.0"=subid.all.0, "subid.all.1"=subid.all.1))
    
    #cat("Data files created","\n")
    
}
