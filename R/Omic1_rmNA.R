#' Removes all missing data (such as NA, missing, no phenotype data) from omic1
#'
#' @param samplesheet Samplesheet matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param sampledata_folder Omics data folder downloaded from TCGA or using TCGAbiolinks (as shown in the example).
#' @param clinical_samplesheet Clinical data matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param data.type1 Omic1 data type name (as in TCGA).
#' @param sample.type Sample type name (Tumor type such as "Primary Tumor", "Solid Tumor" etc.).
#' @param pattern String to extract omic1 data files.
#' @param omic1.expr.column Column number containing the omic1 expression values in the omic-specific datafile
#' @param samplesheet.sample.type.column Column number containing the sample type in the samplesheet file.
#' @param samplesheet.data.category.column A character vector with one element.
#' @param samplesheet.case.id.column Column number containing the case IDs in the samplesheet file.
#' @param samplesheet.data.type.column Column number containing the data type in the samplesheet file.
#' @param clin.samplesheet.case.id.column Column number containing the case IDs in the clinical file.
#' @param clin.samplesheet.stage.column Column number containing the tumor stage in the clinical file.
#' @param clin.samplesheet.age.column Column number containing the age values in the clinical file.
#'
#' @return A list.
#' @export Omic1_rmNA
#'
Omic1_rmNA=function(samplesheet,sampledata_folder,clinical_samplesheet,data.type1,sample.type, pattern, omic1.expr.column,samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.case.id.column,samplesheet.data.type.column,clin.samplesheet.case.id.column,clin.samplesheet.stage.column,clin.samplesheet.age.column){
    Omic1=omic1(samplesheet,sampledata_folder,data.type1,sample.type,pattern, omic1.expr.column,
    samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.data.type.column)
    B=Omic1$A
    count.NA=which(is.na(B[,1]==TRUE))
    if (length(count.NA)>0) {B=B[-count.NA,]} else {B=B}
    
    # Remove subjects with unknown/ambiguous subtype
    status=dstatus(samplesheet, clinical_samplesheet,data.type1,samplesheet.data.type.column,samplesheet.case.id.column,clin.samplesheet.case.id.column,
    clin.samplesheet.stage.column,clin.samplesheet.age.column)

    tmp0=as.numeric(stats::na.omit(match(as.character(status[which(status[,4]==0),1]),Omic1$sid)))
    tmp1=as.numeric(stats::na.omit(match(as.character(status[which(status[,4]==1),1]),Omic1$sid)))

    B.0=B[,tmp0]
    B.1=B[,tmp1]
    rm(tmp0);rm(tmp1)
    D=dir(sampledata_folder, full.names = TRUE, ignore.case = TRUE,recursive = TRUE,pattern )

    X=utils::read.csv(D[1],sep="\t")
    omic1.names=X[-count.NA,1]
    
    # Remove subjects with missing values
    count.missvalues.0=which(B.0[1,]==-99)
    count.missvalues.1=which(B.1[1,]==-99)

    if (length(count.missvalues.0)>0) {B.0=B.0[,-count.missvalues.0]} else {B.0=B.0}
    if (length(count.missvalues.1)>0) {B.1=B.1[,-count.missvalues.1]} else {B.1=B.1}
    return(list("B.0"=B.0,"B.1"=B.1,"omic1.names"=omic1.names))

}
