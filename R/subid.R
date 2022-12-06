#' Finds the subject IDs
#'
#' @param samplesheet file downloaded from TCGA or prepared using TCGAbiolinks package
#' @param data.type1 Select one data type to avoid subject ID repeats
#' @param samplesheet.data.type.column Indicate the data.type column number
#' @param samplesheet.case.id.column Indicate the case ID column number
#'
#' @return A character vector.
#' @export subid

subid=function(samplesheet,data.type1,samplesheet.data.type.column,samplesheet.case.id.column)
{
    a=utils::read.csv(samplesheet,sep="\t")
    a.id=which(a[,samplesheet.data.type.column]==data.type1)
    a2=a[a.id,]
    return(unique(a2[,samplesheet.case.id.column]))
}
