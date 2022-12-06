#' Reads the clinical data file 
#'
#' @param clinical_samplesheet Clinical data matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param clin.samplesheet.case.id.column Column number containing the case IDs in the clinical file.
#'
#' @return A character vector.
#' @export clin_subid
#'
clin_subid=function(clinical_samplesheet,clin.samplesheet.case.id.column)
{
    clin=utils::read.csv(clinical_samplesheet,sep="\t")
    return(clin[,clin.samplesheet.case.id.column])
}
