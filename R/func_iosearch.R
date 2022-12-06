#' Runs iosearch algorithm to calculate test statistic for each set/pathway
#'
#' @param x a numeric value from the pathway list
#' @param mat1 Omic1 matrix with pathway information.
#' @param mat2 Omic2 matrix with pathway information.
#' @param phen Phenotype matrix.
#' @param sets Array of pathways.
#' @param K Default is 1. Tuning parameter specific to elasticnet package.
#' @param ntopOmics2 Number of top omics1 variables to be selected from each pathway (user-defined).
#' @param ntopOmics1 Number of top omics1 variables to be selected from each pathway (user-defined).
#'
#' @return An array.
#' @export func_iosearch
#'

func_iosearch=function(x,mat1,mat2,phen,sets,K, ntopOmics2, ntopOmics1)
{
    cat(x,"\t")
    full_data=list("G"=t(omic2_module(x,mat2,sets)),"P"=t(omic1_module(x,mat1,sets)),"Y"=phen[,4])
    result=suppressWarnings(as.numeric(iosearch_statistic(full_data,K, ntopOmics2, ntopOmics1)))
    return(result)
}


