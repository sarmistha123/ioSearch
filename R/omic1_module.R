#' Collect omic1 variables belonging to a set/pathway
#'
#' @param x A numeric value.
#' @param mat Omic1 matrix.
#' @param sets Name of pathway.
#'
#' @return A matrix.
#' @export omic1_module
#'
omic1_module=function(x,mat,sets)
{
   tmp=mat[which(mat[,2]==sets[x]),]
   return(tmp[,-(1:5)])
}


