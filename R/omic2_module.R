#' Collect omic2 variables belonging to a set/pathway
#'
#' @param x A numeric value.
#' @param mat Omic2 matrix.
#' @param sets Name of pathway.
#'
#' @return A matrix.
#' @export omic2_module
#'
omic2_module=function(x,mat,sets)
{
   tmp=mat[which(mat[,2]==sets[x]),]
   return(tmp[,-(1:5)])
}


