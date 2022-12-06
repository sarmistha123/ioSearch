#' Runs iosearch algorithm to find top omics variables for significant set(s)/pathway(s)
#'
#' @param x A character vector with one element.
#' @param mat1 Omic1 matrix with pathway information.
#' @param mat2 Omic2 matrix with pathway information.
#' @param phen Phenotype matrix.
#' @param sets Array of pathways.
#' @param K Default is 1. Tuning parameter specific to elasticnet package.
#' @param ntopOmics2 Number of top omics1 variables to be selected from each pathway (user-defined).
#' @param ntopOmics1 Number of top omics1 variables to be selected from each pathway (user-defined).
#'
#' @return A list.
#' @export f
#'
f=function(x,mat1,mat2,phen,sets,K, ntopOmics2, ntopOmics1)
{
    cat(x,"\t")
    full_data=list("G"=t(omic2_module(x,mat2,sets)),"P"=t(omic1_module(x,mat1,sets)),"Y"=phen[,4])
    
    X=as.data.frame(apply(full_data$G,2,as.numeric))
    M=as.data.frame(apply(full_data$P,2,as.numeric))
    
    y=full_data$Y
    if (ncol(X)>1)
    {
        full.X.scaled=scale(X, scale = FALSE)
        Sx=svd(t(full.X.scaled)%*%full.X.scaled)
        Ax=Sx$v
        Wx=full.X.scaled%*%Ax
        PCx=Wx[,1]
        # select top omics 2 using pca
        X_cormat=t(full.X.scaled)%*%full.X.scaled
        S=elasticnet::spca(X_cormat, K=K,type="Gram",sparse="varnum",trace=FALSE,para=min(ntopOmics2,ncol(X)))
        A=as.matrix(S$loadings)
        
        reordered_variable.X=NULL
        for (i in 1:K)
        {
            reordered_variable.X=c(reordered_variable.X,which(abs(A[,i])>0))
        }
        # In case a variable appears with non-zero loading in two or more PCs
        reordered_variable.X=unique(reordered_variable.X)
     
    } else {
        PCx=X
        full.X.scaled=scale(X, scale = FALSE)
        reordered_variable.X=1
    }

    if (ncol(M)>1)
    {
        full.M.scaled=scale(M, scale = FALSE)
        Sm1=svd(t(full.M.scaled)%*%full.M.scaled)
        Am1=Sm1$v
        Wm1=full.M.scaled%*%Am1
        PCm1=Wm1[,1]
        # select top omics 1 using pca
        M_cormat=t(full.M.scaled)%*%full.M.scaled
        Sm=elasticnet::spca(M_cormat, K=K,type="Gram",sparse="varnum",trace=FALSE,para=min(ntopOmics1,ncol(M)))
        Am=as.matrix(Sm$loadings)
        
        reordered_variable.M=NULL
        for (i in 1:K)
        {
            reordered_variable.M=c(reordered_variable.M,which(abs(Am[,i])>0))
        }
        # In case a variable appears with non-zero loading in two or more PCs
        reordered_variable.M=unique(reordered_variable.M)
       
    }else {
        PCm1=M
        full.M.scaled=scale(M, scale = FALSE)
        reordered_variable.M=1
    }

    fit.glm=suppressWarnings(stats::glm(y~as.matrix(PCm1)+as.matrix(PCx),family=stats::binomial))
    #summary(fit.glm)


    tmp.omic2=mat2[which(mat2[,2]==sets[x]),]
    tmp.omic1=mat1[which(mat1[,2]==sets[x]),]

    return(list("A"=tmp.omic2[reordered_variable.X,1:5], "B"=tmp.omic1[reordered_variable.M,1:5]))
}

