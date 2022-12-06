#' Calculates ioSearch test statistic
#'
#' @param full_data A list of omics1 and omics2 and phenotype data matrices
#' @param K Default is 1. Tuning parameter specific to elasticnet package.
#' @param ntopOmics2 Number of top omics1 variables to be selected from each pathway (user-defined).
#' @param ntopOmics1 Number of top omics1 variables to be selected from each pathway (user-defined).
#'
#' @return An array.
#' @export iosearch_statistic
#'
iosearch_statistic=function(full_data,K, ntopOmics2, ntopOmics1)
{
   # read data

    X=as.data.frame(apply(full_data$G,2,as.numeric))
    M=as.data.frame(apply(full_data$P,2,as.numeric))
    
    y=full_data$Y
    n=length(y)
   ###############################################################################################################
   # fit GLM to estimate alpha_1 and alpha_2
   if (ncol(X)>1)
   {
       full.X.scaled=scale(X, scale = FALSE)
       Sx=svd(t(full.X.scaled)%*%full.X.scaled)
       Ax=Sx$v
       Wx=full.X.scaled%*%Ax
       PCx=Wx[,1]
       # select top omics2 using pca
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
       # select top omics1 using pca
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
   estimates1.coeff=c(as.numeric(fit.glm$coefficients))
   estimates1.se=as.numeric(summary(fit.glm)$coefficients[,2])
   test.stat1=as.numeric(estimates1.coeff[2]/estimates1.se[2])
   test.stat3=as.numeric(estimates1.coeff[3]/estimates1.se[3])
#############################################################################################
   # Estimate vec(B) with only top omics2 and omics1
        X=as.matrix(as.data.frame(apply(full_data$G,2,as.numeric))[,reordered_variable.X])
        M=as.matrix(as.data.frame(apply(full_data$P,2,as.numeric))[,reordered_variable.M])
        n_pc=length(reordered_variable.X)
        s=length(reordered_variable.M)
        X=scale(X, scale = FALSE)
        S=A=NULL
        S=svd(t(X)%*%X)
        A=S$v
        W=as.matrix((X%*%A)[,1:n_pc])  # PC of OTUs: The req_PC explains >90% information
        t.A=as.matrix((t(A))[1:n_pc,1:n_pc])
        sigma2.g_initial=1
        if (ncol(M)>1)  G=stats::cov(M)  # Variance Covariance matrix of Omics1 (without the scalar product term)
        if (ncol(M)==1)  G=stats::var(M)
        traceG=psych::tr(G)
        solve.Vi=solve.Vi.vector=solve.V=NULL
        Vi=sigma2.g_initial*G
        gamma_hat=matrix(0,nrow=s,ncol=n_pc)
        A.star=matrix(0,nrow=n_pc,ncol=n_pc)
        for(ii in 1:n_pc)
        {
            for(jj in 1:n_pc)
            A.star[ii,jj]=sum(W[,ii]*W[,jj])
        }
        A1=matrixcalc::svd.inverse(A.star)%*%t(W)
        t.A1=t(A1)
        M.prime=t(scale(M, scale = FALSE))
        gamma.hat_matrix=matrix(0,nrow=s,ncol=n_pc)
        for(ii in 1:n_pc)
        {
                gamma.hat_matrix[,ii]=do.call("rbind",lapply(1:s,function(jj){sum(M.prime[jj,]*t.A1[,ii])}))
        }
        t.gamma.hat_matrix=gamma.hat_matrix
        gamma.hat_matrix=t(t.gamma.hat_matrix)
        B.hat.prime=matrix(0,nrow=s,ncol=n_pc)
        for (j in 1:n_pc)
        {
            for (i in 1:s)
            B.hat.prime[i,j]=sum(t.gamma.hat_matrix[i,]*t.A[,j])
        }
        B.hat.prime=as.matrix(B.hat.prime)
        X=as.matrix(X[,1:n_pc])
        A=as.matrix(A[1:n_pc,1:n_pc])
        int.term0=function(i)
        {
            int.term0.1=NULL
            for (j in 1:s)
            int.term0.1=c(int.term0.1,sum(X[i,]*B.hat.prime[j,]))
            return(int.term0.1)
        }
        int.term=NULL
        int.term=do.call("c",lapply(1:n,int.term0))
        int.term2=matrixcalc::vec(scale(M, scale = FALSE))-int.term
        int.term3=sum(int.term2^2)
        sigma2.g_hat=int.term3/(n*traceG)
        Vi=sigma2.g_hat*G
        Svd.inv=matrixcalc::svd.inverse(A.star)
        var_beta_hat=kronecker(A%*% Svd.inv %*%t(A),Vi)
        fastsvd=corpcor::fast.svd(var_beta_hat)
        if(length(fastsvd$d)==1)
        {
            inv.fastsvd=t(fastsvd$v)%*%(1/rep(fastsvd$d,length(fastsvd$v)))%*%(fastsvd$u)
            test.stat2=sum(matrixcalc::vec(B.hat.prime)^2)*inv.fastsvd
        }
        if(length(fastsvd$d)>1)
        {
            inv.fastsvd=fastsvd$v%*%(diag(1/fastsvd$d))%*%t(fastsvd$u)
            test.stat2=as.matrix(t(matrixcalc::vec(B.hat.prime)))%*%inv.fastsvd%*%as.matrix(matrixcalc::vec(B.hat.prime))
        }
        test.stat=c(test.stat1,test.stat2,n_pc*s,test.stat3)
        return(test.stat)
}
