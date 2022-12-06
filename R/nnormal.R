#'  Test for normality in Groups 1 and 2 separately after transfromation and removing outliers
#'
#' @param samplesheet Samplesheet matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param sampledata_folder Omics data folder downloaded from TCGA or using TCGAbiolinks (as shown in the example).
#' @param clinical_samplesheet Clinical data matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param data.type1 Omic1 data type name (as in TCGA).
#' @param sample.type Sample type name (Tumor type such as "Primary Tumor", "Solid Tumor" etc.).
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
#' @return A list.
#' @export nnormal
#'
nnormal=function(samplesheet,sampledata_folder,clinical_samplesheet,data.type1,sample.type,pattern, omic1.expr.column,samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.case.id.column,samplesheet.data.type.column,clin.samplesheet.case.id.column,clin.samplesheet.stage.column,clin.samplesheet.age.column)
{
    B=Omic1_rmNA(samplesheet,sampledata_folder,clinical_samplesheet,data.type1,sample.type, pattern, omic1.expr.column,samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.case.id.column,samplesheet.data.type.column, clin.samplesheet.case.id.column,clin.samplesheet.stage.column,clin.samplesheet.age.column)

    b=NULL
    B.1=B$B.1
    for (i in 1:nrow(B.1))
    {
        tmp.0=as.numeric(B.1[i,])
        if(max(tmp.0)<0) tmp.0=-tmp.0
        tmp=(tmp.0-min(tmp.0)+1)/(max(tmp.0)+1)
        y = unique(tmp)
        #hist(y,breaks = 12)
        result = MASS::boxcox(y~1, lambda = seq(-5,5,0.1),plotit=FALSE)
        mylambda = result$x[which.max(result$y)]
        if (mylambda!=0) {y2 = (y^mylambda-1)/mylambda} else {y2=log(y)}
        y3=setdiff(y2,graphics::boxplot(y2,plot=F)$out) # remove the outliers from the transformed values
        b=c(b,stats::shapiro.test(y3)$p.value)
    }
    b.1=stats::p.adjust(b,"bonferroni")
    non.normal.omic1=NULL
    non.normal.omic1=which(b.1<0.05)
    b=NULL
    B.0=B$B.0
    for (i in 1:nrow(B.0))
    {
        tmp.0=as.numeric(B.0[i,])
        if(max(tmp.0)<0) tmp.0=-tmp.0
        tmp=(tmp.0-min(tmp.0)+1)/(max(tmp.0)+1)
        y = unique(tmp)
        #hist(y,breaks = 12)
        result = MASS::boxcox(y~1, lambda = seq(-5,5,0.1),plotit=FALSE)
        mylambda = result$x[which.max(result$y)]
        #cat(mylambda,"\t")
        if (mylambda!=0) {y2 = (y^mylambda-1)/mylambda} else {y2=log(y)}
        y3=setdiff(y2,graphics::boxplot(y2,plot=F)$out)
        #hist(y2)
        b=c(b,stats::shapiro.test(y3)$p.value)
    }
    b.1=stats::p.adjust(b,"bonferroni")
    non.normal.omic1=c(non.normal.omic1,which(b.1<0.05))

    # remove the omic1 that are not normally distributed and normalize the raw data
    tmp=NULL
    tmp=sort(unique(non.normal.omic1))
    if (length(tmp)>0)
    {
        B.0=B.0[-tmp,]
        B.1=B.1[-tmp,]
    } else{
        B.0=B.0
        B.1=B.1
    }
    omic1.names=B$omic1.names
    if (length(tmp)>0) {omic1.names=omic1.names[-tmp]} else {omic1.names=omic1.names}


    bct.0=function(i)
        {
            tmp.0=as.numeric(B.0[i,])
            if(max(tmp.0)<0) tmp.0=-tmp.0
            y=(tmp.0-min(tmp.0)+1)/(max(tmp.0)+1)
            result = MASS::boxcox(y~1, lambda = seq(-5,5,0.1),plotit=FALSE)
            mylambda = result$x[which.max(result$y)]
            if (mylambda!=0) {y2 = (y^mylambda-1)/mylambda} else {y2=log(y)}
            return(y2)
        }

    bct.1=function(i)
        {
            tmp.0=as.numeric(B.1[i,])
            if(max(tmp.0)<0) tmp.0=-tmp.0
            y=(tmp.0-min(tmp.0)+1)/(max(tmp.0)+1)
            result = MASS::boxcox(y~1, lambda = seq(-5,5,0.1),plotit=FALSE)
            mylambda = result$x[which.max(result$y)]
            if (mylambda!=0) {y2 = (y^mylambda-1)/mylambda} else {y2=log(y)}
            return(y2)
        }
    B.00=B.11=NULL
    B.00=do.call("rbind",lapply(1:nrow(B.0),bct.0))
    B.11=do.call("rbind",lapply(1:nrow(B.1),bct.1))

    # Paste the omic1 names in row and Subj IDs in col
    rownames(B.00)=omic1.names
    rownames(B.11)=omic1.names
    colnames(B.00)=colnames(B.0)
    colnames(B.11)=colnames(B.1)
    return(list("B.00"=B.00,"B.11"=B.11))
}
