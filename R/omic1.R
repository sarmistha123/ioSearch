#' Read the data and create a matrix of Omic1 (eg. protein) names by subjID
#'
#' @param samplesheet Samplesheet matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param sampledata_folder Omics data folder downloaded from TCGA or using TCGAbiolinks (as shown in the example).
#' @param data.type1 Omic1 data type name (as in TCGA).
#' @param sample.type Sample type name (Tumor type such as "Primary Tumor", "Solid Tumor" etc.).
#' @param pattern String to extract omic1 data files.
#' @param omic1.expr.column Column number containing the omic1 expression values in the omic-specific datafile
#' @param samplesheet.sample.type.column Column number containing the sample type in the samplesheet file.
#' @param samplesheet.data.category.column Column number containing the data category in the samplesheet file.
#' @param samplesheet.data.type.column Column number containing the data type in the samplesheet file.
#'
#' @return A list.
#' @export omic1
#'
omic1=function(samplesheet,sampledata_folder,data.type1, sample.type, pattern, omic1.expr.column,
samplesheet.sample.type.column,samplesheet.data.category.column,samplesheet.data.type.column)
{
    cat("Reading input data files...","\n")
    a=utils::read.csv(samplesheet,sep="\t")
    
    a1=a[which(a[,samplesheet.sample.type.column]==sample.type),]
    a.id=which(a1[,samplesheet.data.type.column]==data.type1)
    D=dir(sampledata_folder, full.names = TRUE, ignore.case = TRUE,recursive = TRUE,pattern )
    
    tmp.count=length(unlist(strsplit(D[1],split="/")))
    
    nr.omic1=nrow(utils::read.csv(D[1],sep="\t"))
    cat("Collating data on omic1 for all individuals...","\n")
    A=matrix(NA,nrow=nr.omic1)
    A=A[,-1]
    sid=NULL
    for (i in 1:(length(a.id)))
    {
        sid.tmp=paste(unlist(strsplit(unlist(strsplit(D[i],split="/"))[tmp.count],split="-"))[1:3],collapse=".")
        sid=c(sid,sid.tmp)
        fname=paste0(c("z",i),collapse="")
        X=assign(fname,utils::read.csv(D[i],sep="\t"))
        if (nrow(X)==nr.omic1) {A=cbind(A,X[,omic1.expr.column])} else
        A=cbind(A,rep(-99,nr.omic1))
    
    }
    #length(unique(sid))
    colnames(A)=sid
    return(list("A"=A,"sid"=sid))

}
