#' Create omic2 data matrix for both disease groups by collecting information from all subject-specific omic1 files.
#'
#' @param samplesheet Samplesheet matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param sampledata_folder Omics data folder downloaded from TCGA or using TCGAbiolinks (as shown in the example).
#' @param pheno Phenotype file created from Omic1_rmNoOmic2().
#' @param subject.id.0 Group1 subject IDs file created from Omic1_rmNoOmic2().
#' @param subject.id.1 Group2 subject IDs file created from Omic1_rmNoOmic2().
#' @param data.type1 Omic1 data type name (as in TCGA).
#' @param data.type2 Omic2 data type name (as in TCGA).
#' @param sample.type Sample type name (Tumor type such as "Primary Tumor", "Solid Tumor" etc.).
#' @param samplesheet.data.category.column Column number containing the data category in the samplesheet file.
#' @param gene.type Gene type of interest (eg. "protein_coding"). To be found in omics datafile. Could be an array as well.
#' @param pattern String to extract omic2 data files.
#' @param gene.type.column Column number containing the gene type in the omic-specific file.
#' @param samplesheet.sample.type.column Column number containing the sample type in the samplesheet file.
#' @param samplesheet.case.id.column Column number containing the case IDs in the samplesheet file.
#' @param samplesheet.data.type.column Column number containing the data type in the samplesheet file.
#' @param omic2.expr.column Column number containing the omic2 expression values in the omic-specific datafile
#'
#' @return Individual matrices for omic2 expression data for each disease group.
#' @export omic2
#'
omic2=function(samplesheet,sampledata_folder,pheno,subject.id.0,subject.id.1,data.type1,data.type2,sample.type,gene.type,pattern, gene.type.column ,omic2.expr.column, samplesheet.data.category.column, samplesheet.sample.type.column, samplesheet.case.id.column,samplesheet.data.type.column)
{
    #status0=utils::read.csv(pheno,sep="\t")
    status0=pheno
    cat("Reading input data files...","\n")
    a=utils::read.csv(samplesheet,sep="\t")
    D=dir(sampledata_folder, full.names = TRUE, ignore.case = TRUE,recursive = TRUE, pattern)
    
    subid.all.0=subject.id.0
    subid.all.1=subject.id.1
    

    #subid.all.0=utils::read.csv(subject.id.0,sep="\t",header=F)
    #subid.all.1=utils::read.csv(subject.id.1,sep="\t",header=F)
    
    #subid.all.0=subid.all.0[,1]
    #subid.all.1=subid.all.1[,1]
    
    a1=a[which(a[,samplesheet.data.category.column]==data.type2),]
    a1=a1[which(a1[,samplesheet.sample.type.column]==sample.type),]
    a1=a1[!duplicated(a1[,samplesheet.case.id.column]),]
    a1[,samplesheet.case.id.column]=gsub("-", ".", a1[,samplesheet.case.id.column])
    
    tmp.0=NULL
    for (i in 1:length(subid.all.0)) tmp.0=c(tmp.0,which(a1[,samplesheet.case.id.column]==subid.all.0[i]))
    tmp.1=NULL
    for (i in 1:length(subid.all.1)) tmp.1=c(tmp.1,which(a1[,samplesheet.case.id.column]==subid.all.1[i]))

    if (length(tmp.0)>0) {a2.0=a1[tmp.0,]} else {cat("Not enough Group1 omic1 samples","\n")}
    if (length(tmp.1)>0) {a2.1=a1[tmp.1,]} else {cat("Not enough Group2 omic1 samples","\n")}
    
   
    
    tmp.count=length(unlist(strsplit(D[1],split="/")))-1
    
    k=NULL
    for(i in 1:length(D))
    {
        tmp.file=unlist(strsplit(D[i],split="/"))[tmp.count]
        if (length(which(a2.1[,1]==tmp.file))>0) k=c(k,i)
        if (length(which(a2.0[,1]==tmp.file))>0) k=c(k,i)
    }


    # Retain omic2 files with clin and omic1 data
    subid.names=matrix(0,ncol=2)
    subid.names=subid.names[-1,]
    for (i in k)
    {
        tmp=tmp0=tmp1=NULL
        tmp=unlist(strsplit(D[i],split="/"))[tmp.count+1]
        tmp0=which(a2.0[,2]==tmp)
        tmp1=which(a2.1[,2]==tmp)
        if (length(tmp0)>0) subid.names=rbind(subid.names,c(a2.0[tmp0,samplesheet.case.id.column],i))
        if (length(tmp1)>0) subid.names=rbind(subid.names,c(a2.1[tmp1,samplesheet.case.id.column],i))
    }
    
    X=utils::read.table(D[1],sep="\t",header=T)
    X=X[which(X[,gene.type.column]==gene.type),]
    nr.omic2=nrow(X)
    cat("Collating data on omic2 for all individuals...","\n")
    cat("It might take a few minutes due to the large size of the files...","\n")


    X=NULL
    A=matrix(NA,nrow=nr.omic2)
    A=A[,-1]
    pb = utils::txtProgressBar(min = 1, max = length(k), initial = 1)
    stepi = 1
    
    
    for (i in k)
    {
        Sys.sleep(1)
        utils::setTxtProgressBar(pb,stepi)
        stepi = stepi + 1
        fname=paste0(c("e",i),collapse="")
        X=assign(fname,utils::read.table(D[i],sep="\t",header=T))
        X=X[which(X[,gene.type.column]==gene.type),]
        if (nrow(X)==nr.omic2) {A=cbind(A,X[,omic2.expr.column])} else
        A=cbind(A,rep(-99,nr.omic2))
    }
    
    
    close(pb)


    colnames(A)=subid.names[,1]
    X=utils::read.table(D[1],sep="\t",header=T)
    X=X[which(X[,gene.type.column]==gene.type),]
    rownames(A)=X[,2]

    subid.0=subid.1=NULL
    subid.0=status0[which(status0[,4]==0),1]
    subid.1=status0[which(status0[,4]==1),1]

    id0=id1=NULL
    for (i in 1:length(subid.0)) id0=c(id0,which(colnames(A)==subid.0[i]))
    for (i in 1:length(subid.1)) id1=c(id1,which(colnames(A)==subid.1[i]))
    
    BB.0=BB.1=NULL
    if (length(id0)>0) {BB.0=A[,id0]} else {cat("Not enough Group1 omic2 samples","\n")}
    if (length(id1)>0) {BB.1=A[,id1]} else {cat("Not enough Group2 omic2 samples","\n")}
    BB.00=BB.11=NULL
    BB.00=cbind(X[,2],BB.0)
    BB.11=cbind(X[,2],BB.1)

    rownames(BB.00)=c()
    rownames(BB.11)=c()
    
    cat("Creating omic2 data files to be used as ioSearch function inputs...","\n")

    BB.00=BB.00[-which(duplicated(BB.00[,1])=="TRUE"),]
    BB.11=BB.11[-which(duplicated(BB.11[,1])=="TRUE"),]

    #utils::write.table(BB.00,"omic2_Group1.txt",quote=F,row.names=F,col.names=T,sep="\t")
    #utils::write.table(BB.11,"omic2_Group2.txt",quote=F,row.names=F,col.names=T,sep="\t")

    #a0=utils::read.table("omic2_Group1.txt",header=T,sep="\t")
    #a1=utils::read.table("omic2_Group2.txt",header=T,sep="\t")

    #a0=a0[-which(duplicated(a0[,1])=="TRUE"),]
    #a1=a1[-which(duplicated(a1[,1])=="TRUE"),]

    #utils::write.table(a0,"omic2_Group1.txt",quote=F,row.names=T,col.names=T,sep="\t")
    #utils::write.table(a1,"omic2_Group2.txt",quote=F,row.names=T,col.names=T,sep="\t")

    
    return(list("omic2_Group1"=as.data.frame(BB.00), "omic2_Group2"= as.data.frame(BB.11)))

    
    #cat("Data files created","\n")

    
}
