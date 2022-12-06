#' Divides the disease status into binary group and creates a phenoytype matrix
#'
#' @param samplesheet Samplesheet matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param clinical_samplesheet Clinical data matrix downloaded from TCGA or constructed with TCGAbiolinks (as shown in the example).
#' @param data.type1 Omic1 data type name (as in TCGA).
#' @param samplesheet.data.type.column Column number containing the data type in the samplesheet file.
#' @param samplesheet.case.id.column Column number containing the case IDs in the samplesheet file.
#' @param clin.samplesheet.case.id.column Column number containing the case IDs in the clinical file.
#' @param clin.samplesheet.stage.column Column number containing the tumor stage in the clinical file.
#' @param clin.samplesheet.age.column Column number containing the age values in the clinical file.
#'
#' @return A matrix.
#' @export dstatus
#'
dstatus=function(samplesheet, clinical_samplesheet,data.type1,samplesheet.data.type.column,samplesheet.case.id.column,clin.samplesheet.case.id.column,
clin.samplesheet.stage.column,clin.samplesheet.age.column)
{
    clin=utils::read.csv(clinical_samplesheet,sep="\t")
    status=clin[as.numeric(stats::na.omit(match(subid(samplesheet,data.type1,samplesheet.data.type.column,samplesheet.case.id.column),clin_subid(clinical_samplesheet,clin.samplesheet.case.id.column)))),c(clin.samplesheet.case.id.column,clin.samplesheet.stage.column,clin.samplesheet.age.column)]
    
    id=id.1=id.2=id.3=id.4=id.5=id.6=NULL
    id=c(which(status[,2]==""),which(status[,2]=="'--"),which(status[,2]=="Stage X"))
    if(length(id)==0) {status=status} else
    {status=status[-id,]}  # remove Stages blank and X

status=cbind(status,rep(1,nrow(status)))
id.1=which(status[,2]=="Stage IB");  if (length(id.1)>0) {status[id.1,4]<-0} else {status=status}
id.2=which(status[,2]=="Stage IA");  if (length(id.2)>0) {status[id.2,4]<-0} else {status=status}
id.3=which(status[,2]=="Stage I");   if (length(id.3)>0) {status[id.3,4]<-0} else {status=status}
id.4=which(status[,2]=="Stage IIB"); if (length(id.4)>0) {status[id.4,4]<-0} else {status=status}
id.5=which(status[,2]=="Stage IIA"); if (length(id.5)>0) {status[id.5,4]<-0} else {status=status}
id.6=which(status[,2]=="Stage II");  if (length(id.6)>0) {status[id.6,4]<-0} else {status=status}
status[,1]=gsub("-", ".", status[,1])
colnames(status)[4]<-"Group"
#utils::write.table(status,"status.txt",quote=F,row.names=T,col.names=T,sep="\t")
return(status)

}
