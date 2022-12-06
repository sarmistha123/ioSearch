#' Collates pathway information with omic2 data values for all subjects
#'
#' @param pathway.genes Read the information obtained from find_pathwayGenes().
#' @param omic2.expr_Group1 Read omic1 data matrix for all Group 1 subjects obtained from omic2().
#' @param omic2.expr_Group2 Read omic1 data matrix for all Group 2 subjects obtained from omic2().
#'
#' @return A matrix with pathway information and omic2 expression values for all subjects.
#' @export pathwayOmic2_expr
#'
pathwayOmic2_expr=function(pathway.genes,omic2.expr_Group1,omic2.expr_Group2)
{
    cat("Combining information from Pathway and Omic2...","\n")
    #rdata_omic2.0=utils::read.table(omic2.expr_Group1,header=T)
    #rdata_omic2.1=utils::read.table(omic2.expr_Group2,header=T)
    rdata_omic2.0=omic2.expr_Group1
    rdata_omic2.1=omic2.expr_Group2
    rdata_omic2=cbind(rdata_omic2.0,rdata_omic2.1[,-1])
    rdata_omic2.names=rdata_omic2[,1]
    sub=colnames(rdata_omic2)[-1]

    cat("Might take a few minutes due to large file size...","\n")
    rdata_omic2.tmp=do.call("rbind",lapply(1:nrow(rdata_omic2),function(x) {
          #cat(x,"\t")
          tmp<-as.numeric(replace(rdata_omic2[x,-1],which(rdata_omic2[x,-1]==0),1))
          return(log(tmp))
          
          
      }))

    rownames(rdata_omic2.tmp)=rdata_omic2.names
    colnames(rdata_omic2.tmp)=sub
    rdata_omic2=rdata_omic2.tmp
    rdata_omic2=cbind(rownames(rdata_omic2),rdata_omic2)
    rdata_omic2=as.data.frame(rdata_omic2)

    # no genes are duplicated in the above matrix

    #pathway_genes=utils::read.csv(pathway.genes)
    pathway_genes=pathway.genes
    rdata_omic2_kegg=merge(pathway_genes,rdata_omic2,by.x="GeneSymbol",by.y="V1")  # merge by Ensembl IDs (contains kegg information with gene expression values for the patients)
    #utils::write.csv(rdata_omic2_kegg,"pathway_omic2.expr.csv",row.names=F)
    return(rdata_omic2_kegg)
    cat("Done","\n")

}
