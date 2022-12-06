#' Collates pathway information with omic1 (that requires antibody map) data values for all subjects
#'
#' @param pathway.genes Read the information obtained from find_pathwayGenes().
#' @param antibody_map Read the antibody map information directly downloaded from TCGA.
#' @param omic1.expr_Group1 Read omic1 data matrix for all Group 1 subjects obtained from Omic1_rmNoOmic2().
#' @param omic1.expr_Group2 Read omic1 data matrix for all Group 2 subjects obtained from Omic1_rmNoOmic2().
#'
#' @return A matrix with pathway information and omic1 expression values for all subjects.
#' @export pathwayOmic1_expr
#'
pathwayOmic1_expr=function(pathway.genes, antibody_map,omic1.expr_Group1,omic1.expr_Group2)
{
    cat("Combining information from Pathway and Omic1...","\n")
    antibody.map=utils::read.table(url(antibody_map),sep="\t",header=T)
    antibody.map.agid=antibody.map[,1]
    #pathway_genes=utils::read.csv(pathway.genes)
    pathway_genes=pathway.genes
    #rdata_omic1.0=utils::read.table(omic1.expr_Group1,header=T)
    rdata_omic1.0=omic1.expr_Group1
    #rdata_omic1.1=utils::read.table(omic1.expr_Group2,header=T)
    rdata_omic1.1=omic1.expr_Group2
    
    rn=rownames(rdata_omic1.0)
    antibody.map.genes=NULL
    for(i in 1:length(rn))
    antibody.map.genes=c(antibody.map.genes,antibody.map[which(antibody.map.agid==rn[i]),3])

    rdata_omic1=cbind(rdata_omic1.0,rdata_omic1.1)
    rdata_omic1=as.data.frame(rdata_omic1)
    rdata_omic1=cbind(antibody.map.genes,rdata_omic1)

    # duplicate genes are removed (some AGIDs are mapped to more than one gene)
    rdata_omic1=rdata_omic1[-which(duplicated(rdata_omic1[,1])=="TRUE"),]

    rdata_omic1_kegg=merge(pathway_genes,rdata_omic1,by.x="GeneSymbol",by.y="antibody.map.genes")  # merge by Ensembl IDs (contains kegg information with gene expression values for the patients)
    #utils::write.csv(rdata_omic1_kegg,"pathway_omic1.expr.csv",row.names=F)
    return(as.data.frame(rdata_omic1_kegg))
    cat("Done","\n")

}

