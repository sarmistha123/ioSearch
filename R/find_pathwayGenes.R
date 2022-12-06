#' Extracting pathway information for each gene from KEGG database
#'
#'
#' @return A matrix with genes available from subject-specific files and the pathways they are involved in.
#' @export find_pathwayGenes
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom KEGGREST keggLink keggList
#' @importFrom org.Hs.eg.db org.Hs.eg.db

#' @examples
#' ### Extracting pathway information for each gene from KEGG database
#'
#' # library(AnnotationDbi)
#' # library(KEGGREST)
#' # library(org.Hs.eg.db)
#' # find_pathwayGenes()
#'
#' # pathway.genes=find_pathwayGenes()
#'
#' ## DO NOT RUN
find_pathwayGenes=function()
{
    #require(AnnotationDbi)
    #require(org.Hs.eg.db)
    #require(KEGGREST)
    
    #db <- available.packages(repos=BiocManager::repositories(), type = "win.binary")
    #tools::package_dependencies("AnnotationDbi", db, recursive = TRUE)
    #tools::package_dependencies("org.Hs.eg.db", db, recursive = TRUE)
    #tools::package_dependencies("KEGGREST", db, recursive = TRUE)
    #require(KEGGREST)
    #require(org.Hs.eg.db)
    #require(AnnotationDbi)
    #BiocInstaller::biocLite("AnnotationDbi")
    #BiocInstaller::biocLite("org.Hs.eg.db")
    #BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi","KEGGREST"))
    
    
    #requireNamespace("AnnotationDbi")
    #requireNamespace("KEGGREST")
    #requireNamespace("BiocManager")
    #requireNamespace("org.Hs.eg.db")
    
   
    cat("Extracting pathway information from KEGG database...","\n")
    `%>%` <- magrittr::`%>%`
    `mutate` <- dplyr::mutate
    eg <- NULL
    .<- NULL
    `mapIds` <- AnnotationDbi::mapIds


    hsa_path_eg  <- KEGGREST::keggLink("pathway", "hsa") %>%  dplyr::tibble(pathway = ., eg = sub("hsa:", "", names(.)))
#AnnotationDbi
    hsa_kegg_anno <- hsa_path_eg %>% mutate(symbol = mapIds(org.Hs.eg.db::org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),ensembl = mapIds(org.Hs.eg.db::org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID"))
    #hsa_kegg_anno

    hsa_pathways <- KEGGREST::keggList("pathway", "hsa") %>% dplyr::tibble(pathway = names(.), description = .)

    full_table=dplyr::left_join(hsa_kegg_anno, hsa_pathways)
    full_table[,2]=paste0("hsa",as.numeric(unlist(full_table[,2])),sep="")


    colnames(full_table)=paste(c("KEGG_PathwayID","KEGG_GeneID","GeneSymbol","EnsemblID","PathwayDescription"))
    #utils::write.csv(full_table,"pathway_genes.csv",quote=F,row.names=F)
    return(full_table)
    cat("Done","\n")
    

}
