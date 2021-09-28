#' MSigdb signature scoring for single-cells
#' @title scMsigdbScoring
#' @param object Seurat object
#' @param species Species name. Such as "human" or "mouse"
#' @param category String of MSigDB category to use, 
#'   default category='H' stands for 'hallmark gene sets'
#' @param genesets Vector of genesets under the specific category of 
#'   MSigDB, such as "APOPTOSIS". Default genesets=NULL, which means 
#'   all genesets under the specified category will be used 
#' @return Seurat object
#'
#' @author Wei Zhang
#' @export
#'
#' @examples
#' data("H3N2_small")
#' x <- scMsigdbScoring(object = H3N2_small, 
#'   category = "H",
#'   genesets = c("HALLMARK_INFLAMMATORY_RESPONSE")
#' )
#' head(x)
#'
scMsigdbScoring <- function(object = NULL, species = "human", 
                            category = NULL, genesets = NULL) {
  if (is.null(category)){
    category <- "H"
  }
    
  if (category=="H" & (species == "human" | species == "Homo sapiens")){
    MSIGDB <- MSIGDB_H
  }else{
    all_gene_sets <- msigdbr::msigdbr(species=species, category=category)
    all_gene_sets <- as.data.frame(x[, c("gs_cat", "gs_name", "gene_symbol")])
    all_gene_sets <- aggregate(. ~ gs_cat + gs_name, data=x, FUN="c")
    MSIGDB <- list()
    MSIGDB[[category]] <- list()
    for (i in seq_len(nrow(all_gene_sets))){
      m[[y$gs_cat[i]]][[y$gs_name[i]]] <- unlist(y$gene_symbol[i])
    }
  }

  if (length(names(MSIGDB)) == 0) {
    stop("MSigDB has no valid entries!")
  }

  if (is.null(geneSets)) {
    geneSets = names(MSIGDB)
  } else {
    geneSets <- intersect(geneSets, names(MSIGDB))
    if (length(geneSets) == 0) {
      stop("geneSets ", geneSets, "does not exist in ", MSigDB, "!")
    }
  }
  for (item in geneSets) {
    if (item %in% names(object[[]])) {
      object[[item]] <- NULL
    }
    features <- MSIGDB[[item]]
    #features <- features[grep("^$", features, invert = TRUE)]
    features <- list(Score = features)
    object.hallmark <- Seurat::AddModuleScore(object = object, 
      features = features, name = item, ctrl = 
        min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1))))
    hallmark.columns <- grep(pattern = item, x = colnames(x = object.hallmark[[]]), value = TRUE)
    hallmark.scores <- object.hallmark[[hallmark.columns]]
    rm(object.hallmark)
    colnames(x = hallmark.scores) <- c(item)
    object[[colnames(x = hallmark.scores)]] <- hallmark.scores
  }
  reductions <- intersect(c("pca", "tsne", "umap"), names(object))
  for (reduction in reductions) {
    meta_ids <- gsub("coord", toupper(reduction), c("coord_1", "coord_2"))
    coord <- Seurat::Embeddings(object = object, reduction = reduction)[, c(1, 2)]
    object <- Seurat::AddMetaData(object = object, metadata = coord, col.name = meta_ids)
  }
  return(object)
}

