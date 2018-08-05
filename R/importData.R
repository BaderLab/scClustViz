#' Read in data from a Seurat object automatically
#'
#' Loads the necessary data from a Seurat object for use in both the cluster-wise 
#' differential expression testing function, as well as in the Shiny app itself.
#'
#' @param inD A Seurat object containing slots as outlined in Details.
#'
#' @param convertGeneIDs A logical indicating whether to convert rownames of
#'   \code{inD@@data} to official gene symbols (suggested for ease of data
#'   interpretation).
#'
#' @param mart A biomaRt \code{mart} object for use during the gene ID conversion with 
#'   both mart and dataset specified. See example code below. (only required if
#'   \code{convertGeneIDs = TRUE}).
#'
#' @param geneRowNames The biomaRt descriptor for the current gene name IDs (rownames of
#'   @@data, ie. "ensembl_gene_id"). See \code{biomaRt::listAttributes(mart)} for biomaRt
#'   descriptors. (only required if \code{convertGeneIDs = TRUE}).
#'
#' @param speciesSymbol The biomaRt descriptor for the official gene symbol to be
#'   converted to (ie. "mgi_symbol" for mouse or "hgnc_symbol" for human). (only required
#'   if \code{convertGeneIDs = TRUE}).
#'   
#' @return The function returns a list containing input data necessary for both the 
#'   cluster-wise differential expression testing function and the Shiny app itself. 
#'   The list contains the following elements:
#'   \describe{
#'     \item{nge}{The normalized gene expression matrix.} 
#'     \item{md}{The metadata dataframe, not including cluster assignments.} 
#'     \item{cl}{The cluster assignment dataframe, containing cluster assignments for each
#'       resolution tested. The columns will be sorted in order of increasing resolution 
#'       (k, number of clusters).}
#'     \item{dr_clust}{The cell embeddings used in the clustering, from PCA.}
#'     \item{dr_viz}{The cell embeddings used for visualization in 2D, from tSNE.} 
#'   }
#'
#' @section Seurat object slots: The following slots are expected in the Seurat object. If
#'   you're using Seurat v1.x, the equivalent slots are expected (this code takes
#'   advantage of \code{UpdateSeuratObject} to find the relevant data in older Seurat
#'   objects.)
#'   \describe{
#'     \item{@@data}{Holds the normalized gene expression matrix.}
#'     \item{@@meta.data}{Holds the metadata, including cluster assignments. Cluster
#'       assignment columns of the metadata should be titled with their resolution 
#'       parameters, as is the default in Seurat (ie. "res.0.8"). If no cell cycle 
#'       prediction method was used (determined by searching the metadata headings for 
#'       "cycle|phase|G2M"), Seurat has a function to predict cell cycle phase from 
#'       expression of canonical marker genes. These genes are stored as HGNC symbols, 
#'       but if your data is mouse it will try case-insensitive matches to homologues (in 
#'       which case you will see a warning in AddModuleScore indicating that it attempted 
#'       to match case). As a consequence of this, these predictions will only run if the
#'       row names of the data are or have been converted to official gene symbols.}
#'     \item{@@dr$pca@@cell.embeddings}{Holds the results of the PCA run by Seurat. The
#'       cell embeddings are used for the silhouette plot in the Shiny app. If Seurat v2.x 
#'       or greater was used, only the PC dimensions used in clustering will be considered
#'       in silhouette calculations. If an alternative dimensionality reduction method was
#'       used prior to clustering, use \code{readFromManual} to manually specify the 
#'       desired cell embeddings.}
#'     \item{@@dr$tsne@@cell.embeddings}{Holds the results of the tSNE run by Seurat. The
#'       cell embeddings are used for cell visualizations in the Shiny app. If an 
#'       alternative 2D projection method was used, use \code{readFromManual} to manually 
#'       specify the desired cell embeddings.}
#'   }
#'
#' @examples 
#' \dontrun{
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object,
#'                                        convertGeneIDs=F)
#'  rm(your_seurat_object) 
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#'  
#'  
#'  # Assuming your data is from mouse and the rownames are currently ensembl IDs.
#'  my_mart <- biomaRt::useMart(mart="ensembl",
#'                              dataset="mmusculus_gene_ensembl")
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object,
#'                                        convertGeneIDs=T,
#'                                        mart=my_mart,
#'                                        geneRowNames="ensembl_gene_id",
#'                                        speciesSymbol="mgi_symbol")
#'  rm(your_seurat_object) 
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#'  
#'  
#'  # Assuming your data is from human and the rownames are currently NCBI gene IDs.
#'  my_mart <- biomaRt::useMart(mart="ensembl",dataset="hsapiens_gene_ensembl")
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object,
#'                                        convertGeneIDs=T,
#'                                        mart=my_mart,
#'                                        geneRowNames="entrezgene",
#'                                        speciesSymbol="hgnc_symbol")
#'  rm(your_seurat_object) 
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#' }
#'
#' @family importData functions
#' 
#' @seealso https://satijalab.org/seurat/ for more information on the Seurat package, and 
#'   \code{\link{readFromManual}} for loading data by manually passing the requisite data 
#'   objects.
#' 
#' @export

readFromSeurat <- function(inD,convertGeneIDs=FALSE,
                           mart=NULL,geneRowNames=NULL,speciesSymbol=NULL) {
  temp_env <- ls(".GlobalEnv")
  if (class(inD) != "seurat") {
    stop("inD must be a Seurat object, otherwise use 'readFromManual()' to read in data.")
  }
  inD <- Seurat::UpdateSeuratObject(inD) 
  # In case your Seurat object is from an older version of Seurat
  
  out <- list()

  if (convertGeneIDs) {
    if (!requireNamespace("biomaRt",quietly=T)) {
      stop("To convert gene IDs, the bioconductor biomaRt package must be installed.")
    } 
    if (class(mart) == "Mart") {
      stop("The argument 'mart' must be a 'Mart' object generated by biomaRt::useMart().")
    }
    if (is.null(geneRowNames)) {
      stop("geneRowNames must be specified.")
    }
    if (is.null(speciesSymbol)) {
      stop("speciesSymbol must be specified.")
    }
    e2g <- biomaRt::getBM(attributes=c(geneRowNames,speciesSymbol),
                          mart=mart,filters=geneRowNames,
                          values=rownames(inD@data))
    e2g <- e2g[e2g[,speciesSymbol] != "",] 
    # removing unmapped gene symbols from conversion table
    print(paste(sum(duplicated(e2g[,geneRowNames])),
                geneRowNames,"mapped to multiple",speciesSymbol))
    # Arbitrarily picking one mapping for the above duplicates,
    # since these generally map to predicted genes anyway.
    e2g <- e2g[!duplicated(e2g[,geneRowNames]),] 
    rownames(e2g) <- e2g[,geneRowNames]
    inD@data <- inD@data[e2g[,geneRowNames],] # removing unmapped genes from data
    print(paste(sum(duplicated(e2g[,speciesSymbol])),
                speciesSymbol,"mapped to multiple",geneRowNames))
    # Going to collapse these by summing UMI counts between duplicated rows.
    temp_r <- inD@data[e2g[,speciesSymbol] %in% 
                    e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],]
    inD@data <- inD@data[!e2g[,speciesSymbol] %in% 
                 e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],]
    # Removed duplicated rows from data, saved as separate object
    rownames(inD@data) <- e2g[rownames(inD@data),speciesSymbol] 
    # renamed rows in data as gene symbols
    temp_r <- t(sapply(e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],function(X) 
      colSums(temp_r[e2g[,geneRowNames][e2g[,speciesSymbol] == X],])))
    # Collapsed by summing each duplicated gene symbol's row
    inD@data <- rbind(inD@data,temp_r)  
    # added those data back to matrix
  }
  out[["nge"]] <- inD@data
  
  out[["md"]] <- inD@meta.data[,!grepl("res\\.[0-9]",colnames(inD@meta.data))]  
  # metadata for cells (dataframe of cells), not including cluster assignments.
  
  if (is.data.frame(inD@meta.data[,grepl("res\\.[0-9]",colnames(inD@meta.data))])) {
    cl <- data.frame(lapply(inD@meta.data[,grepl("res\\.[0-9]",colnames(inD@meta.data))],as.factor))
  } else {
    cl <- data.frame(inD@meta.data[,grepl("res\\.[0-9]",colnames(inD@meta.data))])
    colnames(cl) <- grep("res\\.[0-9]",colnames(inD@meta.data),value=T)
  }
  rownames(cl) <- rownames(inD@meta.data) 
  out[["cl"]] <- cl[order(sapply(cl,function(X) length(levels(X))))]
  # cluster assignments per clustering resolution (dataframe: cells x cluster labels as factors)
  
  if (length(inD@calc.params) == 0) {
    out[["dr_clust"]] <- inD@dr$pca@cell.embeddings
  } else {
    out[["dr_clust"]] <- inD@dr$pca@cell.embeddings[,inD@calc.params$RunTSNE$dims.use]  
  }
  # cell embeddings in low-dimensional space used for clustering distances
  
  out[["dr_viz"]] <- inD@dr$tsne@cell.embeddings  
  # cell embeddings in 2D space for visualization (usually tSNE) (matrix: cells x coordinates)
  return(out)
}


#' Read in data manually by passing the requested objects
#'
#' Creates the data object expected by the cluster-wise differential expression testing 
#' function and the Shiny app. The user must provide the input objects as arguments.
#'
#' @param ngs A matrix (can be a sparse matrix from the Matrix package) of normalized gene
#'   expression per cell.  Genes on the rows, with row names as gene names (suggest using
#'   official gene symbols, will convert if requested - see \code{convertGeneIDs})
#' 
#' @param md A dataframe containing cell metadata, not including cluster assignments.
#' 
#' @param cl A dataframe containing cell cluster assignments, where every column is the 
#'   result of a clustering run with different parameters. The columns will be sorted in 
#'   order of increasing resolution (k, number of clusters).
#' 
#' @param dr_clust A matrix of cell embeddings in the reduced-dimensional space used for 
#'   clustering (ie. PCA), with rows of cells (with rownames), and columns of dimensions. 
#'   This will be used to calculate euclidean distance between cells for the silhouette 
#'   plot, so it will be more relevant if the scale of each dimension is weighted by the 
#'   variance it explains.
#' 
#' @param dr_viz A matrix of cell embeddings in 2D space used for visualization (ie. 
#'   tSNE), with rows of cells (with rownames), and 2 columns of dimensions. 
#'
#' @param convertGeneIDs A logical indicating whether to convert rownames of
#'   \code{inD@@data} to official gene symbols (suggested for ease of data
#'   interpretation).
#'
#' @param mart A biomaRt \code{mart} object for use during the gene ID conversion with 
#'   both mart and dataset specified. See example code in \code{\link{readFromSeurat}}. 
#'   (only required if \code{convertGeneIDs = TRUE}).
#'
#' @param geneRowNames The biomaRt descriptor for the current gene name IDs (rownames of
#'   @@data, ie. "ensembl_gene_id"). See \code{biomaRt::listAttributes(mart)} for biomaRt
#'   descriptors. (only required if \code{convertGeneIDs = TRUE}).
#'
#' @param speciesSymbol The biomaRt descriptor for the official gene symbol to be
#'   converted to (ie. "mgi_symbol" for mouse or "hgnc_symbol" for human). (only required
#'   if \code{convertGeneIDs = TRUE}).
#'
#' @return The function returns a list containing input data necessary for both the 
#'   cluster-wise differential expression testing function and the Shiny app itself. 
#'   The list contains the following elements:
#'   \describe{
#'     \item{nge}{The normalized gene expression matrix.} 
#'     \item{md}{The metadata dataframe, not including cluster assignments.} 
#'     \item{cl}{The cluster assignment dataframe, containing cluster assignments for each
#'       resolution tested. The columns will be sorted in order of increasing resolution 
#'       (k, number of clusters).}
#'     \item{dr_clust}{The cell embeddings used in the clustering, from PCA.}
#'     \item{dr_viz}{The cell embeddings used for visualization in 2D, from tSNE.} 
#'   }
#' 
#' @family importData functions
#' 
#' @seealso \code{\link{readFromSeurat}} for loading data from a Seurat object, and 
#'   examples of how \code{biomaRt} is used to convert gene IDs to gene symbols when 
#'   \code{convertGeneIDs = TRUE}.
#' 
#' @export

readFromManual <- function(nge,md,cl,dr_clust,dr_viz,
                           convertGeneIDs=FALSE,
                           mart=NULL,geneRowNames=NULL,speciesSymbol=NULL) {
  out <- list()
  if (convertGeneIDs) {
    if (!requireNamespace("biomaRt",quietly=T)) {
      stop("To convert gene IDs, the bioconductor biomaRt package must be installed.")
    } 
    if (class(mart) == "Mart") {
      stop("The argument 'mart' must be a 'Mart' object generated by biomaRt::useMart().")
    }
    if (is.null(geneRowNames)) {
      stop("geneRowNames must be specified.")
    }
    if (is.null(speciesSymbol)) {
      stop("speciesSymbol must be specified.")
    }
    e2g <- biomaRt::getBM(attributes=c(geneRowNames,speciesSymbol),
                          mart=mart,filters=geneRowNames,
                          values=rownames(nge))
    e2g <- e2g[e2g[,speciesSymbol] != "",] 
    # removing unmapped gene symbols from conversion table
    print(paste(sum(duplicated(e2g[,geneRowNames])),
                geneRowNames,"mapped to multiple",speciesSymbol))
    # Arbitrarily picking one mapping for the above duplicates,
    # since these generally map to predicted genes anyway.
    e2g <- e2g[!duplicated(e2g[,geneRowNames]),] 
    rownames(e2g) <- e2g[,geneRowNames]
    nge <- nge[e2g[,geneRowNames],] # removing unmapped genes from data
    print(paste(sum(duplicated(e2g[,speciesSymbol])),
                speciesSymbol,"mapped to multiple",geneRowNames))
    # Going to collapse these by summing UMI counts between duplicated rows.
    temp_r <- nge[e2g[,speciesSymbol] %in% 
                         e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],]
    nge <- nge[!e2g[,speciesSymbol] %in% 
                           e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],]
    # Removed duplicated rows from data, saved as separate object
    rownames(nge) <- e2g[rownames(nge),speciesSymbol] 
    # renamed rows in data as gene symbols
    temp_r <- t(sapply(e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],function(X) 
      colSums(temp_r[e2g[,geneRowNames][e2g[,speciesSymbol] == X],])))
    # Collapsed by summing each duplicated gene symbol's row
    nge <- rbind(nge,temp_r)  
    # added those data back to matrix
  }
  out[["nge"]] <- nge
  out[["md"]] <- md
  if (!all(sapply(cl,is.factor))) {
    temp_cells <- rownames(cl)
    cl <- data.frame(lapply(cl,as.factor))
    rownames(cl) <- temp_cells
  }
  out[["cl"]] <- cl[order(sapply(cl,function(X) length(levels(X))))]
  out[["dr_clust"]] <- dr_clust
  out[["dr_viz"]] <- dr_viz
  return(out)
}
