#' Read in data from a Seurat object automatically
#'
#' Loads the necessary data from a Seurat object for use in both the
#' cluster-wise differential expression testing function, as well as in the
#' Shiny app itself.
#'
#' @param inD A Seurat object containing slots as outlined in Details.
#'
#' @return The function returns a list containing input data necessary for both
#'   the cluster-wise differential expression testing function and the Shiny app
#'   itself. The list contains the following elements: 
#'   \describe{ 
#'     \item{nge}{The normalized gene expression matrix.} 
#'     \item{md}{The metadata dataframe, not including cluster assignments.} 
#'     \item{cl}{The cluster assignment dataframe, containing cluster 
#'       assignments for each resolution tested. The columns will be sorted in 
#'       order of increasing resolution (k, number of clusters). Any dashes in 
#'       cluster names will be replaced with underscores to prevent conflicts 
#'       with comparison naming conventions.}
#'     \item{dr_clust}{The cell embeddings used in the clustering, from PCA.}
#'     \item{dr_viz}{The cell embeddings used for visualization in 2D, from 
#'       tSNE.}
#'   }
#'
#' @section Seurat object slots: The following slots are expected in the Seurat
#'   object. If you're using Seurat v1.x, the equivalent slots are expected
#'   (this code takes advantage of \code{UpdateSeuratObject} to find the
#'   relevant data in older Seurat objects.) 
#'   \describe{ 
#'     \item{@@data}{Holds the normalized gene expression matrix.} 
#'     \item{@@meta.data}{Holds the metadata, including cluster assignments. 
#'       \emph{Cluster assignment columns of the metadata should be titled with 
#'       their resolution parameters, as is the default in Seurat (ie. 
#'       "res.0.8").} }
#'     \item{@@dr$pca@@cell.embeddings}{Holds the results of the PCA run by 
#'       Seurat. The cell embeddings are used for the silhouette plot in the 
#'       Shiny app. If Seurat v2.x or greater was used, only the PC dimensions 
#'       used in clustering will be considered in silhouette calculations. If 
#'       an alternative dimensionality reduction method was used prior to 
#'       clustering, use \code{readFromManual} to manually specify the desired 
#'       cell embeddings.} 
#'     \item{@@dr$tsne@@cell.embeddings}{Holds the results of the tSNE run by 
#'       Seurat. The cell embeddings are used for cell visualizations in the 
#'       Shiny app. If an alternative 2D projection method was used, use 
#'       \code{readFromManual} to manually specify the desired cell embeddings.} 
#'   }
#'
#' @examples
#' \dontrun{
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object)
#'  rm(your_seurat_object)
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#'
#'  DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz)
#'
#'  save(data_for_scClustViz,DE_for_scClustViz,
#'       file="for_scClustViz.RData")
#'  # Save these objects so you'll never have to run this slow function again!
#'
#'  runShiny(filePath="for_scClustViz.RData")
#' }
#'
#' @family importData functions
#'
#' @seealso https://satijalab.org/seurat/ for more information on the Seurat
#'   package, and \code{\link{readFromManual}} for loading data by manually
#'   passing the requisite data objects.
#'
#' @export

readFromSeurat <- function(inD) {
  if (class(inD) != "seurat") {
    stop("inD must be a Seurat object, otherwise use 'readFromManual()' to read in data.")
  }
  if (.hasSlot(object, "version")) {
    if (packageVersion("Seurat") < package_version("2.0.0")) {
      inD <- Seurat::UpdateSeuratObject(inD) 
    }
  } else {
    inD <- Seurat::UpdateSeuratObject(inD) 
  }
  
  out <- list()
  out$nge <- inD@data
  
  out$md <- inD@meta.data[,!grepl("res\\.[0-9]",colnames(inD@meta.data))]
  for (X in which(sapply(out$md,is.character))) {
    out$md[[X]] <- as.factor(out$md[[X]])
  }
  # metadata for cells (dataframe of cells), not including cluster assignments.
  
  if (is.data.frame(inD@meta.data[,grepl("res\\.[0-9]",colnames(inD@meta.data))])) {
    cl <- data.frame(lapply(inD@meta.data[,grepl("res\\.[0-9]",colnames(inD@meta.data))],as.factor))
  } else {
    cl <- data.frame(inD@meta.data[,grepl("res\\.[0-9]",colnames(inD@meta.data))],stringsAsFactors=T)
    colnames(cl) <- grep("res\\.[0-9]",colnames(inD@meta.data),value=T)
  }
  if (nrow(cl) < 1) {
    stop("Can't find any clustering results. Try readFromManual to manually specify cluster results.")
  }
  rownames(cl) <- rownames(inD@meta.data)
  for (l in names(cl)) {
    levels(cl[[l]]) <- gsub("-","_",levels(cl[[l]]))
  }
  out$cl <- cl[order(sapply(cl,function(X) length(levels(X))))]
  # cluster assignments per clustering resolution (dataframe: cells x cluster labels as factors)
  
  if (length(inD@calc.params) == 0) {
    out$dr_clust <- inD@dr$pca@cell.embeddings
  } else {
    out$dr_clust <- inD@dr$pca@cell.embeddings[,inD@calc.params$RunTSNE$dims.use]  
  }
  # cell embeddings in low-dimensional space used for clustering distances
  
  out$dr_viz <- inD@dr$tsne@cell.embeddings  
  # cell embeddings in 2D space for visualization (usually tSNE) (matrix: cells x coordinates)
  return(out)
}


#' Read in data manually by passing the requested objects
#'
#' Creates the data object expected by the cluster-wise differential expression
#' testing function and the Shiny app. The user must provide the input objects
#' as arguments.
#'
#' @param nge A matrix (can be a sparse matrix from the Matrix package) of
#'   normalized gene expression per cell.  Genes on the rows, with row names as
#'   gene names (suggest using official gene symbols, will convert if requested
#'   - see \code{convertGeneIDs})
#'
#' @param md A dataframe containing cell metadata, not including cluster
#'   assignments.
#'
#' @param cl A dataframe containing cell cluster assignments, where every column
#'   is the result of a clustering run with different parameters. The columns
#'   will be sorted in order of increasing resolution (k, number of clusters).
#'
#' @param dr_clust A matrix of cell embeddings in the reduced-dimensional space
#'   used for clustering (ie. PCA), with rows of cells (with rownames), and
#'   columns of dimensions. This will be used to calculate euclidean distance
#'   between cells for the silhouette plot, so it will be more relevant if the
#'   scale of each dimension is weighted by the variance it explains.
#'
#' @param dr_viz A nx2 matrix of cell embeddings in two-dimensional space used
#'   for visualization of cells, with rows of cells (with rownames), and 2
#'   columns of dimensions. This is typically a tSNE projection, but any 2D
#'   embedding of cells is accepted.
#'
#' @return The function returns a list containing input data necessary for both
#'   the cluster-wise differential expression testing function and the Shiny app
#'   itself. The list contains the following elements: 
#'   \describe{ 
#'     \item{nge}{The normalized gene expression matrix.} 
#'     \item{md}{The metadata dataframe, not including cluster assignments.} 
#'     \item{cl}{The cluster assignment dataframe, containing cluster 
#'       assignments for each resolution tested. The columns will be sorted in 
#'       order of increasing resolution (k, number of clusters).Any dashes in 
#'       cluster names will be replaced with underscores to prevent conflicts 
#'       with comparison naming conventions.}
#'     \item{dr_clust}{The cell embeddings used in the clustering.}
#'     \item{dr_viz}{The cell embeddings used for visualization in 2D.}
#'   }
#'
#' @examples
#' \dontrun{
#'  ### Reading in data from a SingleCellExperiment class ###
#'  clusterAssignments <- grepl("^Clust",colnames(colData(mySCE)))
#'  # A logical vector separating the cluster assignments from the rest of the
#'  # cell metadata in the colData slot. This is an example that you will have
#'  # to change to reflect your cluster assignment column names.
#'  data_for_scClustViz <- readFromManual(nge=logcounts(mySCE),
#'                                        md=colData(mySCE)[,!clusterAssignments],
#'                                        cl=colData(mySCE)[,clusterAssignments],
#'                                        dr_clust=reducedDim(mySCE,"PCA"),
#'                                        dr_viz=reductedDim(mySCE,"tSNE"))
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#'
#'  DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz)
#'
#'  save(data_for_scClustViz,DE_for_scClustViz,
#'       file="for_scClustViz.RData")
#'  # Save these objects so you'll never have to run this slow function again!
#'
#'  runShiny(filePath="for_scClustViz.RData")
#' }
#'
#' @family importData functions
#'
#' @seealso \code{\link{readFromSeurat}} for loading data from a Seurat object.
#'
#' @export

readFromManual <- function(nge,md,cl,dr_clust,dr_viz) {
  out <- list()
  out$nge <- nge
  out$md <- md
  for (X in which(sapply(out$md,is.character))) {
    out$md[[X]] <- as.factor(out$md[[X]])
  }
  if (!is.data.frame(cl)) {
    out$cl <- data.frame(cl,stringsAsFactors=T)
    names(out$cl) <- deparse(substitute(cl))
  } else {
    out$cl <- cl
    for (X in which(sapply(out$cl,is.character))) {
      out$cl[[X]] <- as.factor(out$cl[[X]])
    }
  }
  for (l in names(out$cl)) {
    levels(out$cl[[l]]) <- gsub("-","_",levels(out$cl[[l]]))
  }
  out$cl <- out$cl[order(sapply(out$cl,function(X) length(levels(X))))]
  out$dr_clust <- dr_clust
  out$dr_viz <- dr_viz
  return(out)
}
