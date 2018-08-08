#' Run the scClustViz Shiny app
#'
#' Performs differential expression testing between clusters for all cluster solutions in
#' order to assess the biological relevance of each cluster solution. Differential
#' expression testing is done using the Wilcoxon rank-sum test implemented in the base R
#' \code{stats} package. For details about what is being compared in the tests, see the
#' "Value" section.
#'
#' @param filePath A character vector giving the relative filepath to an RData file
#'   containing two objects. One must be the list outputted by one of the importData
#'   functions (either \code{\link{readFromSeurat}} or \code{\link{readFromManual}})
#'   containing the data for viewing in the app. The other must be the list outputted by
#'   the \code{\link{clusterWiseDEtest}} function containing differential gene expression
#'   results for viewing in the app. As long as none of the name of the list elements have
#'   been changed, the objects can be named anything you'd like.
#'
#' @param dataTitle Optional. A character vector describing the data. This will be
#'   displayed at the top of the app. If missing it will be inferred from \code{filePath}.
#'
#' @param annotationDB Optional. An \code{AnnotationDb}-derived object containing
#'   identifier mappings for the species in the data. Used to search for genes if rownames
#'   are not converted to gene symbols, and for searching by and displaying gene names in
#'   the app. See examples.
#'
#' @param cellMarkers Optional. If you have canonical marker genes for expected cell
#'   types, list them here (see example code below). The Shiny app will attempt to label
#'   clusters in the tSNE projection by highest median gene expression. Optional.
#'
#' @param exponent The log base of your normalized input data. Seurat normalization uses
#'   the natural log (set this to exp(1)), while other normalization methods generally use
#'   log2 (set this to 2). This is used if you use the function for testing differential
#'   gene expression between custom sets, and should be set to the same parameters as in
#'   \code{clusterWiseDEtest}.
#'
#' @param pseudocount The pseudocount added to all log-normalized values in your input
#'   data. Most methods use a pseudocount of 1 to eliminate log(0) errors. This is used if
#'   you use the function for testing differential gene expression between custom sets,
#'   and should be set to the same parameters as in \code{clusterWiseDEtest}.
#'
#' @param FDR The false discovery rate to use as a threshold for determining statistical
#'   significance of differential expression calculated by the Wilcoxon rank-sum test.
#'   This is used if you use the function for testing differential gene expression between
#'   custom sets, and should be set to the same parameters as in \code{clusterWiseDEtest}.
#'
#' @param threshType Filtering genes for use in differential expression testing can be
#'   done multiple ways. We use an expression ratio filter for comparing each cluster to
#'   the rest of the tissue as a whole, but find that difference in detection rates works
#'   better when comparing clusters to each other. You can set threshType to
#'   \code{"logGER"} to use a gene expression ratio for all gene filtering, or leave it as
#'   default (\code{"dDR"}) to use difference in detection rate as the thresholding method
#'   when comparing clusters to each other. This is used if you use the function for
#'   testing differential gene expression between custom sets, and should be set to the
#'   same parameters as in \code{clusterWiseDEtest}.
#'
#' @param dDRthresh Magnitude of detection rate difference of a gene between clusters to
#'   use as filter for determining which genes to test for differential expression between
#'   clusters. This is used if you use the function for testing differential gene
#'   expression between custom sets, and should be set to the same parameters as in
#'   \code{clusterWiseDEtest}.
#'
#' @param logGERthresh Magnitude of gene expression ratio for a gene between clusters to
#'   use as filter for determining which genes to test for differential expression between
#'   clusters. This is used if you use the function for testing differential gene
#'   expression between custom sets, and should be set to the same parameters as in
#'   \code{clusterWiseDEtest}.
#'
#' @return The function causes the scClustViz Shiny GUI app to open in a seperate window.
#'
#' @examples
#' \dontrun{
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object,
#'                                        convertGeneIDs=F)
#'  rm(your_seurat_object)
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#'
#'  DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz)
#'
#'  save(data_for_scClustViz,DE_for_scClustViz,
#'       file="for_scClustViz.RData")
#'  # Save these objects so you'll never have to run this slow function again!
#'
#'  runShiny(filePath="for_scClustViz.RData",annotationDB=org.Mm.eg.db)
#' }
#'
#' @seealso \code{\link{readFromSeurat}} or \code{\link{readFromManual}} for reading in
#'   data to generate the first input object for this function, and
#'   \code{\link{clusterWiseDEtest}} to do the differential expression testing to generate
#'   the second input object for this function.
#'
#' @export

runShiny <- function(filePath,dataTitle,annotationDB,cellMarkers=list(),
                     exponent=2,pseudocount=1,FDR=0.01,
                     threshType="dDR",dDRthresh=0.15,logGERthresh=1) {

  while(T) {
    if (exists(".lastFileCall")) {
      if (names(.lastFileCall) == filePath) {
        if (exists(.lastFileCall[[1]][1]) & exists(.lastFileCall[[1]][2])) {
          break
        } else {
          rm(.lastFileCall)
        }
      } else {
        rm(.lastFileCall)
      }
    } else {
      .lastFileCall <- list(load(filePath))
      names(.lastFileCall) <- filePath
      break
    }
  }
  temp_objNames <- sapply(.lastFileCall[[filePath]],function(X) names(get(X)),simplify=F)
  il <- get(names(which(sapply(temp_objNames,function(X) "nge" %in% X))))
  dl <- get(names(which(sapply(temp_objNames,function(X) "CGS" %in% X))))
  rm(temp_objNames)
  temp_dataPath <- strsplit(filePath,"/|\\\\")
  dataPath <- sub(temp_dataPath[[1]][length(temp_dataPath[[1]])],"",filePath)
  if (dataPath == "") { dataPath <- "./" }
  if (missing("dataTitle")) {
    dataTitle <- sub("\\.[^.]+$","",tail(temp_dataPath[[1]],1))
  }
  rm(temp_dataPath)
  
  
  
  if (length(cellMarkers) < 1) {
    cellMarkersS <- cellMarkersU <- list()
  } else {
    cellMarkersS <- apply(combn(seq_along(cellMarkers),2),2,
                          function(X) do.call(intersect,unname(cellMarkers[X])))
    try(names(cellMarkersS) <- apply(combn(seq_along(cellMarkers),2),2,
                                     function(X) paste(X,collapse="&")),silent=T)
    cellMarkersS <- cellMarkersS[sapply(cellMarkersS,length) > 0]
    cellMarkersU <- lapply(cellMarkers,function(X) X[!X %in% unlist(cellMarkersS)])
  }
  
  demoRegex <- switch(species,mouse="^Actb$",human="^ACTB$")
  
  
}

