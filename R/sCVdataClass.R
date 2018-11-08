sCVparams <- setClass(Class="sCVparams",
                      slots=c(exponent="numeric",
                              pseudocount="numeric",
                              FDRthresh="numeric",
                              threshType="character",
                              dDRthresh="numeric",
                              logGERthresh="numeric"))


#' An S4 class to store analysis data for the scClustViz app
#' 
#' @slot CGS Cluster-wise gene stats.
#' 
#' @import methods

sCVdata <- setClass(Class="sCVdata",
                    slots=c(Clusts="factor", #need to strip '-' from cluster names
                            CGS="list", #list of length(getClusters(sCVdata)) containing dataframes 
                            DEvsRest="list", #list of length(getClusters(sCVdata)) containing dataframes
                            DEcombn="list", #list of length(getClusters(sCVdata)) containing lists of length(getClusters(sCVdata))-1 containing dataframes 
                            DEmarker="list", #list of length(getClusters(sCVdata)) containing dataframes
                            DEdist="matrix", #symmetrical matrix of side length(getClusters(sCVdata))
                            DEneighb="list", #list of length(getClusters(sCVdata)) containing dataframes
                            params="sCVparams"))

# Slot getters ----------
setGeneric("getClusters",function(x) standardGeneric("getClusters"))
setMethod("getClusters","sCVdata",function(x) x@Clusts)

setGeneric("getClustGeneStats",function(x) standardGeneric("getClustGeneStats"))
setMethod("getClustGeneStats","sCVdata",function(x) x@CGS)

setGeneric("getDEvsRest",function(x) standardGeneric("getDEvsRest"))
setMethod("getDEvsRest","sCVdata",function(x) x@DEvsRest)

setGeneric("getDEcombn",function(x) standardGeneric("getDEcombn"))
setMethod("getDEcombn","sCVdata",function(x) x@DEcombn)

setGeneric("getDEmarker",function(x) standardGeneric("getDEmarker"))
setMethod("getDEmarker","sCVdata",function(x) x@DEmarker)

setGeneric("getDEdist",function(x) standardGeneric("getDEdist"))
setMethod("getDEdist","sCVdata",function(x) x@DEdist)

setGeneric("getDEneighb",function(x) standardGeneric("getDEneighb"))
setMethod("getDEneighb","sCVdata",function(x) x@DEneighb)

setGeneric("getParam",function(x,param) standardGeneric("getParam"))
setMethod("getParam","sCVdata",function(x,param) slot(x@params,param))




#setValidity()



# Generics for object access -----

setGeneric("getGeneExpr",function(x) standardGeneric("getGeneExpr"))
setMethod("getGeneExpr","seurat",function(x) Seurat::GetAssayData(x))
setMethod("getGeneExpr","SingleCellExperiment",function(x) SingleCellExperiment::logcounts(x))

setGeneric("getMetadata",function(x) standardGeneric("getMetadata"))
setMethod("getMetadata","seurat",function(x) slot(eb1S,"meta.data"))
setMethod("getMetadata","SingleCellExperiment",function(x) SingleCellExperiment::colData(x))
