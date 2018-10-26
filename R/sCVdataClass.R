#' An S4 class to store analysis data for the scClustViz app
#' 
#' @slot CGS Cluster-wise gene stats.

sCVdata <- setClass(Class="sCVdata",
                    slots=c(cl="data.frame",
                            CGS="data.frame",
                            DEtissue=""
                            # params:
                            exponent="numeric",
                            pseudocount="numeric",
                            FDRthresh="numeric",
                            threshType="character",
                            dDRthresh="numeric",
                            logGERthresh="numeric"))


setGeneric("getCGS",function(x) standardGeneric("getCGS"))
setMethod("getCGS","sCVdata",function(x) x@CGS)

setGeneric("getParam",function(x,param) standardGeneric("getParam"))
setMethod("getParam","sCVdata",function(x,param) slot(x,param))

setValidity()