#' Computes the mean of log-scaled values
#' 
#' Computes the arithmetic mean in linear space of log-scaled values, and returns the result
#' in the same log scale.
#' 
#' Generally a pseudocount of 1 is added to log-scaled values to prevent +/-Inf results. 
#' However, adding a pseudocount of 1 to the log-scaled mean prior to gene expression ratio 
#' calculations in log space skews the result quite dramatically, so instead we add a small 
#' pseudocount to avoid +/- inf results when means are zero, without the same skewing. Adding 
#' a very small (ie 1e-99) number means that means of zero get set to a large negative log-mean, 
#' when it might be more appropriate to have those values fall closer to the smallest non-zero 
#' log-mean. By using a pseudocount of 1 / number of samples in the experiment, we ensure that 
#' log(zero) is smaller than any non-zero log-mean, while still being in the same ballpark.
#' 
#' @param data A numeric vector in log scale for which the arithmetic mean in linear scale 
#' is to be calculated and returned in the same log scale. 
#' 
#' @param ncell The number of samples (cells assuming an scRNAseq experiment) in the data. 
#' This is not necessarily the same as the length of \code{data}. Since its inverse will 
#' be used as the pseudocount value for all calculations, it should be consistent for all 
#' calls of \code{mean.logX}.
#' 
#' @param ex The log base of the data. Traditionally gene expression data is represented 
#' in base 2, although some methods (ie. Seurat's normalization scheme) use the natural log. 
#' (default=2)
#' 
#' @param pc The pseudocount used when converting the data to log-scale initially. (default=1)
#' 
#' @export

meanLogX <- function(data,ncell,ex=2,pc=1) { 
  log(mean(ex^data - pc) + 1/ncell,base=ex)
}


#' Generates a rainbow colour-scale
#' 
#' Generates a slightly desaturated rainbow colour-scale that's a little easier on the
#' eyes than the R default.
#' 
#' @param n Number of colours to include.

rainbow2 <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 60, c = 100)[1:n]
}
