context("Cluster-wise gene stats")
library(scClustViz)

temp_counts <- Matrix::Matrix(c(rpois(30,1),rpois(20,3)),nrow=5,ncol=10,sparse=T)
temp_nge <- Matrix::Matrix(log2(temp_counts + 1),sparse=T)
colnames(temp_nge) <- as.character(seq_len(ncol(temp_nge)))
rownames(temp_nge) <- as.character(seq_len(nrow(temp_nge)))
temp_cl <- factor(c(1,0,1,1,"two","two",0,0,0,0))

temp_DR <- sapply(levels(temp_cl),function(X) 
  apply(temp_counts[,temp_cl == X,drop=F],1,function(Y) 
    sum(Y > 0)/length(Y)))

temp_MGE <- sapply(levels(temp_cl),function(X) 
  log2(Matrix::rowMeans(temp_counts[,temp_cl == X,drop=F] + 1/ncol(temp_nge))))

temp_MDGE <- sapply(levels(temp_cl),function(X) 
  apply(temp_counts[,temp_cl == X,drop=F],1,function(Y) {
    temp <- Y[Y > 0]
    if(length(temp) == 0) { 
      return(0) 
    } else {
      return(log2(mean(temp + 1/ncol(temp_nge))))
    }
  }))


test_that("CalcCGS calculates DR",{
  expect_equal(
    sapply(scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1),
           function(X) X$DR),
    temp_DR
  )
})

test_that("CalcCGS and meanLogX calculate MGE",{
  expect_equal(
    sapply(scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1),
           function(X) X$MGE),
    temp_MGE
  )
})

test_that("CalcCGS and meanLogX calculate MGE",{
  expect_equal(
    sapply(scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1),
           function(X) X$MDGE),
    temp_MDGE
  )
})
