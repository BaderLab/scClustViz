context("DE vs Rest")
library(scClustViz)

temp_counts <- Matrix::Matrix(c(rpois(30,1),rpois(20,3)),nrow=5,ncol=10,sparse=T)
temp_nge <- Matrix::Matrix(log2(temp_counts + 1),sparse=T)
colnames(temp_nge) <- as.character(seq_len(ncol(temp_nge)))
rownames(temp_nge) <- as.character(seq_len(nrow(temp_nge)))
temp_cl <- factor(c(1,0,1,1,"two","two",0,0,0,0))

temp_logGER <- sapply(levels(temp_cl),function(X) 
  log2(Matrix::rowMeans(temp_counts[,temp_cl == X,drop=F] + 1/ncol(temp_nge))) - 
    log2(Matrix::rowMeans(temp_counts[,temp_cl != X,drop=F] + 1/ncol(temp_nge))))

temp_W <- sapply(levels(temp_cl),function(i)
  apply(temp_nge,1,function(X) 
    suppressWarnings(unlist(wilcox.test(X[temp_cl == i],X[temp_cl != i])["statistic"]))
  )
)

temp_FDR <- sapply(levels(temp_cl),function(i)
  p.adjust(
    apply(temp_nge,1,function(X) 
      suppressWarnings(unlist(wilcox.test(X[temp_cl == i],X[temp_cl != i])["p.value"]))
    ),method="fdr")
)


test_that("CalcDEvsRest calculates logGER",{
  expect_equal(
    sapply(scClustViz:::fx_calcESvsRest(
      nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1,DRthresh=0,
      CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
    ),function(X) X$logGER),
    temp_logGER
  )
})

test_that("CalcDEvsRest calculates Wilcoxon test statistic",{
  expect_equivalent(
    sapply(scClustViz:::fx_calcDEvsRest(
      nge=temp_nge,cl=temp_cl,
      deTes=scClustViz:::fx_calcESvsRest(
        nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1,DRthresh=0,
        CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
      )),function(X) X$Wstat),
      temp_W
    )
})

test_that("CalcDEvsRest calculates FDR",{
  expect_equivalent(
    sapply(scClustViz:::fx_calcDEvsRest(
      nge=temp_nge,cl=temp_cl,
      deTes=scClustViz:::fx_calcESvsRest(
        nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1,DRthresh=0,
        CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
      )),function(X) X$FDR),
    temp_FDR
  )
})

