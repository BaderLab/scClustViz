context("pairwise DE")
library(scClustViz)

temp_counts <- Matrix::Matrix(c(rpois(30,1),rpois(20,3)),nrow=5,ncol=10,sparse=T)
temp_nge <- Matrix::Matrix(log2(temp_counts + 1),sparse=T)
colnames(temp_nge) <- as.character(seq_len(ncol(temp_nge)))
rownames(temp_nge) <- as.character(seq_len(nrow(temp_nge)))
temp_cl <- factor(c(1,0,1,1,"two","two",0,0,0,0))

temp_logGER <- apply(combn(levels(temp_cl),2),2,function(i) 
  log2(Matrix::rowMeans(temp_counts[,temp_cl == i[1],drop=F] + 1/ncol(temp_nge))) - 
    log2(Matrix::rowMeans(temp_counts[,temp_cl == i[2],drop=F] + 1/ncol(temp_nge))))

temp_dDR <- apply(combn(levels(temp_cl),2),2,function(i) 
  apply(temp_counts[,temp_cl == i[1],drop=F],1,function(X) sum(X > 0) / length(X)) -
    apply(temp_counts[,temp_cl == i[2],drop=F],1,function(X) sum(X > 0) / length(X)))
    


temp_W <- apply(combn(levels(temp_cl),2),2,function(i) 
  apply(temp_nge,1,function(X) 
    suppressWarnings(unlist(wilcox.test(X[temp_cl == i[1]],X[temp_cl == i[2]])["statistic"]))
  )
)

temp_FDR <- apply(combn(levels(temp_cl),2),2,function(i) 
  p.adjust(
    apply(temp_nge,1,function(X) 
      suppressWarnings(unlist(wilcox.test(X[temp_cl == i[1]],X[temp_cl == i[2]])["p.value"]))
    ),method="fdr")
)


test_that("CalcDEcombn calculates logGER",{
  expect_equivalent(
    sapply(scClustViz:::fx_calcEScombn(
      cl=temp_cl,DRthresh=0,
      CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
    ),function(X) X$logGER),
    temp_logGER
  )
})

test_that("CalcDEcombn calculates dDR",{
  expect_equivalent(
    sapply(scClustViz:::fx_calcEScombn(
      cl=temp_cl,DRthresh=0,
      CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
    ),function(X) X$dDR),
    temp_dDR
  )
})

test_that("CalcDEvsRest calculates Wilcoxon test statistic",{
  expect_equivalent(
    sapply(scClustViz:::fx_calcDEcombn(
      nge=temp_nge,cl=temp_cl,
      deMes=scClustViz:::fx_calcEScombn(
        cl=temp_cl,DRthresh=0,
        CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
      )),function(X) X$Wstat),
    temp_W
  )
})

test_that("CalcDEvsRest calculates FDR",{
  expect_equivalent(
    sapply(scClustViz:::fx_calcDEcombn(
      nge=temp_nge,cl=temp_cl,
      deMes=scClustViz:::fx_calcEScombn(
        cl=temp_cl,DRthresh=0,
        CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
      )),function(X) X$FDR),
    temp_FDR
  )
})

