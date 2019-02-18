context("Distance calculations")
library(scClustViz)

temp_counts <- Matrix::Matrix(c(rpois(30,1),rpois(20,3)),nrow=5,ncol=10,sparse=T)
temp_nge <- Matrix::Matrix(log2(temp_counts + 1),sparse=T)
colnames(temp_nge) <- as.character(seq_len(ncol(temp_nge)))
rownames(temp_nge) <- as.character(seq_len(nrow(temp_nge)))
temp_cl <- factor(c(1,0,1,1,"two","two",0,0,0,0))

temp_scoreDist <-apply(combn(levels(temp_cl),2),2,function(i) 
  p.adjust(
    apply(temp_nge,1,function(X) 
      suppressWarnings(unlist(wilcox.test(X[temp_cl == i[1]],X[temp_cl == i[2]])["p.value"]))
    ),method="fdr")
)
temp_scoreDist[is.na(temp_scoreDist)] <- 1
temp_scoreDist <- colSums((-log10(temp_scoreDist))^2)^.5

temp_NN <- c("0"=c("1","two")[which.min(temp_scoreDist[1:2])],
            "1"=c("0","two")[which.min(temp_scoreDist[c(1,3)])],
            "two"=c("0","1")[which.min(temp_scoreDist[2:3])])


test_that("DEdist calculates distance from DE scores",{
  expect_equivalent(
    as.vector(as.dist(
      scClustViz:::fx_calcDist_scoreDE(
        scClustViz:::fx_calcDEcombn(
          nge=temp_nge,cl=temp_cl,
          deMes=scClustViz:::fx_calcEScombn(
            cl=temp_cl,DRthresh=0,
            CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
          )))
      )),
    temp_scoreDist
  )
})

test_that("DEdistNN calculates nearest neighbours",{
  expect_equivalent(
    scClustViz::DEdistNN(
      scClustViz:::fx_calcDist_scoreDE(
        scClustViz:::fx_calcDEcombn(
          nge=temp_nge,cl=temp_cl,
          deMes=scClustViz:::fx_calcEScombn(
            cl=temp_cl,DRthresh=0,
            CGS=scClustViz:::fx_calcCGS(nge=temp_nge,cl=temp_cl,exponent=2,pseudocount=1)
          )))),
    temp_NN
  )
})
