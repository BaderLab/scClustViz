library(scales)
library(RColorBrewer)
library(grDevices)
#setEPS()

mean.logX <- function(data,ex=2,pc=1,pc.out) { log(mean(ex^data - pc) + pc.out,base=ex) }
halfway <- 10

testData <- t(sapply(seq(0,.099,.001),function(prob) rbinom(1e3,100,prob)))
GER <- rowMeans(testData)[halfway] / rowMeans(testData)
testMeans <- rowMeans(testData)
logTestData <- log2(testData + 1)


pseudocount.use <- c(1,1e-99,1/nrow(testData))
testLogMeans <- sapply(pseudocount.use,function(X) apply(logTestData,MARGIN=1,FUN=mean.logX,pc.out=X))
logGER <- apply(testLogMeans,2,function(X) X[halfway] - X)

cairo_ps(file="Fig1a.eps",height=6,width=6,fallback_resolution=600)
layout(matrix(c(3,1,4,2),2),widths=c(6,1),heights=c(1,6))
par(mar=c(3,3,.1,.1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=range(log2(GER)[-1]),ylim=c(min(logGER),max(logGER[,-2])),
     xlab="Log2(Gene Expression Ratio)",ylab="logGER with pseudocount")
abline(0,1)
for (i in 1:ncol(logGER)) {
  points(log2(GER)[-1],logGER[-1,i],pch=c(15,20,17)[i],cex=1.5,
         col=alpha(brewer.pal(3,"Dark2"),c(.5,1,.5))[i])
}
legend("topleft",bty="n",pch=c(15,20,17),col=alpha(brewer.pal(3,"Dark2"),c(.5,1,.5)),
       pt.cex=1.5,cex=1.2,legend=c("1","1e-99","1/#cells"),title="Pseudocount")
par(mar=c(3,.1,.1,1))
plot(x=c(1,1),y=logGER[1,-2],ylim=c(min(logGER),max(logGER[,-2])),
     xaxt="n",yaxt="n",xlab=NA,cex=1.5,
     pch=c(15,17),col=alpha(brewer.pal(3,"Dark2"),.5)[c(1,3)])
axis(1,1,expression(infinity))
par(mar=c(.1,3,1,.1))
plot(x=1,y=logGER[1,2],type="n",xaxt="n",ylab=NA,pch=20,cex=1.5)
par(mar=c(.1,.1,1,1))
plot(x=1,y=logGER[1,2],xaxt="n",yaxt="n",cex=1.5,
     pch=20,col=alpha(brewer.pal(3,"Dark2"),1)[2])
dev.off()


cairo_ps(file="Fig1b.eps",height=3,width=6,fallback_resolution=600)
layout(rbind(1:2),widths=c(6,1))
par(mar=c(3,3,1,.1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=c(min(logGER),max(logGER[,-2])),ylim=c(.5,3.5),
     yaxt="n",ylab="Pseudocount",xlab="logGER with pseudocount")
axis(2,1:3,c("1","1e-99","1/#cells"))
for (i in 1:ncol(logGER)) {
  if (i == 2) { drop2 <- -1 } else { drop2 <- 1:nrow(logGER) }
  boxplot(logGER[drop2,i],col=alpha(brewer.pal(3,"Dark2")[i],.5),
          add=T,at=i,horizontal=T,xaxt="n")
}
par(mar=c(3,.1,1,1))
plot(x=logGER[1,2],y=2,ylim=c(.5,3.5),yaxt="n",xlab=NA)
dev.off()




#### TESTING ####

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

testData <- cbind(rep(0,1000),sapply(as.vector(c(1,2,5) %*% t(c(0.001,0.01,0.1,1,10))),rpois,n=1000))
colMeans(testData)
log2TD <- log2(testData + 1)

meanDF <- data.frame(True=colMeans(testData),
                        PsOne=apply(log2TD,2,mean.logX,ex=2,pc=1,pc.out=1),
                        PsCells=apply(log2TD,2,mean.logX,ex=2,pc=1,pc.out=1/10000),
                        PsSmall=apply(log2TD,2,mean.logX,ex=2,pc=1,pc.out=1e-99))

plot(PsOne~True,data=meanDF,pch=20,col=alpha("black",.5),xlab="True mean",ylab="Calculated Mean")
abline(0,1)
plot(PsOne~True,data=meanDF,xlim=range(meanDF[,1]),ylim=range(meanDF[,2:3]),
     pch=20,col=alpha("black",.5),xlab="True mean",ylab="Calculated Mean")
points(PsCells~True,data=meanDF,
       pch=20,col=alpha("red",0.5))


fcDF <- data.frame(
  meanA=apply(gtools::permutations(ncol(testData),2),1,function(X) mean(testData[,X[1]])),
  meanB=apply(gtools::permutations(ncol(testData),2),1,function(X) mean(testData[,X[2]])),
  meanOneA=apply(gtools::permutations(ncol(testData),2),1,
                 function(X) mean.logX(log2TD[,X[1]],2,1,pc.out=1)),
  meanOneB=apply(gtools::permutations(ncol(testData),2),1,
                 function(X) mean.logX(log2TD[,X[2]],2,1,pc.out=1)),
  meanCellsA=apply(gtools::permutations(ncol(testData),2),1,
                   function(X) mean.logX(log2TD[,X[1]],2,1,pc.out=1/nrow(log2TD))),
  meanCellsB=apply(gtools::permutations(ncol(testData),2),1,
                   function(X) mean.logX(log2TD[,X[2]],2,1,pc.out=1/nrow(log2TD))),
  meanSmallA=apply(gtools::permutations(ncol(testData),2),1,
                               function(X) mean.logX(log2TD[,X[1]],2,1,pc.out=1e-99)),
  meanSmallB=apply(gtools::permutations(ncol(testData),2),1,
                   function(X) mean.logX(log2TD[,X[2]],2,1,pc.out=1e-99))
)
fcDF$TrueDiff <- fcDF$meanA - fcDF$meanB
fcDF$TrueLogFC <- log2(fcDF$meanA) - log2(fcDF$meanB)
fcDF$fcOne <- fcDF$meanOneA - fcDF$meanOneB
fcDF$fcCells <- fcDF$meanCellsA - fcDF$meanCellsB
fcDF$fcSmall <- fcDF$meanSmallA - fcDF$meanSmallB

fcDF[fcDF$TrueLogFC > -7 & fcDF$TrueLogFC < -6,]

plot(fcDF$TrueLogFC,fcDF$fcOne,
     xlim=range(fcDF$TrueLogFC[abs(fcDF$TrueLogFC) < Inf]),
     ylim=range(fcDF[abs(fcDF$TrueLogFC) < Inf,c("fcOne","fcCells","fcSmall")]),
     pch=20,col=alpha("black",.5),xlab="True FC",ylab="Calculated FC")
abline(0,1)
points(x=fcDF$TrueLogFC,fcDF$fcSmall,
       pch=20,col=alpha("blue",.5))
points(x=fcDF$TrueLogFC,fcDF$fcCells,
       pch=20,col=alpha("red",.5))

# For MGE and MDGE, calc with pseudocount of 1 - check to see this impact on MGE~DR plot
# storage of ps cells means?
