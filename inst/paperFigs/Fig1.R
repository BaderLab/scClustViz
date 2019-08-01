library(scClustViz)
#setEPS()

mean.logX <- function(data,ex=2,pc=1,pc.out) { log(mean(ex^data - pc) + pc.out,base=ex) }

temp_lambda <- as.vector(c(1,2,5) %*% t(c(0.001,0.01,0.1,1,10)))
temp_n <- runif(length(temp_lambda),100,1000)

testData <- lapply(seq_along(temp_n),function(i) rpois(temp_n[i],temp_lambda[i]))
tempMeans <- sapply(testData,mean)
if (all(tempMeans > 0)) {
  testData[[length(testData) + 1]] <- rep(0,100)
}
log2TD <- lapply(testData,function(X) log2(X + 1))

permTD <- gtools::permutations(length(testData),2)

fcDF <- data.frame(
  meanA=apply(permTD,1,function(X) mean(testData[[X[1]]])),
  meanB=apply(permTD,1,function(X) mean(testData[[X[2]]])),
  meanOneA=apply(permTD,1,function(X) mean.logX(log2TD[[X[1]]],2,1,pc.out=1)),
  meanOneB=apply(permTD,1,function(X) mean.logX(log2TD[[X[2]]],2,1,pc.out=1)),
  meanCellsA=apply(permTD,1,function(X) mean.logX(log2TD[[X[1]]],2,1,pc.out=1/length(unlist(log2TD)))),
  meanCellsB=apply(permTD,1,function(X) mean.logX(log2TD[[X[2]]],2,1,pc.out=1/length(unlist(log2TD)))),
  meanSmallA=apply(permTD,1,function(X) mean.logX(log2TD[[X[1]]],2,1,pc.out=1e-99)),
  meanSmallB=apply(permTD,1,function(X) mean.logX(log2TD[[X[2]]],2,1,pc.out=1e-99))
)
fcDF$TrueDiff <- fcDF$meanA - fcDF$meanB
fcDF$TrueLogFC <- log2(fcDF$meanA) - log2(fcDF$meanB)
fcDF$fcOne <- fcDF$meanOneA - fcDF$meanOneB
fcDF$fcCells <- fcDF$meanCellsA - fcDF$meanCellsB
fcDF$fcSmall <- fcDF$meanSmallA - fcDF$meanSmallB

fcDF[fcDF$TrueLogFC > -7 & fcDF$TrueLogFC < -6,]

tempMean <- apply(fcDF[,c("meanA","meanB")],1,mean)
tempCut <- cut(tempMean,100,labels=F)

# Fig1a ----
tempPAR <- par(no.readonly=T)

cairo_ps(file="Fig1a.eps",height=6,width=6,fallback_resolution=600)

layout(cbind(2,1,3),widths=c(1,6,1))

par(mar=c(3,0,3,0),mgp=2:0)
plot(fcDF$TrueLogFC,fcDF$fcOne,
     xlim=range(fcDF$TrueLogFC[abs(fcDF$TrueLogFC) < Inf]),
     ylim=range(fcDF$fcOne),
     xlab="True log-ratio of gene abundance",
     #main="Calculating logGER with a pseudocount of 1",
     ylab=NA,yaxt="n",pch=19,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut])
abline(0,1)
temp_x <- seq(from=par("usr")[1] + strwidth(round(min(tempMean),2)) * 3,
              to=par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.3,
              length.out=101)
for (i in 1:100) {
  rect(xleft=temp_x[i],xright=temp_x[i+1],
       ybottom=par("usr")[4] - strheight(round(min(tempMean),2)) * 1.5,
       ytop=par("usr")[4] - strheight(round(min(tempMean),2)) * .5,
       col=colorspace::sequential_hcl(100,palette="Viridis",alpha=1,rev=T)[i],border=NA,xpd=NA)
}
text(temp_x[1],par("usr")[4] - strheight(round(min(tempMean),2)),
     labels=round(min(tempMean),2),pos=2)
text(temp_x[101],par("usr")[4] - strheight(round(min(tempMean),2)),
     labels=round(max(tempMean),2),pos=4)
text(temp_x[51],par("usr")[4] - strheight(round(min(tempMean),2)),
      labels="Mean gene abundance",pos=1)
legend(par("usr")[1],par("usr")[4] - strheight(round(min(tempMean),2)) * 2,
       legend="Pseudocount = 1",pch=19,bty="n")

par(mar=c(3,3,3,0))
plot(rep(1,sum(fcDF$TrueLogFC == -Inf)),fcDF[fcDF$TrueLogFC == -Inf,"fcOne"],
     ylim=range(fcDF$fcOne),
     ylab="Calculated log-ratio of gene abundance",
     xlab=NA,xaxt="n",pch=19,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]])
axis(1,1,expression(-infinity))

par(mar=c(3,0,3,3))
plot(rep(1,sum(fcDF$TrueLogFC == Inf)),fcDF[fcDF$TrueLogFC == Inf,"fcOne"],
     ylim=range(fcDF$fcOne),xlab=NA,xaxt="n",ylab=NA,yaxt="n",pch=19,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut[fcDF$TrueLogFC == Inf]])
axis(1,1,expression(infinity))

dev.off()


# Fig1b ----
cairo_ps(file="Fig1b.eps",height=6,width=6,fallback_resolution=600)
layout(matrix(1:9,nrow=3),widths=c(1,6,1),heights=c(1.5,5,1.5))

par(mar=c(0.5,3,3,0.5),mgp=2:0)
plot(x=NA,y=NA,xlim=c(1,1),ylim=range(fcDF[fcDF$TrueLogFC == Inf,"fcSmall"]),
     xlab=NA,xaxt="n",ylab=NA)

par(mar=c(0.5,3,0.5,0.5))
plot(rep(1,sum(fcDF$TrueLogFC == -Inf)),fcDF[fcDF$TrueLogFC == -Inf,"fcCells"],
     ylim=range(fcDF$fcCells),
     ylab="Calculated log-ratio of gene abundance",
     xlab=NA,xaxt="n",pch=22,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",rev=T)[tempCut[fcDF$TrueLogFC == -Inf]],
     bg=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]])

par(mar=c(3,3,0.5,0.5))
plot(rep(1,sum(fcDF$TrueLogFC == -Inf)),fcDF[fcDF$TrueLogFC == -Inf,"fcSmall"],
     xlab=NA,xaxt="n",ylab=NA,pch=23,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",alpha=1,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]],
     bg=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]])
axis(1,1,expression(-infinity))

par(mar=c(0.5,0.5,3,0.5))
plot(x=NA,y=NA,
     xlim=range(fcDF[abs(fcDF$TrueLogFC) < Inf,"TrueLogFC"]),
     ylim=range(fcDF[fcDF$TrueLogFC == Inf,"fcSmall"]),
     xlab=NA,xaxt="n",ylab=NA,yaxt="n")

par(mar=c(0.5,0.5,0.5,0.5))
plot(fcDF$TrueLogFC,fcDF$fcCells,
     xlab=NA,xaxt="n",ylab=NA,yaxt="n",pch=22,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",alpha=1,rev=T)[tempCut],
     bg=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut])
points(fcDF$TrueLogFC,fcDF$fcSmall,
       pch=23,cex=1.5,
       col=colorspace::sequential_hcl(100,palette="Viridis",alpha=1,rev=T)[tempCut],
       bg=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut])
abline(0,1)
temp_x <- seq(from=par("usr")[1] + strwidth(round(min(tempMean),2)) * 3,
              to=par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.3,
              length.out=101)
for (i in 1:100) {
  rect(xleft=temp_x[i],xright=temp_x[i+1],
       ybottom=par("usr")[4] - strheight(round(min(tempMean),2)) * 1.5,
       ytop=par("usr")[4] - strheight(round(min(tempMean),2)) * .5,
       col=colorspace::sequential_hcl(100,palette="Viridis",alpha=1,rev=T)[i],border=NA,xpd=NA)
}
text(temp_x[1],par("usr")[4] - strheight(round(min(tempMean),2)),
     labels=round(min(tempMean),2),pos=2)
text(temp_x[101],par("usr")[4] - strheight(round(min(tempMean),2)),
     labels=round(max(tempMean),2),pos=4)
text(temp_x[51],par("usr")[4] - strheight(round(min(tempMean),2)),
     labels="Mean gene abundance",pos=1)
legend(par("usr")[1],par("usr")[4] - strheight(round(min(tempMean),2)) * 2,
       legend=c("Pseudocount = 1 / # cells","Pseudocount = 1e-99"),
       pch=c(22,23),col=rep("black",2),pt.bg=scales::alpha(rep("black",2),.5),bty="n")

par(mar=c(3,0.5,0.5,0.5))
plot(x=NA,y=NA,
     xlim=range(fcDF[abs(fcDF$TrueLogFC) < Inf,"TrueLogFC"]),
     ylim=range(fcDF[fcDF$TrueLogFC == -Inf,"fcSmall"]),
     ylab=NA,yaxt="n",xlab="True log-ratio of gene abundance")

par(mar=c(0.5,0.5,3,3))
plot(rep(1,sum(fcDF$TrueLogFC == Inf)),fcDF[fcDF$TrueLogFC == Inf,"fcSmall"],
     xlab=NA,xaxt="n",ylab=NA,yaxt="n",pch=23,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",alpha=1,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]],
     bg=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]])

par(mar=c(0.5,0.5,0.5,3))
plot(rep(1,sum(fcDF$TrueLogFC == Inf)),fcDF[fcDF$TrueLogFC == Inf,"fcCells"],
     ylim=range(fcDF$fcCells),xlab=NA,xaxt="n",ylab=NA,yaxt="n",pch=22,cex=1.5,
     col=colorspace::sequential_hcl(100,palette="Viridis",alpha=1,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]],
     bg=colorspace::sequential_hcl(100,palette="Viridis",alpha=.5,rev=T)[tempCut[fcDF$TrueLogFC == -Inf]])

par(mar=c(3,0.5,0.5,3))
plot(x=NA,y=NA,xlim=c(1,1),ylim=range(fcDF[fcDF$TrueLogFC == -Inf,"fcSmall"]),
     xlab=NA,xaxt="n",ylab=NA,yaxt="n")
axis(1,1,expression(infinity))

dev.off()


# Mean plot ----
par(tempPAR)
par(mar=c(3,3,2,2),mgp=2:0)
meanDF <- fcDF[,c("meanA","meanOneA","meanCellsA","meanSmallA")]
meanDF <- meanDF[!duplicated(meanDF$meanA) & meanDF$meanA > 0,]
names(meanDF) <- c("True","psOne","psCells","psSmall")

plot(psOne~True,data=meanDF,log="x",
     ylim=range(meanDF[,2:4]),
     pch=21,col="black",bg=scales::alpha("black",.5),cex=1.5,
     xlab="True mean (log scale)",ylab="Calculated log2 mean")
points(psCells~True,data=meanDF,
       pch=22,col=rainbow2(2)[2],bg=rainbow2(2,.5)[2],cex=1.5)
points(psSmall~True,data=meanDF,
       pch=23,col=rainbow2(2)[1],bg=rainbow2(2,.5)[1],cex=1.5)
legend("topleft",legend=c("Pseudoount = 1","Pseudocount = 1 / # cells","Pseudocount = 1e-99"),
       bty="n",pch=c(21,23,22),col=c("black",rainbow2(2)),
       pt.bg=c(scales::alpha("black",.5),rainbow2(2,.5)))
