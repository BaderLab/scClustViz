## See https://baderlab.github.io/scRNAseq_meCortex/PipelineV2/pseudocountTest
library(scales)
library(RColorBrewer)

mean.logX <- function(data,ex=2,pc=1,pc.out) { log(mean(ex^data - pc) + pc.out,base=ex) }
halfway <- 10

testData <- t(sapply(seq(0,.099,.001),function(prob) rbinom(1e3,100,prob)))
fc <- rowMeans(testData)[halfway] / rowMeans(testData)
testMeans <- rowMeans(testData)
logTestData <- log2(testData + 1)


pseudocount.use <- c(1,1e-99,1/nrow(testData))
testLogMeans <- sapply(pseudocount.use,function(X) apply(logTestData,MARGIN=1,FUN=mean.logX,pc.out=X))
logFC <- apply(testLogMeans,2,function(X) X[halfway] - X)

pdf(file="../Reports/18_F1000Res_scClustViz/Fig1a.pdf",height=6,width=6)
layout(matrix(c(3,1,4,2),2),widths=c(6,1),heights=c(1,6))
par(mar=c(3,3,.1,.1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=range(log2(fc)[-1]),ylim=c(min(logFC),max(logFC[,-2])),
     xlab="log2(fold-change)",ylab="LogFC with pseudocount")
for (i in 1:ncol(logFC)) {
  points(log2(fc)[-1],logFC[-1,i],pch=19,col=alpha(brewer.pal(3,"Dark2")[i],.5))
}
abline(0,1)
legend("topleft",bty="n",pch=19,col=alpha(brewer.pal(3,"Dark2"),.5),
       legend=c("1","1e-99","1/#cells"),title="Pseudocount")
par(mar=c(3,.1,.1,1))
plot(x=c(1,1),y=logFC[1,-2],ylim=c(min(logFC),max(logFC[,-2])),
     xaxt="n",yaxt="n",xlab=NA,
     pch=19,col=alpha(brewer.pal(3,"Dark2"),.5)[c(1,3)])
axis(1,1,expression(infinity))
par(mar=c(.1,3,1,.1))
plot(x=1,y=logFC[1,2],type="n",xaxt="n",ylab=NA)
par(mar=c(.1,.1,1,1))
plot(x=1,y=logFC[1,2],xaxt="n",yaxt="n",
     pch=19,col=alpha(brewer.pal(3,"Dark2"),.5)[2])
dev.off()


pdf(file="../Reports/18_F1000Res_scClustViz/Fig1b.pdf",height=3,width=6)
layout(rbind(1:2),widths=c(6,1))
par(mar=c(3,3,1,.1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=c(min(logFC),max(logFC[,-2])),ylim=c(.5,3.5),
     yaxt="n",ylab="Pseudocount",xlab="LogFC with pseudocount")
axis(2,1:3,c("1","1e-99","1/#cells"))
for (i in 1:ncol(logFC)) {
  if (i == 2) { drop2 <- -1 } else { drop2 <- 1:nrow(logFC) }
  boxplot(logFC[drop2,i],col=alpha(brewer.pal(3,"Dark2")[i],.5),
          add=T,at=i,horizontal=T,xaxt="n")
}
par(mar=c(3,.1,1,1))
plot(x=logFC[1,2],y=2,ylim=c(.5,3.5),yaxt="n",xlab=NA)
dev.off()
