library(cluster)

rainbow2 <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 60, c = 100)[1:n]
}

mean.log2 <- function(X) { log2(mean(2^X - 1) + 1) }


sd <- readRDS("testData.rds")
sd@ident <- sd@meta.data$res.0.6
names(sd@ident) <- rownames(sd@meta.data)

DN <- apply(sd@data,1,function(X) tapply(X,sd@ident,function(Y) sum(Y > 0)))
DR <- apply(sd@data,1,function(X) tapply(X,sd@ident,function(Y) sum(Y > 0) / length(Y)))
MTC <- apply(sd@data,1,function(X) tapply(X,sd@ident,mean.log2))
MDTC <- apply(sd@data,1,function(X) tapply(X,sd@ident,function(Y) 
  if (all(Y == 0)) { return(0) } else { mean.log2(Y[Y > 0]) }))
MedTC <- apply(sd@data,1,function(X) tapply(X,sd@ident,median))
Q90TC <- apply(sd@data,1,function(X) tapply(X,sd@ident,function(Y) quantile(Y,.9)))

#dist_param <- apply(DN,1,function(X) MASS::fitdistr(X,"negative binomial"))
#par(mar=c(0,0,0,0),mfrow=c(3,3))
#for (i in 1:nrow(DN)) {
#  plot(quantile(DN[1,],seq(.01,.99,.01)),
#       qnbinom(seq(.01,.99,.01),size=dist_param[[1]]$estimate["size"],mu=dist_param[[1]]$estimate["mu"]))
#  lines(range(DN[1,]),range(DN[1,]))     
#}

DRthresh <- apply(combn(1:nrow(DR),2),2,function(X) DR[X[1],] > 0.1 & DR[X[2],] > 0.1)
diffDR <- apply(combn(1:nrow(DR),2),2,function(X) DR[X[1],] - DR[X[2],])
diffDR2 <- apply(combn(1:nrow(DR),2),2,function(X) log2(DR[X[1],]) - log2(DR[X[2],]))
diffVarTC <- apply(combn(1:nrow(DR),2),2,function(Y) apply(sd@data,1,function(X) var(X[sd@ident %in% Y])))
diffMTC <- apply(combn(1:nrow(MTC),2),2,function(X) MTC[X[1],] - MTC[X[2],])
diffMDTC <- apply(combn(1:nrow(MDTC),2),2,function(X) MDTC[X[1],] - MDTC[X[2],])
diffMedTC <- apply(combn(1:nrow(MedTC),2),2,function(X) MedTC[X[1],] - MedTC[X[2],])
diffQ90TC <- apply(combn(1:nrow(MedTC),2),2,function(X) Q90TC[X[1],] - Q90TC[X[2],])
diffU <- apply(combn(levels(sd@ident),2),2,function(Y) 
  pbapply::pbapply(sd@data,1,function(X) wilcox.test(X[sd@ident == Y[1]],X[sd@ident == Y[2]])$p.value))
diffU[is.na(diffU)] <- 1


boxplot(list(MTC=sapply(1:ncol(diffMTC),function(i) cor(abs(diffMTC[,i]),-log10(diffU[,i]))),
             MDTC=sapply(1:ncol(diffMDTC),function(i) cor(abs(diffMDTC[,i]),-log10(diffU[,i]))),
             Q90TC=sapply(1:ncol(diffMTC),function(i) cor(abs(diffQ90TC[,i]),-log10(diffU[,i]))),
             Var=sapply(1:ncol(diffMTC),function(i) cor(diffVarTC[,i],-log10(diffU[,i]))),
             DR=sapply(1:ncol(diffMTC),function(i) cor(abs(diffDR[,i]),-log10(diffU[,i])))))

#par(mar=c(3,3,2,1),mgp=2:0,mfrow=c(1,1))
par(mar=c(0,0,0,0),mfcol=c(2,2))
for (i in c(1,8)) {
  plot(abs(diffMTC[,i]),-log10(diffU[,i]),pch=".",cex=2,
       xaxt="n",yaxt="n")
       #main=paste(combn(1:nrow(DR),2)[,i],collapse="-"),
       #xlab=expression(Delta~"Detection Rate"),ylab=expression(Log[2]~"Fold Change"))
  abline(h=-log10(c(0.01,0.001,0.0001)),lty=3:1,col="red")
  abline(v=.25,lty=3,col="blue")
  legend("bottomright",legend=paste("logFC",paste(combn(1:nrow(DR),2)[,i],collapse="-")),bty="n")
  box(col=rainbow2(ncol(diffDR))[i])

  plot(abs(diffDR[,i]),-log10(diffU[,i]),pch=".",cex=2,
       xaxt="n",yaxt="n")
  #main=paste(combn(1:nrow(DR),2)[,i],collapse="-"),
  #xlab=expression(Delta~"Detection Rate"),ylab=expression(Log[2]~"Fold Change"))
  abline(h=-log10(c(0.01,0.001,0.0001)),lty=3:1,col="red")
  abline(v=.1,lty=3,col="blue")
  legend("bottomright",legend=paste("diffDR",paste(combn(1:nrow(DR),2)[,i],collapse="-")),bty="n")
  box(col=rainbow2(ncol(diffDR))[i])
}

par(mar=c(3,3,1,1),mgp=2:0,mfrow=c(1,1),las=3)
barplot(sapply(1:ncol(diffDR),function(i) {
  bot <- which(abs(diffU[,i]) < 0.0001)
  top <- which(abs(diffQ90TC[,i]) > .2)
  return(c(length(setdiff(bot,top)),length(intersect(bot,top)),length(setdiff(top,bot))))
}),col=c("blue","purple","red"),names.arg=apply(combn(1:nrow(DR),2),2,function(X) paste(X,collapse="-")),
legend.text=c("pVal < 0.01","both",expression(Delta~"Q90 > 0.1")),args.legend=list(x="topleft",bty="n",horiz=T))

