testData <- sapply(seq(0,.1,.001),function(prob) rbinom(1e6,100,prob))
fc <- colMeans(testData) / colMeans(testData)[10]

logTestData <- log(testData + 1)


pseudocount.use <- 1 # Seurat default
testLogMeans <- apply(X=logTestData,MARGIN=2,FUN=function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
logFC <- testLogMeans - testLogMeans[10]

hist(testLogMeans,breaks=20)
plot(log(fc),logFC,pch=19,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),
     xlab="log(fold-change)",ylab="LogFC with Seurat's pseudocount",main="Seurat Default")
abline(0,1)


pseudocount.use <- 0 # Set to zero in "FindMarkers" arguments
testLogMeans <- apply(X=logTestData,MARGIN=2,FUN=function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
logFC <- testLogMeans - testLogMeans[10]

hist(testLogMeans,breaks=20)
plot(log(fc),logFC,pch=19,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),
     xlab="log(fold-change)",ylab="LogFC with no pseudocount",main="After setting pseudocount.use=0")
abline(0,1)


pseudocount.use <- 1e-99 # Set in "FindMarkers" arguments
testLogMeans <- apply(X=logTestData,MARGIN=2,FUN=function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
logFC <- testLogMeans - testLogMeans[10]

hist(testLogMeans,breaks=20)
plot(log(fc),logFC,pch=19,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),
     xlab="log(fold-change)",ylab="LogFC with no pseudocount",main="After setting pseudocount.use=1e-99")
abline(0,1)

pseudocount.use <- 1/nrow(logTestData) # Set in "FindMarkers" arguments
testLogMeans <- apply(X=logTestData,MARGIN=2,FUN=function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use))
logFC <- testLogMeans - testLogMeans[10]

hist(testLogMeans,breaks=20)
plot(log(fc),logFC,pch=19,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),
     xlab="log(fold-change)",ylab="LogFC with no pseudocount",main="After setting pseudocount.use=1e-99")
abline(0,1)






test <- rbind(rep(0,500),c(1,(rep(0,499))),c(2,(rep(0,499))),c(1,1,(rep(0,498))))
apply(X=test,MARGIN=1,FUN=function(x) log(x = mean(x = expm1(x = x)) + 1))
apply(X=test,MARGIN=1,FUN=function(x) log(x = mean(x = expm1(x = x)) + 0))
apply(X=test,MARGIN=1,FUN=function(x) log(x = mean(x = expm1(x = x)) + 1e-99))
apply(X=test,MARGIN=1,FUN=function(x) log(x = mean(x = expm1(x = x)) + 1/ncol(test)))
