library(pbapply)
library(RColorBrewer)


######### DropSeq ########
if (!require(MouseCortex)) {
  devtools::install_github("BaderLab/MouseCortex")
  library(MouseCortex)
}

#### Data setup ####
load(system.file("e13/e13.RData",package="MouseCortex"))

res <- "res.0.8"
combos <- combn(levels(ds$cl[,res]),2)
colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
deM_dDR <- apply(combos,2,function(i) dd$CGS[[res]][[i[1]]]$DR - dd$CGS[[res]][[i[2]]]$DR)
rownames(deM_dDR) <- rownames(dd$CGS[[res]][[i[1]]])
deM_logGER <- apply(combos,2,function(i) dd$CGS[[res]][[i[1]]]$MTC - dd$CGS[[res]][[i[2]]]$MTC)
rownames(deM_logGER) <- rownames(dd$CGS[[res]][[i[1]]])

pVal_unfiltered <- pbsapply(colnames(combos),function(i)
  apply(ds$nge,1,function(X) 
    wilcox.test(X[ds$cl[,res] == combos[1,i]],
                X[ds$cl[,res] == combos[2,i]])$p.value),simplify=T)
pVal_unfiltered[is.na(pVal_unfiltered)] <- 1


#### Actual experiment (DropSeq) ####

#par(mar=c(3,3,3,1),mgp=2:0)
#plot(as.vector(deM_logGER),-log10(as.vector(pVal_unfiltered)),pch=".",
#     ylab="-log10(p-value)",xlab="Log2(Gene Expression Ratio)")
#plot(as.vector(deM_dDR),-log10(as.vector(pVal_unfiltered)),pch=".",
#     ylab="-log10(p-value)",xlab=expression(Delta~"Detection Rate"))

TPR <- function(thresh,method) {
  sum(abs(as.vector(method)) > thresh & as.vector(pVal_unfiltered) < 0.01) / sum(as.vector(pVal_unfiltered) < 0.01)
}

FPR <- function(thresh,method) {
  sum(abs(as.vector(method)) > thresh & !(as.vector(pVal_unfiltered) < 0.01)) / sum(!as.vector(pVal_unfiltered) < 0.01)
}

Pr <- function(thresh,method) {
  sum(abs(as.vector(method)) > thresh & as.vector(pVal_unfiltered) < 0.01) / sum(abs(as.vector(method)) > thresh)
}

dDR <- data.frame(TPR=sapply(seq(0,1,.01),TPR,method=deM_dDR),
                  FPR=sapply(seq(0,1,.01),FPR,method=deM_dDR),
                  Pr=sapply(seq(0,1,.01),Pr,method=deM_dDR))
rownames(dDR) <- seq(0,1,.01)
GER <- data.frame(TPR=sapply(seq(0,10,.05),TPR,method=deM_logGER),
                  FPR=sapply(seq(0,10,.05),FPR,method=deM_logGER),
                  Pr=sapply(seq(0,10,.05),Pr,method=deM_logGER))
rownames(GER) <- seq(0,10,.05)

bp <- brewer.pal(3,"Dark2")

par(mar=c(3,3,3,1),mgp=2:0)
plot(Pr~TPR,data=dDR,type="l",col=bp[1],ylab="Precision",xlab="Recall",lwd=2)
lines(Pr~TPR,data=GER,type="l",col=bp[2],lwd=2)
points(Pr~TPR,data=dDR[c("0.1","0.15","0.2"),],pch=19,col=bp[1])
text(dDR[c("0.1","0.15","0.2"),"TPR"],dDR[c("0.1","0.15","0.2"),"Pr"],c("0.1","0.15","0.2"),col=bp[1],pos=4)
legend("top",bty="n",lwd=2,col=bp[1:2],horiz=T,inset=c(0,-.12),xpd=NA,
       legend=c(expression(Delta~"Detection Rate"),"Gene Expression Ratio"),title="Threshold type")




######### 10X Chromium ########
have10X <- require(HumanLiver) 
#10X data from here (private repo until related paper is published in Sept/Oct 2018, ask author for permission)

if (have10X) {
  #### Data setup ####
  load(system.file("liver/HumanLiver.RData",package="HumanLiver"))
  
  res <- "res.0.8"
  combos <- combn(levels(ds$cl[,res]),2)
  colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
  deM_dDR <- apply(combos,2,function(i) dd$CGS[[res]][[i[1]]]$DR - dd$CGS[[res]][[i[2]]]$DR)
  rownames(deM_dDR) <- rownames(dd$CGS[[res]][[i[1]]])
  deM_logGER <- apply(combos,2,function(i) dd$CGS[[res]][[i[1]]]$MTC - dd$CGS[[res]][[i[2]]]$MTC)
  rownames(deM_logGER) <- rownames(dd$CGS[[res]][[i[1]]])
  
  pVal_unfiltered <- pbsapply(colnames(combos),function(i)
    apply(ds$nge,1,function(X) 
      wilcox.test(X[ds$cl[,res] == combos[1,i]],
                  X[ds$cl[,res] == combos[2,i]])$p.value),simplify=T)
  pVal_unfiltered[is.na(pVal_unfiltered)] <- 1
  
  #### Actual experiment ####
  par(mar=c(3,3,3,1),mgp=2:0)
  
  plot(as.vector(deM_logGER),-log10(as.vector(pVal_unfiltered)),pch=".",
       ylab="-log10(p-value)",xlab="Log2(Gene Expression Ratio)")
  plot(as.vector(deM_dDR),-log10(as.vector(pVal_unfiltered)),pch=".",
       ylab="-log10(p-value)",xlab=expression(Delta~"Detection Rate"))
  
  
  TPR <- function(thresh,method) {
    sum(abs(as.vector(method)) > thresh & as.vector(pVal_unfiltered) < 0.01) / sum(as.vector(pVal_unfiltered) < 0.01)
  }
  
  FPR <- function(thresh,method) {
    sum(abs(as.vector(method)) > thresh & !(as.vector(pVal_unfiltered) < 0.01)) / sum(!as.vector(pVal_unfiltered) < 0.01)
  }
  
  Pr <- function(thresh,method) {
    sum(abs(as.vector(method)) > thresh & as.vector(pVal_unfiltered) < 0.01) / sum(abs(as.vector(method)) > thresh)
  }
  
  dDR2 <- data.frame(TPR=sapply(seq(0,1,.01),TPR,method=deM_dDR),
                     FPR=sapply(seq(0,1,.01),FPR,method=deM_dDR),
                     Pr=sapply(seq(0,1,.01),Pr,method=deM_dDR))
  rownames(dDR2) <- seq(0,1,.01)
  GER2 <- data.frame(TPR=sapply(seq(0,10,.05),TPR,method=deM_logGER),
                     FPR=sapply(seq(0,10,.05),FPR,method=deM_logGER),
                     Pr=sapply(seq(0,10,.05),Pr,method=deM_logGER))
  rownames(GER2) <- seq(0,10,.05)
  
  bp <- brewer.pal(3,"Dark2")
  
  plot(TPR~FPR,data=dDR2,type="l",col=bp[1],ylab="True Positive Rate",xlab="False Positive Rate",lwd=2)
  lines(TPR~FPR,data=GER2,type="l",col=bp[2],lwd=2)
  legend("top",bty="n",lwd=2,col=bp[1:2],horiz=T,inset=c(0,-.12),xpd=NA,
         legend=c(expression(Delta~"Detection Rate"),"Gene Expression Ratio"),title="Threshold type")
  
  pdf(file="../Reports/18_F1000Res_scClustViz/Fig2b.pdf",width=6,height=6)
  par(mar=c(3,3,3,1),mgp=2:0)
  plot(Pr~TPR,data=dDR2,type="l",col=bp[1],ylab="Precision",xlab="Recall",lwd=2)
  lines(Pr~TPR,data=GER2,type="l",col=bp[2],lwd=2)
  points(Pr~TPR,data=dDR2[c("0.1","0.15","0.2"),],pch=19,col=bp[1])
  text(dDR2[c("0.1","0.15","0.2"),"TPR"],dDR2[c("0.1","0.15","0.2"),"Pr"],c("0.1","0.15","0.2"),col=bp[1],pos=4)
  legend("top",bty="n",lwd=2,col=bp[1:2],horiz=T,inset=c(0,-.12),xpd=NA,
         legend=c(expression(Delta~"Detection Rate"),"Gene Expression Ratio"),title="Threshold type")
  dev.off()
}


######## Combined Figure ########

pdf(file="../Reports/18_F1000Res_scClustViz/Fig2.pdf",width=6,height=6)
par(mar=c(3,3,4,1),mgp=2:0)
plot(x=NULL,y=NULL,xlim=0:1,ylim=0:1,ylab="Precision",xlab="Recall")

lines(Pr~TPR,data=dDR,type="l",col=bp[1],lwd=2,lty=2)
points(Pr~TPR,data=dDR[c("0.1","0.15","0.2"),],pch=19,col=bp[1])
text(dDR[c("0.1","0.15","0.2"),"TPR"],dDR[c("0.1","0.15","0.2"),"Pr"],c("0.1","0.15","0.2"),col=bp[1],pos=4)
lines(Pr~TPR,data=GER,type="l",col=bp[1],lwd=2,lty=1)
points(Pr~TPR,data=GER[c("0.25","1","1.5"),],pch=19,col=bp[1])
text(GER[c("0.25","1","1.5"),"TPR"],GER[c("0.25","1","1.5"),"Pr"],c("0.25","1","1.5"),col=bp[1],pos=1)

lines(Pr~TPR,data=dDR2,type="l",col=bp[2],lwd=2,lty=2)
points(Pr~TPR,data=dDR2[c("0.1","0.15","0.2"),],pch=19,col=bp[2])
text(dDR2[c("0.1","0.15","0.2"),"TPR"],dDR2[c("0.1","0.15","0.2"),"Pr"],c("0.1","0.15","0.2"),col=bp[2],pos=1)
lines(Pr~TPR,data=GER2,type="l",col=bp[2],lwd=2,lty=1)
points(Pr~TPR,data=GER2[c("0.25","1","1.5"),],pch=19,col=bp[2])
text(GER2[c("0.25","1","1.5"),"TPR"],GER2[c("0.25","1","1.5"),"Pr"],c("0.25","1","1.5"),col=bp[2],pos=1)

legend("top",bty="n",ncol=2,inset=c(0,-.16),xpd=NA,
       lwd=2,lty=c(1,1,2,2),col=bp[c(1,2,1,2)],title="Threshold type",
       legend=c("Gene Expression Ratio (DropSeq)","Gene Expression Ratio (10X)",
                expression(Delta~"Detection Rate (DropSeq)"),expression(Delta~"Detection Rate (10X)")))
dev.off()
