######## User-defined variables ########

dataPath <- "meCortex/e13/e13_Cortical_Only_forViz.RData"
##  ^ Point this to the output file from PrepareInputs.R
##  If you set a default resolution in the Shiny app, it will save to the same directory.

vizScriptPath <- "./" 
##  ^ Point this to the directory in which the "app.R" Shiny script resides

species <- "mouse" 
##  ^ Set species ("mouse"/"human").  
##  If other, add the annotation database from Bioconductor to the egDB <- switch() expression below.


#### List known cell-type markers ####
cellMarkers <- list("Cortical precursors"=c("Mki67","Sox2","Pax6","Pcna","Nes","Cux1","Cux2"),
                    "Interneurons"=c("Gad1","Gad2","Npy","Sst","Lhx6","Tubb3","Rbfox3","Dcx"),
                    "Cajal-Retzius neurons"="Reln",
                    "Intermediate progenitors"="Eomes",
                    "Projection neurons"=c("Tbr1","Satb2","Fezf2","Bcl11b","Tle4",
                                           "Nes","Cux1","Cux2","Tubb3","Rbfox3","Dcx"))
#cellMarkers <- list()
##  ^ If you have canonical marker genes for expected cell types, list them here 
##  (see example above from mouse embryonic cortex).  The Shiny app will attempt 
##  to label clusters in the tSNE projection by highest median gene expression.
##  Otherwise leave the list blank (uncomment line above).


#### Variables for differential expression analysis ####
exponent <- 2  
##  ^ log base of your normalized input data.  
##  Seurat defaults to natural log (set this to exp(1)), 
##  other methods are generally log2 (set this to 2).
pseudocount <- 1 
##  ^ pseudocount added to all log-normalized values in your input data.  
##  Most methods use a pseudocount of 1 to eliminate log(0) errors.

#threshType <- "logGER"  # use a fold-change-based threshold for filtering genes prior to DE testing
threshType <- "dDR"     # use a difference in detection rate threshold for filtering 
##  Filtering genes for use in differential expression testing can be done multiple ways.
##  We use a fold-change filter for comparing each cluster to the tissue as a whole, but find that
##  difference in detection rates works better when comparing clusters to each other.  You can set
##  threshType to "logGER" to use fold-change for all gene filtering if you'd prefer.

logGERthresh <- 1  # magnitude of mean log-expression fold change between clusters to use as filter.
dDRthresh <- 0.15 # magnitude of detection rate difference between clusters to use as filter.
WRSTalpha <- 0.01 # significance level for DE testing using Wilcoxon rank sum test


########################################



######## Code to run the Shiny app ########
library(markdown)
library(shiny)
library(cluster)
library(gplots)
library(scales)
library(viridis)
library(RColorBrewer)
library(TeachingDemos)

egDB <- switch(species,
               mouse={ requireNamespace(org.Mm.eg.db); "org.Mm.eg.db" },
               human={ requireNamespace(org.Hs.eg.db); "org.Hs.eg.db" },
               stop("
Set species please!  
If not mouse/human, add your species' annotation database from Bioconductor:  
source('https://bioconductor.org/biocLite.R')
biocLite('org.Xx.eg.db')
"))

meanLogX <- function(data,ex=exponent,pc=pseudocount) { log(mean(ex^data - pc) + 1/ncol(nge),base=ex) }
rainbow2 <- function(n,a=1) {
  require(scales)
  hues = seq(15, 375, length = n + 1)
  alpha(hcl(h = hues, l = 60, c = 100)[1:n],a)
}

if (length(cellMarkers) < 1) {
  cellMarkersS <- cellMarkersU <- list()
} else {
  cellMarkersS <- apply(combn(seq_along(cellMarkers),2),2,
                        function(X) do.call(intersect,unname(cellMarkers[X])))
  try(names(cellMarkersS) <- apply(combn(seq_along(cellMarkers),2),2,
                                   function(X) paste(X,collapse="&")),silent=T)
  cellMarkersS <- cellMarkersS[sapply(cellMarkersS,length) > 0]
  cellMarkersU <- lapply(cellMarkers,function(X) X[!X %in% unlist(cellMarkersS)])
}

demoRegex <- switch(species,mouse="^Actb$",human="^ACTB$")


load(dataPath)
temp_dataPath <- strsplit(dataPath,"/|\\\\")
dataPath <- sub(temp_dataPath[[1]][length(temp_dataPath[[1]])],"",dataPath)
if (dataPath == "") { dataPath <- "./" }
dataTitle <- sub("\\..+$|_forViz\\..+$","",temp_dataPath[[1]][length(temp_dataPath[[1]])])
rm(temp_dataPath)

for (selDEfile in grep(paste0("^",dataTitle,".+selDE.+RData$"),list.files(dataPath),value=T)) {
  temp <- load(paste0(dataPath,selDEfile))
  cl <- cbind(cl,new_cl)
  CGS <- append(CGS,new_CGS)
  deTissue <- append(deTissue,new_deTissue)
  deMarker <- append(deMarker,new_deMarker)
  rm(list=temp)
}

if (file.exists(paste0(dataPath,dataTitle,"_savedRes.RData"))) {
  load(paste0(dataPath,dataTitle,"_savedRes.RData"))
} else {
  savedRes <- NULL
}

if (!file.exists(paste0(dataPath,"intro.md"))) {
  write(paste0(dataTitle,": You can add to this preamble by editting ",dataPath,"intro.md"),
        file=paste0(dataPath,"intro.md"))
}

silDist <- dist(dr_clust)
##  ^ precalculating distances in reduced dimensionality space for the silhouette plot.

for (l in names(CGS)) {
  for (i in names(CGS[[l]])) {
    CGS[[l]][[i]]$MTCrank <- rank(CGS[[l]][[i]]$MTC,ties.method="min")/nrow(CGS[[l]][[i]])
    if (i == "Unselected") { next }
    CGS[[l]][[i]]$cMu <- rownames(CGS[[l]][[i]]) %in% unlist(cellMarkersU)
    CGS[[l]][[i]]$cMs <- rownames(CGS[[l]][[i]]) %in% unlist(cellMarkersS)
    CGS[[l]][[i]]$overCut <- CGS[[l]][[i]]$MTC > mean(CGS[[l]][[i]]$MTC)
    CGS[[l]][[i]]$genes <- rownames(CGS[[l]][[i]])
  }
}

if (length(cellMarkers) < 1) {
  clusterID <- sapply(names(CGS),function(X) sapply(CGS[[X]],function(Z) return("")),simplify=F)
} else if (!any(unlist(cellMarkers) %in% rownames(nge))) {
  warning(paste("None of the provided cellMarkers are found in the data",
                "(check your gene IDs against rownames in your data)."))
  clusterID <- sapply(names(CGS),function(X) sapply(CGS[[X]],function(Z) return("")),simplify=F)
} else {
  clusterID <- sapply(CGS,function(Z) {
    temp <- names(cellMarkers)[sapply(Z,function(Y) 
      which.max(sapply(cellMarkers,function(X) median(Y$MTC[rownames(Y) %in% X]))))]
    names(temp) <- names(Z)
    temp[names(temp) == "Unselected"] <- "Unselected"
    return(temp)
  },simplify=F)
}

#### Run the Shiny App!  ####
runApp(vizScriptPath)
