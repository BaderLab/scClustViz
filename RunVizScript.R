######## User-defined variables ########

#### Change this to reflect dataset of interest ####
dataPath <- "D:/Dropbox/GDB/SamBrain/code_outputs/180117_DS37/"

#### Point this to the directory in which the "app.R" Shiny script resides ####
vizScriptPath <- "D:/Dropbox/GDB/code/ShinyViz/VizScript_v1/" 

#### Set species ####
species <- "mouse" 

#### List known cell-type markers ####
cellMarkers <- list("pNSCs"=c("Pou5f1","Sox2","Lifr"),
                    "dNSCs"=c("Gfap","Sox2","Prom1","Egfr"),
                    "Neuronal precursors"=c("Ascl1","Myt1l"),
                    "Post-mitotic neurons"=c("Map2","Tubb3","Dcx"),
                    "Oligodendrocytes"=c("Olig2","Olig1","Pdgfra","Foxo4","Mbp"),
                    "Astrocytes"=c("Gfap","Slc1a3","S100b"))

########################################



######## Code to run the Shiny app ########
library(shiny)
library(Seurat) #see http://satijalab.org/seurat/install.html
library(cluster)
library(gplots)
library(viridis)
library(RColorBrewer)
library(scales)
library(TeachingDemos)
library(vioplot)

if (species == "mouse") {
  library(org.Mm.eg.db) # from Bioconductor
  egDB <- "org.Mm.eg.db"
} else if (species == "human") {
  library(org.Hs.eg.db) # from Bioconductor
  egDB <- "org.Hs.eg.db"
} else { print("Set species please!") }

cellMarkersS <- apply(combn(seq_along(cellMarkers),2),2,function(X) do.call(intersect,unname(cellMarkers[X])))
names(cellMarkersS) <- apply(combn(seq_along(cellMarkers),2),2,function(X) paste(X,collapse="&"))
cellMarkersS <- cellMarkersS[sapply(cellMarkersS,length) > 0]
cellMarkersU <- lapply(cellMarkers,function(X) X[!X %in% unlist(cellMarkersS)])

if (!exists("eb1S")) {
  load(paste0(dataPath,"eb1S.RData"))
}

if (!exists("cycScores")) {
  load(paste0(dataPath,"cycScores.RData"))
}

if (file.exists(paste0(dataPath,"savedRes.RData"))) {
  load(paste0(dataPath,"savedRes.RData"))
} else {
  savedRes <- NULL
}

gg_colour_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

runApp(vizScriptPath)
