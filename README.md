# scClustViz
An interactive R Shiny tool for visualizing single-cell RNAseq clustering results from common analysis pipelines.  Its main goal is two-fold: **A:** to help select a biologically appropriate resolution or K from clustering results by assessing differential expression between the resulting clusters; and **B:** help annotate cell types and identify marker genes.  
[Our preprint is on F1000Research!](https://f1000research.com/articles/7-1522/v2)  

### Contents:
- [Example Output](#example-output)
- [Usage](#usage)
  - [Installation](#installation)
  - [Basic Usage](#basic-usage)
  - [Iterative Clustering With scClustViz](#iterative-clustering-with-scclustviz)
  - [Use Your Own Differential Expression Results](#use-your-own-differential-expression-results)
- [Data Packages](#data-packages)
  - [Embryonic Mouse Cerebral Cortex](#embryonic-mouse-cerebral-cortex)
  - [Human Liver Atlas](#human-liver-atlas)
  - [Make Your Own Data Package!](#make-your-own-data-package)
- [Citation](#citation)
- [Contact](#contact)

# Example Output
Before installing a package it's always nice to *see* what it is, so here's a summary of our use case scenario from [the paper](https://f1000research.com/articles/7-1522/v2), using the mouse embryonic cerebral cortex data package [outlined below](#embryonic-mouse-cerebral-cortex).
First, you can use metrics such as silhouette width or the number of differentially expressed genes between nearest neighbouring clusters to select the optimal clustering resolution.
![Differentially expressed genes between nearest neighbouring clusters for all clustering results](/inst/paperFigs/Fig2a.png) ![Silhouette plot of selected clustering result](/inst/paperFigs/Fig2b.png)
...


# Usage
## Installation
Install scClustViz using devtools:
```r
# install devtools
install.packages("devtools")

# install scClustViz
devtools::install_github("BaderLab/scClustViz")
```
(If you're on linux and getting errors running `devtools::install_github`, make sure RCurl is working - you might need to install libcurl4-openssl-dev).  

## Basic Usage
Following normalization, dimensionality reduction (include 2D cell embedding), and clustering using a workflow of your choice, scClustViz can be used to do differential expression testing (using the Wilcoxon rank-sum test) to both assess different clustering solutions and explore your results.  First, run the DE testing as follows:
```r
library(scClustViz)

# if using Seurat, this regex can grab 
# the metadata columns representing cluster results:
your_cluster columns <- grepl("res[.0-9]+$",
                              names(getMD(your_scRNAseq_data_object)))
your_cluster_results <- getMD(
  your_scRNAseq_data_oubject[your_cluster_columns]
)


sCVdata_list <- CalcAllSCV(
  inD=your_scRNAseq_data_object,
  clusterDF=your_cluster_results,
  assayData=NULL, #specify assay slot of data
  DRforClust="pca",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.1, #gene filter - minimum detection rate
  testAll=F, #stop testing clusterings when no DE between clusters
  FDRthresh=0.05,
  calcSil=T, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=T
)

save(your_scRNAseq_data_object,sCVdata_list,
     file="for_scClustViz.RData")
# This file can now be shared so anyone 
# can view your results with the Shiny app!
```
Once the previous setup step has been performed once and the output saved, you can explore the data in the interactive Shiny interface by simply pointing it to the saved file:
```r
# Lets assume this is data from an embryonic mouse cerebral cortex:
# (This is the call wrapped by MouseCortex::viewMouseCortex("e13"))
runShiny(
  filePath="for_scClustViz.RData",
  
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  
  annotationDB="org.Mm.eg.db",
  # This is an optional argument, but will add annotations.
  
  cellMarkers=list("Cortical precursors"=c("Mki67","Sox2","Pax6",
                                           "Pcna","Nes","Cux1","Cux2"),
                   "Interneurons"=c("Gad1","Gad2","Npy","Sst","Lhx6",
                                    "Tubb3","Rbfox3","Dcx"),
                   "Cajal-Retzius neurons"="Reln",
                   "Intermediate progenitors"="Eomes",
                   "Projection neurons"=c("Tbr1","Satb2","Fezf2",
                                          "Bcl11b","Tle4","Nes",
                                          "Cux1","Cux2","Tubb3",
                                          "Rbfox3","Dcx")
  ),
  # This is a list of canonical marker genes per expected cell type.
  # The app uses this list to automatically annotate clusters.
  
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)

```

## Iterative Clustering With scClustViz
Incorporating the scClustViz cluster assessment metric into your analysis 
pipeline is simply a matter of running the differential expression testing after 
every clustering run, instead of post-hoc. This allows you to systematically 
increase the resolution or K parameter of the clustering algorithm until 
statistically significant differential expression between nearest neighbour 
clusters is lost. An example using the Seurat clustering method is shown here.
```r
DE_bw_clust <- TRUE
seurat_resolution <- 0
sCVdata_list <- list()

while(DE_bw_clust) {
  seurat_resolution <- seurat_resolution + 0.2
  # ^ Iteratively incrementing resolution parameter

  your_seurat_obj <- Seurat::FindClusters(your_seurat_obj,
                                          resolution=seurat_resolution)
  # ^ Calculate clusters using method of choice.

  curr_sCVdata <- CalcSCV(
    inD=your_seurat_obj,
    cl=your_seurat_obj@ident, #factor containing cluster assignments
    assayData=NULL, #specify assay slot of data
    DRforClust="pca", #reduced dimensions for silhouette calc
    exponent=exp(1), #log base of normalized data
    pseudocount=1,
    DRthresh=0.1, #gene filter - minimum detection rate
    calcSil=T, #use cluster::silhouette to calc silhouette widths
    calcDEvsRest=T,
    calcDEcombn=T
  )

  DE_bw_NN <- sapply(DEneighb(curr_sCVdata,0.05),length)
  # ^ counts # of DE genes between neighbouring clusters at 5% FDR

  if (min(DE_bw_NN) < 1) { DE_bw_clust <- FALSE }
  # ^ If no DE genes between nearest neighbours, don't loop again.

  sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
  # Add sCVdata object to list with an appropriate name.
}

save(your_seurat_obj,sCVdata_list,
     file="for_scClustViz.RData")

runShiny(filePath="for_scClustViz.RData")
# ^ see ?runShiny for detailed argument list
```

## Use Your Own Differential Expression Results
scClustViz uses the wilcoxon rank-sum test for its differential expression testing.  You can provide your own DE results from a testing method of your choice instead, skipping sCV's testing steps. In both `CalcAllSCV` and `CalcSCV` there are arguments `calcDEvsRest` and `calcDEcombn`, which can be set to false to skip those differential expression calculations.  You can then use `DEvsRest(your_sCVdata_object) <- your_DE_dataframe_list` and `DEcombn(your_sCVdata_object) <- your_DE_dataframe_list` to pass your results into the sCVdata objects. DEvsRest represents differential expression tests between each cluster and the remaining cells, and should be a named list of data frames where each name refers to the tested cluster (see `?CalcDEvsRest` for details). DEcombn represents differential expression tests between all pairwise combinations of clusters, and should be a named list of data frames were each name refers to the cluster pair, with cluster names separated by "-" (see `?CalcDEcombn` for details). In both cases, data frames must contain variables `logGER` (an effect size measure: gene expression ratio in log space, often referred to as logFC) and `FDR` (significance measure: false discovery rate), as well as `dDR` (an effect size measure: difference in detection rate) for `DEcombn`.

# Data Packages
The following data packages can be used to explore the features of scClustViz. You can also follow the vignette below to build your own data package to easily share your analysed scRNAseq data with collaborators and the public.

## Embryonic Mouse Cerebral Cortex
The data from the 2017 Cell Reports paper [Developmental Emergence of Adult Neural Stem Cells as Revealed by Single-Cell Transcriptional Profiling](https://doi.org/10.1016/j.celrep.2017.12.017) by Yuzwa *et al.* are available to explore by installing the R package [MouseCortex](https://github.com/BaderLab/MouseCortex). These are DropSeq data from timepoints spanning neurogenesis and filtered for cortically-derived cells, processed on an earlier version of the pipeline outlined below (using scran for normalization and Seurat for clustering) and imported into scClustViz using the steps outlined above.

Install MouseCortex using devtools as follows:
```r
# install devtools
install.packages("devtools")

# install MouseCortex (demo data from Yuzwa et al, Cell Reports 2017)
devtools::install_github("BaderLab/MouseCortex") 
# this takes a minute or two

# install mouse gene annotations from bioconductor (optional)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
```
Then run the scClustViz Shiny app to view your dataset of choice! 
There's a wrapper function in the MouseCortex package that handles the call to scClustViz, so it's nice and simple. 
If you're interested, `?runShiny` has example code showing the function call used by the wrapper function.
```r
library(MouseCortex)
viewMouseCortex("e13")
```

## Human Liver Atlas
The data from the 2018 Nature Communications paper [Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations](https://doi.org/10.1038/s41467-018-06318-7) by MacParland *et al.* are available to explore by installing the R package [HumanLiver](https://github.com/BaderLab/HumanLiver). These are 10X Chromium data from the livers of 5 human donors, processed on the pipeline outlined below (using scran for normalization and Seurat for clustering) and imported into scClustViz using the steps outlined above.

Install HumanLiver using devtools as follows:
```r
# install devtools
install.packages("devtools")

# install HumanLiver 
# (R data package for MacParland et al., Nat Commun 2018)
devtools::install_github("BaderLab/HumanLiver") 
# this takes a minute or two

# install human gene annotations from bioconductor (optional)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
```
Then run the scClustViz Shiny app to view your dataset of choice! 
There's a wrapper function in the MouseCortex package that handles the call to scClustViz, so it's nice and simple. 
If you're interested, `?runShiny` has example code showing the function call used by the wrapper function.
```r
library(HumanLiver)
viewHumanLiver()
```

## Make Your Own Data Package!
Building an R package is a relatively easy task thanks to RStudio and the roxygen2 and devtools packages. The following vignette will show you how to take your saved output from the scClustViz setup and share it as an R package on github as seen in the data packages above. It is entirely based on the invaluable book [R packages](http://r-pkgs.had.co.nz/) by Hadley Wickham.  
First, you must have generated your input file for the `runShiny` command in scClustViz by following the steps in the [usage guide](#scclustviz-usage-guide) above.  
Then, create a new project in RStudio, selecting "New directory" -> "R package" and making sure to check "Create a git repository". If you haven't already set up git/github in RStudio, check out [this blogpost](https://www.r-bloggers.com/rstudio-and-github/) for an explanation. If you only want to make a package to share with colleagues, you can skip github and simply send them the bundled package when you're done.  
Once you've opened your new package in RStudio, make sure to have both "Use devtools package functions" and "Generate documentation with Roxygen" selected under "Project Options" -> "Build Tools". Also, delete the existing NAMESPACE file, since Roxygen will create a new one when you build the package.  
You're now ready to build your package.  First, make a folder in the package directory called "inst", and put your input file for *runShiny* there. All files in "inst" become part of the root directory of the package after installation, so it's best to store your data in a folder within inst.
```r
dir.create("inst/packageData/",recursive=T)
save(data_for_scClustViz,DE_for_scClustViz,
     file="inst/packageData/MyDataTitle.RData")
```
If you'd like a default resolution to load when the user views your data in scClustViz, now's the time to save that.
```r
runShiny("inst/packageData/MyDataTitle.RData")
```
Save your selected cluster resolution as default in the app. It will be saved as "inst/packageData/MyDataTitle_savedRes.RData". You will also see a file called "inst/packageData/MyDataTitle_intro.md". This is a markdown file that stores the text displayed at the top of the scClustViz GUI. You can edit it to say what you want (perhaps a link to the paper the data is from, and maybe the abstract?).  
Now all you need to do is write the wrapper function to call *runShiny*. Here is an example R script (overwrite R/HelloWorld.R) to save in the "R" directory of the package.
```r
#' View MyData data in the scClustViz Shiny app
#'
#' A wrapper function to view the \code{MyData} dataset in the
#' \code{scClustViz} Shiny app.
#'
#' @param outPath Default = "./" (the working directory). Specify the 
#'   directory used to save/load any analysis files you generate while 
#'   exploring the \code{MyData} data.
#'
#' @return The function causes the scClustViz Shiny GUI app to open in a
#'   seperate window.
#'
#' @examples
#'   viewMyData()
#'
#' @seealso \url{https://baderlab.github.io/scClustViz} for information on
#'   \code{scClustViz}.
#'
#' @export

viewMyData <- function(outPath="./",imageFileType="pdf") {
  filePath <- system.file("packageData/MyDataTitle.RData",
                          package="MyDataPackage")
  cellMarkers <- list()
  # If you have a list of cell-type marker genes for you data,
  # add them here!
  
  # Change "org.Hs.eg.db" to the appropriate AnnotationDbi object for your 
  # data. This way if your user has the library installed, it will be used, 
  # otherwise it will be skipped without causing any errors.
  if (require("org.Hs.eg.db",quietly=T)) {
    annotationDB <- org.Hs.eg.db
    scClustViz::runShiny(filePath=filePath,
                         outPath=outPath,
                         cellMarkers=cellMarkers,
                         annotationDB=annotationDB,
                         imageFileType=imageFileType)

  } else {
    scClustViz::runShiny(filePath=filePath,
                         outPath=outPath,
                         cellMarkers=cellMarkers,
                         imageFileType=imageFileType)
  }
}
```
Now that you have a wrapper function, all that's left to do is fix up the DESCRIPTION file. The most important entries for functionality in the file are the following:
```
Suggests: org.Hs.eg.db
Imports: scClustViz
Remotes: BaderLab/scClustViz
```
Change "org.Hs.eg.db" to the appropriate AnnotationDbi library. This lets the user know that they would benefit from having it installed. More importantly, `Imports: scClustViz` tells R devtools to install scClustViz when installing your package. Since scClustViz isn't in CRAN, the line `Remotes: BaderLab/scClustViz` lets devtools know where to find it.  
Now that everything's ready, use the "Install and Restart" button in RStudio or hit Ctrl+Shift+B to build and install the package locally. You should now be able to use the wrapper command to open scClustViz with your data. If you're happy with everything, it's time to push to github!  
First you must [create a new repository on github](https://help.github.com/articles/creating-a-new-repository/) for your package. Then it's as simple as pushing your first commit (commands here are in the bash shell):
```sh
# Set the remote to the github account:
git remote add origin https://github.com/YourGithubAccount/MyDataPackage.git 

# Stage your directory
git add .

# Make your first commit
git commit -m "MyData is now an R package!"

# Push your first commit to github 
# (could be slow, since you're uploading data files)
git push -u origin master
```
Now all you need to do is edit the README file to tell the world how to install and run your package:
```r
devtools::install_github("YourGithubAccount/MyDataPackage")
library(MyDataPackage)
viewMyData()
```


# Citation
Innes BT and Bader GD. scClustViz – Single-cell RNAseq cluster assessment and visualization [version 2; peer review: 1 approved, 1 approved with reservations]. F1000Research 2019, 7:1522 (doi: [10.12688/f1000research.16198.2](https://doi.org/10.12688/f1000research.16198.2))


# Contact
You can [contact me](http://www.baderlab.org/BrendanInnes) for questions about this repo.  For general scRNAseq questions, do what I do and [ask the Toronto single-cell RNAseq working group on Slack](http://bit.ly/scRNAseqTO)!  



