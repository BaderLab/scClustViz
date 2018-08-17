# scClustViz
An interactive R Shiny tool for visualizing single-cell RNAseq clustering results from the *Seurat* R package or any other analysis pipeline.  Its main goal is two-fold: **A:** to help select a biologically appropriate resolution or K from clustering results by assessing differential expression between the resulting clusters; and **B:** help annotate cell types and identify marker genes.

-   [Quick Start](#quick-start)
-   [scClustViz Usage Guide](#scclustviz-usage-guide)  
    -   [Read in data](#read-in-data)  
    -   [Differential expression testing](#differential-expression-testing)  
    -   [Run the Shiny app](#run-the-shiny-app)  
-   [Demo with Embryonic Mouse Cerebral Cortex Data](#demo-with-embryonic-mouse-cerebral-cortex-data)
-   [scRNAseq analysis pipeline](#scrnaseq-analysis-pipeline)  
-   [Contact](#contact)  

## Quick Start
Install scClustViz using devtools:
```{r}
# install devtools
install.packages("devtools")

# install scClustViz
devtools::install_github("BaderLab/scClustViz")
```
Load data from your Seurat analysis for differential expression testing and visualization in scClustViz:
```{r}
library(scClustViz)

data_for_scClustViz <- readFromSeurat(your_seurat_object)
rm(your_seurat_object)
# All the data scClustViz needs is in 'data_for_scClustViz'.

DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz,exponent=exp(1))

save(data_for_scClustViz,DE_for_scClustViz,
     file="for_scClustViz.RData")
# Save these objects so you'll never have to run this slow function again!

runShiny(filePath="for_scClustViz.RData")
```

## scClustViz Usage Guide
scClustViz takes the output object from your single-cell analysis pipeline of choice, and runs differential expression testing for all the clustering solutions generated during your analysis to generate a cluster assessment metric used in the visualization tool. The visualization tool itself is an R Shiny app that generates a variety of figures designed to help assess clustering results, and identify clusters and their marker genes.

### Read in data
scClustViz assumes you have tried a variety of parameterizations when clustering the cells from your scRNAseq data, and want to decide which clustering solution you should use (if you haven't yet clustered your data, or are interested in an example of integrating the differential expression metric used in this tool to systematically test different clustering resolutions, see the example [pipeline below](#scrnaseq-analysis-pipeline)).  
To read in your data from a Seurat object (check the documentation to ensure your object meets requirements), you can run:
```{r}
data_for_scClustViz <- readFromSeurat(your_seurat_object)
```
If your data isn't in a Seurat object, or otherwise doesn't fit the requirements for `readFromSeurat` you can run `readFromManual`, which allows you to manually add all the required components of your analysis to the object scClustViz uses for the differential expression testing. See its man page (`?readFromManual`) for details, or use the example here using a hypothetical SingleCellExperiment class from Bioconductor as the input:
```{r}
# A logical vector separating the cluster assignments from the rest of the
# cell metadata in the colData slot. This is an example that you will have
# to change to reflect your cluster assignment column names.
clusterAssignments <- grepl("^Clust",colnames(colData(mySCE)))

data_for_scClustViz <- readFromManual(nge=logcounts(mySCE),
                                      md=colData(mySCE)[,!clusterAssignments],
                                      cl=colData(mySCE)[,clusterAssignments],
                                      dr_clust=reducedDim(mySCE,"PCA"),
                                      dr_viz=reductedDim(mySCE,"tSNE"))
# All the data scClustViz needs is in 'data_for_scClustViz'.
```
### Differential expression testing
*A more thorough explanation of the DE testing scheme and how to bypass it (structure of the output lists in case you want to replace it with your own DE method/results) will be here soon. For now, see `?clusterWiseDEtest`*
```{r}
DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz,
                                       # Stop once DE is lost between nearest neighbouring clusters
                                       testAll=FALSE,
                                       # Normalized data is in log2 space
                                       exponent=2,
                                       # Pseudocount of 1 was added to log-normalized data
                                       pseudocount=1,
                                       # False discovery rate threshold of 1%
                                       FDRthresh=0.01,
                                       # Use difference in detection rate to filter genes for testing
                                       threshType="dDR",
                                       # Genes with at least 15% detection rate difference will be tested
                                       dDRthresh=0.15
                                       )

# Save the results of the preprocessing for use in the Shiny app!
save(data_for_scClustViz,DE_for_scClustViz,file="for_scClustViz.RData")
```

### Run the Shiny app
Finally, its time to run the app. Running this function will open the Shiny UI in a separate window.  Have fun exploring your data!
```{r}
runShiny(filePath="for_scClustViz.RData")
```

## Demo with Embryonic Mouse Cerebral Cortex Data
The data from the 2017 Cell Reports paper [Developmental Emergence of Adult Neural Stem Cells as Revealed by Single-Cell Transcriptional Profiling](https://doi.org/10.1016/j.celrep.2017.12.017) by Yuzwa *et al.* are available to explore by installing the R package [MouseCortex](https://github.com/BaderLab/MouseCortex). These are DropSeq data from timepoints spanning neurogenesis and filtered for cortically-derived cells, processed on an earlier version of the pipeline outlined below (using scran for normalization and Seurat for clustering) and imported into scClustViz using the steps outlined above.

Install scClustViz and MouseCortex using devtools as follows:
```{r}
# install devtools
install.packages("devtools")
# install scClustViz
devtools::install_github("BaderLab/scClustViz")
# install MouseCortex (demo data from Yuzwa et al, Cell Reports 2017)
devtools::install_github("BaderLab/MouseCortex") # this takes a few minutes
# install mouse cell annotations from bioconductor (optional)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
```
Then run the scClustViz Shiny app to view your dataset of choice! Replace "e13" in the filename with e11, e15, or e17 to view your timepoint of interest.
```{r}
library(scClustViz)
runShiny(
         # Load input file (E13.5 data) from package directory.
         filePath=system.file("e13cortical.RData",package="MouseCortex"),
         # Save any further analysis performed in the app to the 
         # working directory rather than library directory.
         outPath=".",
         # This is an optional argument, but will add annotations.
         annotationDB="org.Mm.eg.db",
         # This is a list of canonical marker genes per expected cell type.
         # The app uses this list to automatically annotate clusters.
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
                          )
         
         )

```


## scRNAseq analysis pipeline
*Currently being updated to use the functions from scClustViz - check back soon.*

## Contact
You can [contact me](http://www.baderlab.org/BrendanInnes) for questions about this repo.  For general scRNAseq questions, do what I do and [ask the Toronto single-cell RNAseq working group on Slack](http://bit.ly/scRNAseqTO)!  



