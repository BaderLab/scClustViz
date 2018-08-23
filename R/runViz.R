#' Run the scClustViz Shiny app
#'
#' Performs differential expression testing between clusters for all cluster
#' solutions in order to assess the biological relevance of each cluster
#' solution. Differential expression testing is done using the Wilcoxon rank-sum
#' test implemented in the base R \code{stats} package. For details about what
#' is being compared in the tests, see the "Value" section.
#'
#' @param filePath A character vector giving the relative filepath to an RData
#'   file containing two objects. One must be the list outputted by one of the
#'   importData functions (either \code{\link{readFromSeurat}} or
#'   \code{\link{readFromManual}}) containing the data for viewing in the app.
#'   The other must be the list outputted by the \code{\link{clusterWiseDEtest}}
#'   function containing differential gene expression results for viewing in the
#'   app. As long as none of the name of the list elements have been changed,
#'   the objects can be named anything you'd like. Note that any files generated
#'   by the Shiny app (ie. saving the selected cluster solution, saving custom
#'   set DE testing results) will be saved/loaded in the same directory as the
#'   input file.
#'
#' @param outPath Optional. If you'd like to save/load any analysis files
#'   to/from a different directory than the input directory (for example, if
#'   you're using data from a package), specify that directory here.
#'
#' @param cellMarkers Optional. If you have canonical marker genes for expected
#'   cell types, list them here (see example code below). Note that the gene
#'   names must match rownames of your data (ie. use ensembl IDs if your gene
#'   expression matrix rownames are ensembl IDs). The Shiny app will attempt to
#'   label clusters in the tSNE projection by highest median gene expression.
#'
#' @param annotationDB Optional. An AnnotationDbi object for your data's species
#'   (ie. org.Mm.eg.db / org.Hs.eg.db for mouse / human respectively). If
#'   present, gene names will be shown in gene-specific figures, official gene
#'   symbols (instead of your rownames) will be displayed in figures, and gene
#'   searches performed using both official gene symbols and your rownames. If
#'   the gene IDs in your data aren't official gene symbols, using this argument
#'   will make the visualization tool much more useful.
#'
#' @param rownameKeytype Optional. A character vector indicating the
#'   AnnotationDbi keytype (see \code{AnnotationDbi::keytypes(annotationDB)})
#'   that represents your rownames. If the annotationDB argument is present and
#'   this is missing, the function will assume the rownames are official gene
#'   symbols. If less than 80% of rownames map to official gene symbols, the
#'   function will try to predict the appropriate keytype of the rownames (this
#'   takes a bit of time).
#'
#' @param exponent Default = Taken from \code{clusterWiseDEtest} output. The log
#'   base of your normalized input data. Seurat normalization uses the natural
#'   log (set this to exp(1)), while other normalization methods generally use
#'   log2 (set this to 2). This is used if you use the function for testing
#'   differential gene expression between custom sets, and is set automatically
#'   to match the parameters used in \code{clusterWiseDEtest}.
#'
#' @param pseudocount Default = Taken from \code{clusterWiseDEtest} output. The
#'   pseudocount added to all log-normalized values in your input data. Most
#'   methods use a pseudocount of 1 to eliminate log(0) errors. This is used if
#'   you use the function for testing differential gene expression between
#'   custom sets, and is set automatically to match the parameters used in
#'   \code{clusterWiseDEtest}.
#'
#' @param FDRthresh Default = Taken from \code{clusterWiseDEtest} output. The
#'   false discovery rate to use as a threshold for determining statistical
#'   significance of differential expression calculated by the Wilcoxon rank-sum
#'   test. This is used if you use the function for testing differential gene
#'   expression between custom sets, and is set automatically to match the
#'   parameters used in \code{clusterWiseDEtest}.
#'
#' @param threshType Default = Taken from \code{clusterWiseDEtest} output.
#'   Filtering genes for use in differential expression testing can be done
#'   multiple ways. We use an expression ratio filter for comparing each cluster
#'   to the rest of the tissue as a whole, but find that difference in detection
#'   rates works better when comparing clusters to each other. You can set
#'   threshType to \code{"logGER"} to use a gene expression ratio for all gene
#'   filtering, or leave it as default (\code{"dDR"}) to use difference in
#'   detection rate as the thresholding method when comparing clusters to each
#'   other. This is used if you use the function for testing differential gene
#'   expression between custom sets, and is set automatically to match the
#'   parameters used in \code{clusterWiseDEtest}.
#'
#' @param dDRthresh Default = Taken from \code{clusterWiseDEtest} output.
#'   Magnitude of detection rate difference of a gene between clusters to use as
#'   filter for determining which genes to test for differential expression
#'   between clusters. This is used if you use the function for testing
#'   differential gene expression between custom sets, and is set automatically
#'   to match the parameters used in \code{clusterWiseDEtest}.
#'
#' @param logGERthresh Default = Taken from \code{clusterWiseDEtest} output.
#'   Magnitude of gene expression ratio for a gene between clusters to use as
#'   filter for determining which genes to test for differential expression
#'   between clusters. This is used if you use the function for testing
#'   differential gene expression between custom sets, and is set automatically
#'   to match the parameters used in \code{clusterWiseDEtest}.
#'
#' @return The function causes the scClustViz Shiny GUI app to open in a
#'   seperate window.
#'
#' @examples
#' \dontrun{
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object)
#'  rm(your_seurat_object)
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#'
#'  DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz)
#'
#'  save(data_for_scClustViz,DE_for_scClustViz,
#'       file="for_scClustViz.RData")
#'  # Save these objects so you'll never have to run this slow function again!
#'
#'  runShiny(filePath="for_scClustViz.RData")
#'
#'  ### Using example data from the MouseCortex package ###
#'  devtools::install_github("BaderLab/MouseCortex")
#'  viewMouseCortex("e13") 
#'  # MouseCortex has a wrapper function which calls the following:
#'  runShiny(system.file("e13/e13.RData",package="MouseCortex"),
#'           # Load input file (E13.5 data) from package directory.
#'           outPath="./",
#'           # Save any further analysis performed in the app to the
#'           # working directory rather than library directory.
#'           annotationDB="org.Mm.eg.db",
#'           # This is an optional argument, but will add annotations.
#'           cellMarkers=list("Cortical precursors"=c("Mki67","Sox2","Pax6",
#'                                                    "Pcna","Nes","Cux1","Cux2"),
#'                            "Interneurons"=c("Gad1","Gad2","Npy","Sst","Lhx6",
#'                                             "Tubb3","Rbfox3","Dcx"),
#'                            "Cajal-Retzius neurons"="Reln",
#'                            "Intermediate progenitors"="Eomes",
#'                            "Projection neurons"=c("Tbr1","Satb2","Fezf2",
#'                                                   "Bcl11b","Tle4","Nes",
#'                                                   "Cux1","Cux2","Tubb3",
#'                                                   "Rbfox3","Dcx")
#'                            )
#'           # This is a list of canonical marker genes per expected cell type.
#'           # The app uses this list to automatically annotate clusters.
#'           )
#' }
#'
#' @seealso \code{\link{readFromSeurat}} or \code{\link{readFromManual}} for
#'   reading in data to generate the first input object for this function, and
#'   \code{\link{clusterWiseDEtest}} to do the differential expression testing
#'   to generate the second input object for this function.
#'
#' @import shiny
#' @importFrom scales alpha
#' @importFrom viridis viridis
#'
#' @export

runShiny <- function(filePath,outPath,
                     cellMarkers=list(),
                     annotationDB,rownameKeytype,
                     exponent,pseudocount,FDRthresh,
                     threshType,dDRthresh,logGERthresh) {
  # ^ Load data from file ------------------------------------------------------------------
  while(T) {
    if (exists(".lastFileCall")) {
      if (names(.lastFileCall) == filePath) {
        if (exists(.lastFileCall[[1]][1]) & exists(.lastFileCall[[1]][2])) {
          break
        } else {
          rm(.lastFileCall,envir=.GlobalEnv)
        }
      } else {
        rm(.lastFileCall,envir=.GlobalEnv)
      }
    } else {
      .lastFileCall <<- list(load(filePath,envir=.GlobalEnv))
      names(.lastFileCall) <<- filePath
      break
    }
  }
  # The above weird-ass loop (or weird ass-loop if you prefer) checks to see if
  # the file has already been loaded (if this function has been run previously
  # this session), otherwise loads the file.
  
  temp_objNames <- sapply(.lastFileCall[[filePath]],function(X) names(get(X)),simplify=F)
  for (L in names(temp_objNames)) {
    for(N in names(get(L))) {
      assign(N,get(L)[[N]])
    }
  }
  rm(temp_objNames)
  # Unpacks the two objects from the input file, which are lists of objects
  # needed in the Shiny app, and saves the objects in the function environment
  # under the names the shiny app expects.
  
  # Load parameters from clusterWiseDEtest output
  if (missing(exponent)) { exponent <- params$exponent }
  if (missing(pseudocount)) { pseudocount <- params$pseudocount }
  if (missing(FDRthresh)) { FDRthresh <- params$FDRthresh }
  if (missing(threshType)) { threshType <- params$threshType }
  if (missing(dDRthresh)) { dDRthresh <- params$dDRthresh }
  if (missing(logGERthresh)) { logGERthresh <- params$logGERthresh }

  cl <- cl[names(deNeighb)]
  # Ensures that only clusters that were tested for differential expression are
  # displayed. This prevents a whole pile of errors.
  
  # ^^ dataPath & dataTitle --------------------------------------------------------------
  temp_dataPath <- strsplit(filePath,"/|\\\\")
  dataPath <- sub(temp_dataPath[[1]][length(temp_dataPath[[1]])],"",filePath)
  if (dataPath == "") { dataPath <- "./" }
  dataTitle <- sub("\\.[^.]+$","",tail(temp_dataPath[[1]],1))
  rm(temp_dataPath)
  if (!missing(outPath)) {
    if (!grepl("[/\\]$",outPath)) { outPath <- paste0(outPath,"/") }
  }
  
  # Seperates the file name (which becomes the dataTitle) from the path (which
  # becomes dataPath). This is used when saving and loading various things in
  # the app (default cluster solution, custom set DE results).
  
  # ^^ Load saved comparisons (if any) ---------------------------------------------------
  if (!missing(outPath)) { #Load from both dataPath and outPath if outPath exists.
    for (selDEfile in grep(paste0("^",dataTitle,".+selDE.+RData$"),list.files(outPath),value=T)) {
      temp <- load(paste0(outPath,selDEfile))
      cl <- cbind(cl,new_cl)
      CGS <- append(CGS,new_CGS)
      deTissue <- append(deTissue,new_deTissue)
      deMarker <- append(deMarker,new_deMarker)
      rm(list=temp)
    }
  }
  for (selDEfile in grep(paste0("^",dataTitle,".+selDE.+RData$"),list.files(dataPath),value=T)) {
    temp <- load(paste0(dataPath,selDEfile))
    cl <- cbind(cl,new_cl)
    CGS <- append(CGS,new_CGS)
    deTissue <- append(deTissue,new_deTissue)
    deMarker <- append(deMarker,new_deMarker)
    rm(list=temp)
  }
  #### MEMORY USAGE SHENANIGANS? ####
  # Load in the custom DE results from the same directory as the input file and
  # add them to the relevant lists. THIS MIGHT DOUBLE MEMORY USAGE since the
  # objects are now being modified from the versions stored in the list objects
  # stored in the global environment, which might mean assigning them to a
  # separate address in memory. CHECK ON THIS!
  
  # ^^ Load default cluster solution (if any) --------------------------------------------
  savedRes <- NULL
  if (file.exists(paste0(dataPath,dataTitle,"_savedRes.RData"))) {
    load(paste0(dataPath,dataTitle,"_savedRes.RData"))
  }
  if (!missing(outPath)) {
    if (file.exists(paste0(outPath,dataTitle,"_savedRes.RData"))) {
      load(paste0(outPath,dataTitle,"_savedRes.RData"))
    }
  }
  # Load the default cluster solution from the user-provided filepath (outPath)
  # preferentially, otherwise load from the same directory as the input file. If
  # neither exist, set to NULL.  This is actually accomplished by overwriting
  # the initial NULL up to two times if they both exist, but since it's a single
  # character string it's not appreciably slower.
  
  # ^^ Generate blank preamble if none saved ---------------------------------------------
  introPath <- paste0(dataPath,dataTitle,"_intro.md")
  if (!missing(outPath)) {
    if (file.exists(paste0(outPath,dataTitle,"_intro.md"))) {
      introPath <- paste0(outPath,dataTitle,"_intro.md")
    }
  }
  # Generate a section of preamble text in markdown that the user can edit.
  if (!file.exists(introPath)) {
    write(paste0(dataTitle,": You can add to this preamble by editting ",introPath),
          file=introPath)
  }
  
  # Done loading, now set dataPath so that things saved go to the right place.
  if (!missing(outPath)) { dataPath <- outPath }
  
  # ^ Map rownames to gene symbol ----------------------------------------------------
  if (!missing(annotationDB)) {
    if (is.character(annotationDB)) {
      require(annotationDB,quietly=T,character.only=T)
      annotationDB <- get(annotationDB)
    }
    if (missing(rownameKeytype)) {
      rownameKeytype <- "SYMBOL"
    }
    if (sum(rownames(nge) %in% keys(annotationDB,rownameKeytype)) / nrow(nge) < 0.8) {
      print("Less than 80% of rownames map to official gene symbols.")
      print("Automatically determining keytype from rownames...")
      temp_keyMatch <- pbapply::pbsapply(AnnotationDbi::keytypes(annotationDB),function(X) 
        sum(rownames(nge) %in% AnnotationDbi::keys(annotationDB,X)))
      rownameKeytype <- names(which.max(temp_keyMatch))
      print(paste0("Keytype '",rownameKeytype,"' matched ",
                   max(temp_keyMatch),"/",nrow(nge)," rownames."))
    }
    if (rownameKeytype != "SYMBOL") {
      symbolMap <- mapIds(annotationDB,rownames(nge),"SYMBOL",rownameKeytype,multiVals="first")
    }  
  }
  
  # ^ Cell type annotation from cellMarkers (and other quick calculations for Shiny) -----
  silDist <- dist(dr_clust)
  # Calculate distance in reduced dimensionality space for the silhouette plot.
  
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
  
  for (l in names(CGS)) {
    for (i in names(CGS[[l]])) {
      if (exists("symbolMap")) {
        CGS[[l]][[i]]$genes <- symbolMap[rownames(CGS[[l]][[i]])]
        CGS[[l]][[i]]$genes[is.na(CGS[[l]][[i]]$genes)] <- 
          rownames(CGS[[l]][[i]])[is.na(CGS[[l]][[i]]$genes)]
      } else {
        CGS[[l]][[i]]$genes <- rownames(CGS[[l]][[i]])
      }
      CGS[[l]][[i]]$MTCrank <- rank(CGS[[l]][[i]]$MTC,ties.method="min")/nrow(CGS[[l]][[i]])
      if (i == "Unselected") { next }
      CGS[[l]][[i]]$cMu <- 
        rownames(CGS[[l]][[i]]) %in% unlist(cellMarkersU) |
        CGS[[l]][[i]]$genes %in% unlist(cellMarkersU)
      CGS[[l]][[i]]$cMs <- 
        rownames(CGS[[l]][[i]]) %in% unlist(cellMarkersS) |
        CGS[[l]][[i]]$genes %in% unlist(cellMarkersS)
      CGS[[l]][[i]]$overCut <- CGS[[l]][[i]]$MTC > mean(CGS[[l]][[i]]$MTC)
    }
  }
  
  if (length(cellMarkers) < 1) {
    clusterID <- sapply(names(CGS),function(X) sapply(CGS[[X]],function(Z) return("")),simplify=F)
  } else if (
    if (exists("symbolMap")) {
      !any(unlist(cellMarkers) %in% c(symbolMap,names(symbolMap)))
    } else {
      !any(unlist(cellMarkers) %in% rownames(nge))
    }
  ) {
    warning(paste("None of the provided cellMarkers are found in the data",
                  "(check your gene IDs against rownames in your data)."))
    clusterID <- sapply(names(CGS),function(X) sapply(CGS[[X]],function(Z) return("")),simplify=F)
  } else {
    clusterID <- sapply(CGS,function(Z) {
      temp <- names(cellMarkers)[sapply(Z,function(Y) 
        which.max(sapply(cellMarkers,function(X) median(Y$MTC[Y$genes %in% X]))))]
      names(temp) <- names(Z)
      temp[names(temp) == "Unselected"] <- "Unselected"
      return(temp)
    },simplify=F)
  }
  
  
  
  
  # UI -------------------------------------------------------------------------------------
  ui <- fixedPage(
    fixedRow(
      titlePanel(paste("scClustViz -",dataTitle)),
      includeMarkdown(introPath)
    ),
    hr(),
    
    # ^ Clustering Solution Selection ------------------------------------------------------
    fixedRow(
      titlePanel("Clustering Solution Selection"),
      p(paste(
        "Here you can compare the results of clustering at different resolutions to",
        "determine the appropriate clustering solution for your data. You can see the",
        "cluster solutions represented as boxplots on the left, where each boxplot",
        "represents the number of genes differentially expressed between each cluster",
        "and its nearest neighbour, or marker genes per cluster. The cluster selected",
        "in the pulldown menu is highlighted in red, and the silhouette plot for that",
        "cluster is shown on the right."
      )),
      p(paste(
        "A silhouette plot is a horizontal barplot where each bar is a cell, grouped by",
        "cluster. The width of each bar represents the difference between mean distance",
        "to other cells within the cluster and mean distance to cells in the nearest",
        "neighbouring cluster. Distance is Euclidean in reduced dimensional space.",
        "Positive silhouettes indicate good cluster cohesion."
      )),
      p(paste(
        "Once you've selected an appropriate cluster solution (we suggest picking one",
        "where all nearest neighbouring clusters have differentially expressed genes",
        "between them), click 'View clusters at this resolution' to proceed. If you",
        "want to save this cluster solution as the default for next time, click 'Save",
        "this resolution as default'. All figures can be downloaded as PDFs by clicking",
        "the buttons next to each figure."
      )),
      h1()
    ),
    fixedRow(
      column(6,
             fixedRow(column(6,uiOutput("resSelect"),align="left"),
                      column(6,align="right",
                             actionButton("go","View clusters at this resolution",icon("play"),
                                          style="color: #fff; background-color: #008000"),
                             uiOutput("saveButton")
                      )
             ),
             radioButtons("deType",NULL,
                          list("# of DE genes to nearest neighbouring cluster"="deNeighb",
                               "# of marker genes per cluster"="deMarker"),inline=T),
             plotOutput("cqPlot",height="500px")),
      column(6,plotOutput("sil",height="600px"))
    ),
    fixedRow(
      column(6,downloadButton("cqPlotSave","Save as PDF"),align="left"),
      column(6,downloadButton("silSave","Save as PDF"),align="right")
    ),
    hr(),
    
    # ^ Dataset and Cluster Metadata Inspection --------------------------------------------
    fixedRow(
      titlePanel("Dataset and Cluster Metadata Inspection"),
      p(paste(
        "Here you can explore your dataset as a whole: cluster assignments for all",
        "cells; metadata overlays for cell projections; and figures for comparing",
        "both numeric and categorical metadata. The top two figures show cells",
        "projected into 2D space, where proximity indicates transcriptional similarity.",
        "On the left you can see cluster assignments and the nearest neighbours used in",
        "the differential expression calculations. If cell type marker genes were",
        "provided in RunVizScript.R, it will also show predicted cell type annotations.",
        "On the right you can add a metadata overlay to the cell projection. Below",
        "you can view relationships in the metadata as a scatterplot or compare clusterwise",
        "distributions of metadata as bar- or box-plots. If you select a cluster of interest",
        "(by clicking on a cell in the top-left plot, or from the list two sections down)",
        "it will be highlighted for comparison in these figures."
      )),
      strong(paste(
        "You can select any cluster for further assessment by clicking on a cell",
        "from that cluster in the top-left figure."
      )),
      h1()
    ),
    fixedRow(
      column(6,
             if (length(cellMarkers) > 0 & !all(unlist(clusterID) == "")) {
               radioButtons("tsneLabels","Labels:",inline=T,
                            choices=list("Cluster numbers"="cn",
                                         "Cluster annotations"="ca",
                                         "Cluster annotations (label all)"="can"))
             } else {
               radioButtons("tsneLabels","Labels:",inline=T,
                            choices=list("Cluster numbers"="cn"))
             },
             checkboxInput("nnArrow",value=F,width="100%",
                           label="Show nearest neighbouring clusters by # of DE genes.")
      ),
      
      column(4,selectInput("tsneMDcol",label="Metadata:",width="100%",choices=colnames(md),
                           selected=grep("phase",colnames(md),value=T,ignore.case=T)[1])),
      column(2,uiOutput("tsneMDlog"))
    ),
    fixedRow(
      column(6,plotOutput("tsne",height="580px",click="tsneClick")),
      column(6,plotOutput("tsneMD",height="580px"))
    ),
    fixedRow(
      column(6,align="left",downloadButton("tsneSave","Save as PDF")),
      column(6,align="right",downloadButton("tsneMDSave","Save as PDF"))
    ),
    hr(),
    
    fixedRow(
      column(2,selectInput("mdScatterX","X axis:",choices=colnames(md),selected="total_counts")),
      column(2,selectInput("mdScatterY","Y axis:",choices=colnames(md),selected="total_features")),
      column(2,uiOutput("scatterLog")),
      
      column(3,selectInput("mdFactorData","Metadata:",choices=colnames(md),
                           selected=grep("phase",colnames(md),value=T,ignore.case=T)[1])),
      column(3,uiOutput("mdFactorOpts"))
    ),
    fixedRow(
      column(6,plotOutput("mdScatter",height="560px")),
      column(6,plotOutput("mdFactor",height="560px"))
    ),
    fixedRow(
      column(6,align="left",downloadButton("mdScatterSave","Save as PDF")),
      column(6,align="right",downloadButton("mdFactorSave","Save as PDF"))
    ),
    hr(),
    
    # ^ Differentially Expressed Genes per Cluster -----------------------------------------
    fixedRow(
      titlePanel("Differentially Expressed Genes per Cluster"),
      p(HTML(paste0(
        "Here you can explore the significantly differentially expressed genes per ",
        "cluster. '<b>DE vs Rest</b>' refers to positively differentially expressed genes ",
        "when comparing a cluster to the rest of the cells as a whole. '<b>Marker genes</b>' ",
        "refers to genes positively differentially expressed versus all other clusters ",
        "in a series of pairwise tests. '<b>DE vs neighbour</b>' refers to genes positively ",
        "differentially expressed versus the nearest neighbouring cluster, as measured ",
        "by number of differentially expressed genes between clusters. In all cases, ",
        "Wilcoxon rank-sum tests are used, with a ",FDRthresh * 100,"% false detection ",
        "rate threshold."
      ))),
      p(paste(
        "The dotplot is generated using the differentially expressed genes from the test",
        "and number of genes selected below. A dotplot is a modified heatmap where each",
        "dot encodes both detection rate and average gene expression in detected cells",
        "for a gene in a cluster. Darker colour indicates higher average gene expression",
        "from the cells in which the gene was detected, and larger dot diameter indicates",
        "that the gene was detected in greater proportion of cells from the cluster.",
        "Differentially expressed gene lists can be downloaded as tab-separated text files",
        "by selecting the test type and cluster, and clicking 'Download gene list'.",
        "Genes used in the dotplot can be viewed in the gene expression plots below as well."
      )),
      h1()
      
    ),
    
    fixedRow(
      column(2,uiOutput("heatDEtype")),
      column(6,uiOutput("DEgeneSlider")),
      column(2,uiOutput("DEclustSelect")),
      column(2,downloadButton("deGeneSave","Download gene list"),
             downloadButton("heatmapSave","Save as PDF"),align="right")
    ),
    fixedRow(plotOutput("dotplot",height="600px")),
    hr(),
    
    # ^ Gene Expression Distributions per Cluster ------------------------------------------
    fixedRow(
      titlePanel("Gene Expression Distributions per Cluster"),
      p(paste("Here you can investigate the expression of individual genes per cluster and",
              "across all clusters. The first plot shows mean expression of genes in a cluster",
              "as a function of their detection rate and transcript count when detected. The",
              "x-axis indicates the proportion of cells in the cluster in which each gene was",
              "detected (transcript count > 0), while the y-axis shows the mean normalized",
              "transcript count for each gene from the cells in the cluster in which that gene",
              "was detected. You can select the cluster to view from the menu below, and genes",
              "can be labelled in the figure based on the cell-type markers provided in",
              "RunVizScipt.R, the differentially expressed genes from the selected cluster in",
              "the above heatmap, or by searching for them in the box below the figure.")),
      p(paste("Clicking on the first plot will populate the list of genes near the point clicked,",
              "which can be found above the next figure. By selecting a gene from this list,",
              "you can compare the expression of that gene across all clusters in the second figure.",
              "This list can also be populated using the gene search feature. Plotting options",
              "for the second figure include the option to overlay normalized transcript count",
              "from each cell in the cluster over their respective boxplots ('Include scatterplot'),",
              "and the inclusion of the percentile rank of that gene's expression per cluster as",
              "small triangles on the plot using the right y-axis ('Include gene rank').")),
      h1()
    ),
    fixedRow(
      column(3,uiOutput("genePlotClustSelect")),
      column(9,if (length(cellMarkers) > 0) {
        radioButtons("cgLegend",inline=T,label="Highlighted genes:",
                     choices=c("Cell-type markers"="markers",
                               "Top DE genes (from heatmap)"="heatmap",
                               "Gene symbols from search box below"="search"))
      } else {
        radioButtons("cgLegend",inline=T,label="Highlighted genes:",
                     choices=c("Top DE genes (from heatmap)"="heatmap",
                               "Gene symbols from search box below"="search"))
      })
    ),
    fixedRow(align="right",
             plotOutput("clusterGenes",height="600px",click="cgClick"),
             downloadButton("clusterGenesSave","Save as PDF")
    ),
    
    # ^ Gene expression comparison ---------------------------------------------------------
    fixedRow(
      column(3,radioButtons("searchType",label="Search by:",
                            choices=c("Gene list (comma-separated)"="comma",
                                      "Regular expression"="regex"))),
      column(8,uiOutput("geneSearchBox")),
      column(1,actionButton("GOIgo","Search",icon=icon("search")))
    ),tags$style(type='text/css', "button#GOIgo { margin-top: 25px;  margin-left: -25px; }"),
    fixedRow(
      column(3,radioButtons("boxplotGene",inline=F,
                            label="Genes of interest (to populate list):",
                            choices=c("From click on plots above or below"="click",
                                      "From gene search"="search"))),
      column(4,uiOutput("cgSelect")),
      column(5,checkboxGroupInput("bxpOpts",label="Figure options:",
                                  selected=c("sct","rnk","notch"),inline=T,
                                  choices=list("Include scatterplot"="sct",
                                               "Include gene rank"="rnk",
                                               "Show notch"="notch")))
    ),
    fixedRow(plotOutput("geneTest",height="500px"),
             downloadButton("geneTestSave","Save as PDF")
    ),
    hr(),
    
    # ^ Distribution of genes of interest --------------------------------------------------
    fixedRow(
      titlePanel("Cell Distribution of Genes of Interest"),
      p(paste("Here you can overlay gene expression values for individual genes of interest",
              "on the cell projection. Search for your gene using the search box below,",
              "then select your gene(s) of interest from the dropdown 'Select genes' menu.",
              "You can select multiple genes, but note that for each cell only the gene",
              "expression of the gene with the highest expression in that cell will be displayed.",
              "You have the option to include the cluster labels from the first cell projection",
              "figure in these plots, and to colour the clusters themselves. There are two",
              "copies of this figure for ease of comparison between genes of interest.")),
      h1()
    ),
    fixedRow(
      column(2,radioButtons("searchType1",label="Search by:",
                            choices=c("Gene list"="comma",
                                      "Regular expression"="regex"))),
      column(3,uiOutput("geneSearchBox1")),
      column(1,actionButton("GOI1go","Search",icon=icon("search"))),
      
      column(2,radioButtons("searchType2",label="Search by:",
                            choices=c("Gene list"="comma",
                                      "Regular expression"="regex"))),
      column(3,uiOutput("geneSearchBox2")),
      column(1,actionButton("GOI2go","Search",icon=icon("search")))
    ),tags$style(type='text/css', paste("button#GOI1go { margin-top: 25px; margin-left: -25px; }",
                                        "button#GOI2go { margin-top: 25px; margin-left: -25px; }")),
    
    fixedRow(
      column(3,
             radioButtons("plotClust1",inline=T,label="Plot:",selected="goi",
                          choices=list("Clusters"="clust","Gene expression overlay"="goi")),
             checkboxInput("plotLabel1",label="Include cluster labels (style as above)",value=T)
      ),
      column(3,uiOutput("GOI1select")),
      
      column(3,
             radioButtons("plotClust2",inline=T,label="Plot:",selected="goi",
                          choices=list("Clusters"="clust","Gene expression overlay"="goi")),
             checkboxInput("plotLabel2",label="Include cluster labels (style as above)",value=T)
      ),
      column(3,uiOutput("GOI2select"))
    ),
    
    fixedRow(
      column(6,strong("If multiple genes are selected, the max expression per cell will be displayed")),
      column(6,strong("If multiple genes are selected, the max expression per cell will be displayed"))
    ),
    fixedRow(
      column(6,plotOutput("goiPlot1",height="600px")),
      column(6,plotOutput("goiPlot2",height="600px"))
    ),
    fixedRow(
      column(6,align="left",downloadButton("goiPlot1Save","Save as PDF")),
      column(6,align="right",downloadButton("goiPlot2Save","Save as PDF"))
    ),
    hr(),
    
    # ^ Cluster comparison -----------------------------------------------------------------
    fixedRow(
      titlePanel("Cluster/Set Comparison of Gene Statistics"),
      p(HTML(paste("Here you can directly compare gene expression statistics between clusters.",
                   "Any clusters from the currently selected cluster solution can be compared,",
                   "and you can switch cluster solutions from the menu here for convenience.",
                   "The stats that can be compared are mean normalized transcript count per",
                   "cluster (<b>Mean gene expression</b>), proportion of cells in a cluster in",
                   "which each gene was detected (<b>Detection rate</b>), and mean normalized",
                   "transcript count in cells of the cluster in which the gene was detected",
                   "(<b>Mean detected gene expression</b>). Genes can be labelled based on",
                   "differential expression from the heatmap above, or using the gene search",
                   "feature above."))),
      p(paste("The most different genes in the current comparison can also be labelled. This",
              "calculation can simply be subtracting the gene stat of the x-axis cluster from",
              "that of the y-axis, or distance (residual) from the line of best fit. The latter",
              "calculation may be of value if there is concern that a technical factor such as",
              "library size is confounding a direct comparison between clusters. In either case,",
              "the resulting values can be downloaded as a ranked list where positive values are",
              "higher in the cluster on the y-axis, and negative values are higher in the x-axis",
              "cluster. Since this list ranks all genes in the experiment, it could be used as an",
              "input for GSEA.")),
      h1()
    ),
    fixedRow(
      column(7,plotOutput("setScatter",height="640px",click="scatterClick")),
      column(5,
             fixedRow(
               column(10,uiOutput("resSelect2")),
               column(2,actionButton("go2","View",icon("play"),
                                     style="color: #fff; background-color: #008000"))
             ),
             fixedRow(column(12,uiOutput("saveButton2"))),
             fixedRow(
               column(6,uiOutput("setScatterY")),
               column(6,uiOutput("setScatterX"))
             ),
             fixedRow(
               column(8,radioButtons("scatterInput",label="Gene stat to display:",
                                     choices=c("Mean gene expression"="MTC",
                                               "Detection rate"="DR",
                                               "Mean detected gene expression"="MDTC"))),
               column(4,radioButtons("scatterLine",label="Difference calculation:",
                                     choices=c("Subtraction"="sub","From line of best fit"="lbf")))
             ),
             fixedRow(
               column(8,radioButtons("diffLabelType",label="Label genes by:",
                                      choices=c("Most different by calculation"="diff",
                                                "Top DE genes (from heatmap)"="de",
                                                "Genes symbols from search box above"="search"))),
               column(4,checkboxGroupInput("scatterLabelAngle",label="Plot options:",
                                           choices=c("Flip label angle"="flip")))
               ),
             fixedRow(column(12,uiOutput("diffLabelSelect"))),
             fixedRow(
               column(4,downloadButton("setScatterSave","Save as PDF")),
               column(4,downloadButton("setComparisonSave","Download ranked list"))
             )
      )
    ),tags$style(type='text/css',paste("button#go2 { margin-top: 25px;  margin-left: -25px; }",
                                       "button#updateForViz2 { margin-top: -25px; }",
                                       "button#setComparisonSave { margin-left: -25px; }")),
    
    hr(),
    
    # ^ Custom sets for DE -----------------------------------------------------------------
    fixedRow(titlePanel("Manually Select Cells for DE Testing")),
    fixedRow(
      column(8,plotOutput("tsneSelDE",brush="tsneBrush",height="750px")),
      column(4,
             p(paste("Here you can select cells to further explore using the figures above.",
                     "Click and drag to select cells, and use the buttons below to add them",
                     "to a set of cells. When your sets are ready, name the comparison and",
                     "click the 'Calculate differential gene expression' button. Once the",
                     "calculation is done the comparison will be added to the cluster list",
                     "at the top of the page and the current cluster solution will be updated",
                     "to show this comparison. The comparison can be saved by clicking 'Save",
                     "this comparison to disk' next to either cluster solution menu.")),
             hr(),
             selectInput("tsneSelDEcol","Metadata overlay:",choices=c("",colnames(md))),
             hr(),
             column(6,htmlOutput("textSetA"),
                    actionButton("addCellsA","Set A: Add Cells",icon("plus"),
                                 style="color: #fff; background-color: #a50026"),
                    actionButton("removeCellsA","Set A: Remove Cells",icon("minus"),
                                 style="color: #a50026; background-color: #fff; border-color: #a50026")
             ),
             column(6,htmlOutput("textSetB"),
                    actionButton("addCellsB","Set B: Add Cells",icon("plus"),
                                 style="color: #fff; background-color: #313695"),
                    actionButton("removeCellsB","Set B: Remove Cells",icon("minus"),
                                 style="color: #313695; background-color: #fff; border-color: #313695")
             ),
             htmlOutput("textOverlap"),
             hr(),
             textInput("DEsetName","Short name for this comparison:",
                       placeholder="A-z0-9 only please"),
             actionButton("calcDE","Calculate differential gene expression",icon("play")),
             hr(),
             span(textOutput("calcText"),style="color:red")
      )
    ),
    h1()
  )
  
  

  # Server -------------------------------------------------------------------------------
  server <- function(input,output,session) {
    d <- reactiveValues(cl=cl,CGS=CGS,
                        clusterID=clusterID,
                        deTissue=deTissue,
                        deMarker=deMarker)
    
    clustCols <- function(res) {
      if (grepl("^Comp",res)) {
        c(RColorBrewer::brewer.pal(3,"PRGn")[c(1,3)],"grey80")
      } else if (length(levels(d$cl[,res])) <= 8) {
        RColorBrewer::brewer.pal(length(levels(d$cl[,res])),"Dark2")[1:length(levels(d$cl[,input$res]))]
      } else {
        rainbow2(length(levels(d$cl[,res])))
      }
    }
    
    
    # ^ Clustering Solution Selection ------------------------------------------------------
    # ^^ Inter-cluster DE boxplots -------------------------------------------------------
    numClust <- sapply(cl[!grepl("^Comp",colnames(cl))],function(X) length(levels(X)))
    clustList <- reactive({ 
      temp <- as.list(colnames(d$cl)) 
      names(temp)[seq_along(numClust)] <- paste0(unlist(temp)[seq_along(numClust)],
                                                 ": ",numClust," clusters")
      if (length(temp) > length(numClust)) {
        names(temp)[setdiff(seq_along(temp),seq_along(numClust))] <- 
          paste0("Comparison: ",
                 sub("Comp.","",fixed=T,
                     x=unlist(temp)[setdiff(seq_along(temp),seq_along(numClust))]))  
      }
      return(temp)
    })
    output$resSelect <- renderUI({
      if (is.null(res())) { temp_sel <- savedRes} else { temp_sel <- res() }
      selectInput("res","Resolution:",choices=clustList(),selected=temp_sel)
    })
    output$saveButton <- renderUI({
      if (grepl("^Comp",input$res)) {
        actionButton("updateForViz","Save this comparison to disk",icon("save"))
      } else {
        actionButton("save","Save this resolution as default",icon("bookmark"))
      }
    })
    numClust <- numClust[numClust > 1]
    
    plot_cqPlot <- function() {
      numDEgenes <- lapply(get(input$deType)[!grepl("^Comp",names(get(input$deType)))],
                           function(X) sapply(X,nrow))
      toplim <- c(21,max(unlist(numDEgenes)) + 20)
      botlim <- c(-1,21)
      
      if (grepl("^Comp",input$res)) {
        par(mar=c(3,3.5,1,1))
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,paste("Press 'View clusters at this resolution'",
                         "to view the comparison",
                         sub("Comp.","",input$res,fixed=T),sep="\n"))
      } else {
        par(mar=c(0.2,3.5,1,1),mgp=2:0,mfrow=c(2,1))
        plot(x=numClust,y=sapply(numDEgenes,median),type="l",
             xlim=range(numClust)+c(-.5,.5),ylim=toplim,yaxs="i",xaxt="n",ylab=NA)
        abline(h=seq(0,max(unlist(numDEgenes)),10),lty=3,col=alpha(1,0.3))
        for (i in names(numDEgenes)[names(numDEgenes) != input$res]) {
          boxplot(numDEgenes[[i]],add=T,at=numClust[i],yaxt="n")
        }
        if (any(names(numDEgenes) == input$res)) {
          boxplot(numDEgenes[[input$res]],add=T,at=numClust[input$res],border="red")
        }
        
        par(mar=c(3,3.5,0.2,1),mgp=2:0)
        plot(x=numClust,y=sapply(numDEgenes,median),type="l",
             xlim=range(numClust)+c(-.5,.5),ylim=botlim,yaxs="i",xlab="Number of clusters",ylab=NA)
        abline(h=seq(0,max(unlist(numDEgenes)),10),lty=3,col=alpha(1,0.3))
        for (i in names(numDEgenes)[names(numDEgenes) != input$res]) {
          boxplot(numDEgenes[[i]],add=T,at=numClust[i],yaxt="n")
        }
        if (any(names(numDEgenes) == input$res)) {
          boxplot(numDEgenes[[input$res]],add=T,at=numClust[input$res],border="red")
        }
        mtext(switch(input$deType,
                     "deMarker"="Positive DE genes per cluster to all other clusters",
                     "deNeighb"="Positive DE genes per cluster to nearest cluster")
              ,side=2,line=2.5,at=botlim[2],xpd=NA)
      }
    }
    
    output$cqPlot <- renderPlot({
      print(plot_cqPlot())
    })
    
    output$cqPlotSave <- downloadHandler(
      filename="cqPlot.pdf",
      content=function(file) {
        pdf(file,width=7,height=6)
        print(plot_cqPlot())
        dev.off()
      }
    )
    
    # ^^ Silhouette plot -----------------------------------------------------------------
    plot_sil <- function() {
      tempSil <- cluster::silhouette(as.integer(d$cl[,input$res]),dist=silDist)
      par(mar=c(4.5,.5,1.5,1.5),mgp=2:0)
      if (length(tempSil) <= 1) {
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,paste("Silhouette plot cannot be computed",
                         "with less than two clusters.",sep="\n"))
      } else {
        plot(tempSil,beside=T,border=NA,main=NA,col=clustCols(input$res),do.n.k=T)
      }
    }
    
    output$sil <- renderPlot({
      print(plot_sil())
    })
    
    output$silSave <- downloadHandler(
      filename="sil.pdf",
      content=function(file) {
        pdf(file,width=6,height=7)
        print(plot_sil())
        dev.off()
      }
    )
    
    # ^^ Resolution selection buttons ----------------------------------------------------
    res <- reactiveVal() 
    observeEvent(input$go,res(input$res),ignoreNULL=F)
    observeEvent(input$go2,res(input$res2),ignoreNULL=F)
    
    observeEvent(input$save,{
      savedRes <<- input$res 
      # <<- updates variable outside scope of function. In this case, that's the enclosing
      # function (runShiny), where savedRes was set.
      save(savedRes,file=paste0(dataPath,dataTitle,"_savedRes.RData"))
    })
    
    
    # ^ Dataset and Cluster Metadata Inspection --------------------------------------------
    clusts <- reactive(d$cl[,res()])
    
    # ^^ Cell-type tSNE ####
    plot_tsne_labels <- function() {
      if (input$tsneLabels == "ca") {
        temp_labelNames <- sapply(unique(d$clusterID[[res()]]),function(X) 
          names(which(d$clusterID[[res()]] == X)),simplify=F)
        temp_labels <- apply(dr_viz,2,function(Y) 
          tapply(Y,apply(sapply(temp_labelNames,function(X) clusts() %in% X),1,which),mean))
        if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
        text(temp_labels,labels=names(temp_labelNames),font=2,cex=1.5)
      } else if (input$tsneLabels == "can") {
        temp_labels <- apply(dr_viz,2,function(X) tapply(X,clusts(),mean))
        if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
        text(temp_labels,labels=d$clusterID[[res()]],font=2,cex=1.5)
      } else if (input$tsneLabels == "cn") {
        temp_labels <- apply(dr_viz,2,function(X) tapply(X,clusts(),mean))
        if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
        text(temp_labels,labels=levels(clusts()),font=2,cex=1.5)
      } else {
        legend("center",legend="You changed the label choice names...")
      }
    }
    
    plot_tsne <- function() {
      par(mar=c(3,3,4,1),mgp=2:0)
      plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
           main=paste("tSNE at",res(),"using",ncol(dr_clust),"PCs"),
           xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
      if (any(ci())) {
        points(dr_viz[!ci(),],pch=21,
               col=alpha(clustCols(res())[clusts()],0.2)[!ci()],
               bg=alpha(clustCols(res())[clusts()],0.1)[!ci()])
        points(dr_viz[ci(),],pch=21,
               col=alpha(clustCols(res())[clusts()],1)[ci()],
               bg=alpha(clustCols(res())[clusts()],0.5)[ci()])
      } else {
        points(dr_viz,pch=21,
               col=alpha(clustCols(res())[clusts()],1),
               bg=alpha(clustCols(res())[clusts()],0.5))
      }
      if (hiC() != "") {
        mtext(side=3,line=-1,text=paste("Cluster",hiC(),"-",
                                        d$clusterID[[res()]][hiC()],"-",
                                        sum(clusts() == hiC()),"cells"))
      }
      if (input$nnArrow) {
        temp_nn <- sapply(deNeighb[[res()]],function(X) 
          unique(gsub(pattern="^vs\\.|\\.[A-Za-z]+?$","",colnames(X))),simplify=F)
        temp_labels <- apply(dr_viz,2,function(X) tapply(X,clusts(),mean))
        sapply(names(temp_nn),function(X)
          arrows(lwd=2,col=alpha("black",0.5),length=0.1,
                 x0=temp_labels[X,1],y0=temp_labels[X,2],
                 x1=temp_labels[temp_nn[[X]],1],y1=temp_labels[temp_nn[[X]],2]))
      }
    }
    
    output$tsne <- renderPlot({
      if (length(res()) > 0) {
        print(plot_tsne())
        print(plot_tsne_labels())
      }
    })
    
    output$tsneSave <- downloadHandler(
      filename="tsne.pdf",
      content=function(file) {
        pdf(file,width=7,height=7)
        print(plot_tsne())
        print(plot_tsne_labels())
        dev.off()
      }
    )
    
    # ^^ clusterSelect -------------------------------------------------------------------
    clusterSelect <- reactiveValues(cl=NULL)
    
    observeEvent(input$tsneClick,{ clusterSelect$cl <- input$tsneClick })
    
    cSelected <- reactive({
      t <- nearPoints(as.data.frame(dr_viz),clusterSelect$cl,xvar="tSNE_1",yvar="tSNE_2",threshold=5)
      t2 <- d$cl[rownames(t)[1],res()]
      if (is.na(t2)) { 
        return("") 
      } else if (t2 == "Unselected") {
        return("")
      } else { 
        return(t2) 
      }
    })
    
    hiC <- reactive({ 
      if (length(res()) < 1) {
        return("")
      } else if (input$genePlotClust != "") {
        d$cl[which(d$cl[,res()] == input$genePlotClust)[1],res()]
      } else {
        return(input$genePlotClust)
      }
    })
    
    ci <- reactive({
      if (hiC() == "") {
        rep(F,length(clusts()))
      } else {
        clusts() == hiC()
      }
    })
    
    # ^^ Metadata tSNE overlay -----------------------------------------------------------
    output$tsneMDlog <- renderUI({
      if (!(is.factor(md[,input$tsneMDcol]) | is.character(md[,input$tsneMDcol]))) {
        checkboxGroupInput("tsneMDlog",label="Colour scale",
                           choices=c("Log scale"="log"),width="100%")
      }
    })
    
    plot_tsneMD <- function() {
      if (is.factor(md[,input$tsneMDcol]) | is.character(md[,input$tsneMDcol])) {
        id <- as.factor(md[,input$tsneMDcol])
        if (length(levels(md[,input$tsneMDcol])) <= 8) {
          idcol <- RColorBrewer::brewer.pal(length(levels(md[,input$tsneMDcol])),
                              "Dark2")[1:length(levels(md[,input$tsneMDcol]))]
        } else {
          idcol <- rainbow2(length(levels(md[,input$tsneMDcol])))
        }
      } else {
        if ("log" %in% input$tsneMDlog) {
          id <- cut(log10(md[,input$tsneMDcol]),100)
        } else {
          id <- cut(md[,input$tsneMDcol],100)
        }
        idcol <- viridis(100,d=-1)
      }
      layout(cbind(2:1),heights=c(1,9))
      par(mar=c(3,3,0,1),mgp=2:0)
      plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
           xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
      if (any(ci())) {
        points(dr_viz[!ci(),],pch=21,
               col=alpha(idcol,.1)[id[!ci()]],
               bg=alpha(idcol,0.05)[id[!ci()]])
        points(dr_viz[ci(),],pch=21,
               col=alpha(idcol,.8)[id[ci()]],
               bg=alpha(idcol,0.4)[id[ci()]])
      } else {
        points(dr_viz,pch=21,
               col=alpha(idcol,.8)[id],
               bg=alpha(idcol,0.4)[id])
      }
      plot_tsne_labels()
      if (is.factor(md[,input$tsneMDcol]) | is.character(md[,input$tsneMDcol])) {
        par(mar=c(0,0,0,0))
        plot.new()
        legend("bottom",bty="n",horiz=T,pch=c(NA,rep(21,length(levels(md[,input$tsneMDcol])))),
               legend=c(paste0(input$tsneMDcol,":"),levels(md[,input$tsneMDcol])),
               col=c(NA,idcol),pt.bg=c(NA,alpha(idcol,0.5)))
      } else {
        if ("log" %in% input$tsneMDlog) {
          tempMain <- paste(input$tsneMDcol,"(log scale)")
        } else {
          tempMain <- input$tsneMDcol
        }
        par(mar=c(0,5,3,3))
        barplot(rep(1,100),space=0,col=idcol,xaxt="n",yaxt="n",border=NA,main=tempMain)
        text(x=c(1,100),y=1,pos=c(2,4),xpd=NA,labels=round(range(md[,input$tsneMDcol]),2))
      }
    }
    
    output$tsneMD <- renderPlot({
      if (length(res()) > 0) {
        print(plot_tsneMD())
      }
    })
    
    output$tsneMDSave <- downloadHandler(
      filename="tsneMD.pdf",
      content=function(file) {
        pdf(file,width=7,height=7)
        print(plot_tsneMD())
        dev.off()
      }
    )
    
    # ^^ Metadata Scatterplot ------------------------------------------------------------
    output$scatterLog <- renderUI({
      if ((is.factor(md[,input$mdScatterX]) | is.character(md[,input$mdScatterX])) |
          (is.factor(md[,input$mdScatterY]) | is.character(md[,input$mdScatterY]))) {
        checkboxGroupInput("scatterLog",inline=F,label=NULL,
                           choices=c("Log x axis"="x","Log y axis"="y","Show notch"="notch"),
                           selected="notch")
      } else {
        checkboxGroupInput("scatterLog",inline=F,label=NULL,
                           choices=c("Log x axis"="x","Log y axis"="y"))
      }
    })
    
    plot_mdScatter <- function() {
      if ((is.factor(md[,input$mdScatterX]) | is.character(md[,input$mdScatterX])) &
          (is.factor(md[,input$mdScatterY]) | is.character(md[,input$mdScatterY]))) {
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,"This figure is not designed to compare to categorical variables.")
      } else if (is.factor(md[,input$mdScatterX]) | is.character(md[,input$mdScatterX])) {
        par(mar=c(3,3,2,1),mgp=2:0)
        if (any(ci())) {
          temp1 <- tapply(md[!ci(),input$mdScatterY],as.factor(md[!ci(),input$mdScatterX]),c)
          temp2 <- tapply(md[ci(),input$mdScatterY],as.factor(md[ci(),input$mdScatterX]),c)
          plot(x=NULL,y=NULL,ylim=range(md[,input$mdScatterY]),
               xlim=c(0,length(levels(as.factor(md[,input$mdScatterX]))) * 3),
               log=sub("notch","",paste(input$scatterLog,collapse="")),xaxt="n",
               xlab=input$mdScatterX,ylab=input$mdScatterY)
          boxplot(temp1,add=T,xaxt="n",notch="notch" %in% input$scatterLog,
                  at=seq(1,length(levels(as.factor(md[,input$mdScatterX]))) * 3,by=3))
          boxplot(temp2,add=T,xaxt="n",notch="notch" %in% input$scatterLog,border="red",
                  at=seq(2,length(levels(as.factor(md[,input$mdScatterX]))) * 3,by=3))
          axis(side=1,at=seq(1.5,length(levels(as.factor(md[,input$mdScatterX]))) * 3,by=3),
               labels=names(temp1))
          legend("top",bty="n",xpd=NA,inset=c(0,-.05),pch=0,col="red",
                 legend=paste("Cluster",hiC(),"-",d$clusterID[[res()]][hiC()]))
        } else {
          boxplot(tapply(md[,input$mdScatterY],as.factor(md[,input$mdScatterX]),c),
                  xlab=input$mdScatterX,ylab=input$mdScatterY,
                  log=sub("notch","",paste(input$scatterLog,collapse="")),
                  notch="notch" %in% input$scatterLog)
        }
      } else if (is.factor(md[,input$mdScatterY]) | is.character(md[,input$mdScatterY])) {
        par(mar=c(3,3,2,1),mgp=2:0)
        if (any(ci())) {
          temp1 <- tapply(md[!ci(),input$mdScatterX],as.factor(md[!ci(),input$mdScatterY]),c)
          temp2 <- tapply(md[ci(),input$mdScatterX],as.factor(md[ci(),input$mdScatterY]),c)
          plot(x=NULL,y=NULL,xlim=range(md[,input$mdScatterX]),
               ylim=c(0,length(levels(as.factor(md[,input$mdScatterY]))) * 3),
               log=sub("notch","",paste(input$scatterLog,collapse="")),yaxt="n",
               xlab=input$mdScatterX,ylab=input$mdScatterY)
          boxplot(temp1,add=T,horizontal=T,yaxt="n",notch="notch" %in% input$scatterLog,
                  at=seq(1,length(levels(as.factor(md[,input$mdScatterY]))) * 3,by=3))
          boxplot(temp2,add=T,horizontal=T,yaxt="n",notch="notch" %in% input$scatterLog,border="red",
                  at=seq(2,length(levels(as.factor(md[,input$mdScatterY]))) * 3,by=3))
          axis(side=2,at=seq(1.5,length(levels(as.factor(md[,input$mdScatterY]))) * 3,by=3),
               labels=names(temp1))
          legend("top",bty="n",xpd=NA,inset=c(0,-.05),pch=0,col="red",
                 legend=paste("Cluster",hiC(),"-",d$clusterID[[res()]][hiC()]))
        } else {
          boxplot(tapply(md[,input$mdScatterX],as.factor(md[,input$mdScatterY]),c),
                  horizontal=T,xlab=input$mdScatterX,ylab=input$mdScatterY,
                  log=sub("notch","",paste(input$scatterLog,collapse="")),
                  notch="notch" %in% input$scatterLog)
        }
      } else {
        layout(matrix(c(2,1,0,3),2),c(5,1),c(1,5))
        par(mar=c(3,3,0,0),mgp=2:0,cex=1.1)
        plot(md[!ci(),input$mdScatterX],md[!ci(),input$mdScatterY],
             log=sub("notch","",paste(input$scatterLog,collapse="")),
             pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),
             xlab=input$mdScatterX,ylab=input$mdScatterY)
        points(md[ci(),input$mdScatterX],md[ci(),input$mdScatterY],
               pch=21,col=alpha("red",0.4),bg=alpha("red",0.2))
        if (any(ci())) {
          legend("topleft",bty="n",pch=21,col="red",pt.bg=alpha("red",0.5),
                 legend=paste("Cluster",hiC(),"-",d$clusterID[[res()]][hiC()]))
        }
        if ("x" %in% input$scatterLog) { tempLX <- "x" } else { tempLX <- "" }
        if ("y" %in% input$scatterLog) { tempLY <- "y" } else { tempLY <- "" }
        par(mar=c(0,3,1,0))
        boxplot(tapply(md[,input$mdScatterX],ci(),c),log=tempLX,
                horizontal=T,xaxt="n",yaxt="n",border=c("black","red"))
        par(mar=c(3,0,0,1))
        boxplot(tapply(md[,input$mdScatterY],ci(),c),log=tempLY,
                horizontal=F,xaxt="n",yaxt="n",border=c("black","red"))
      }
    }
    
    output$mdScatter <- renderPlot({
      if (length(res()) > 0) {
        print(plot_mdScatter())
      }
    })
    
    output$mdScatterSave <- downloadHandler(
      filename="mdScatter.pdf",
      content=function(file) {
        pdf(file,width=7,height=7)
        print(plot_mdScatter())
        dev.off()
      }
    )
    
    # ^^ Metadata Factor Barplot ---------------------------------------------------------
    output$mdFactorOpts <- renderUI({
      if (is.factor(md[,input$mdFactorData]) | is.character(md[,input$mdFactorData])) {
        radioButtons("mdFactorRA","Factor counts per cluster:",inline=T,
                     choices=list("Absolute"="absolute","Relative"="relative"))
      } else {
        checkboxGroupInput("mdFactorOpts",inline=T,label="Figure options",
                           choices=c("Log scale"="y","Show notch"="notch"),selected="notch")
      }
    })
    
    plot_mdFactor <- function() {
      if (is.factor(md[,input$mdFactorData]) | is.character(md[,input$mdFactorData])) {
        id <- switch(input$mdFactorRA,
                     "relative"=tapply(md[,input$mdFactorData],clusts(),
                                       function(X) table(X) / length(X)),
                     "absolute"=tapply(md[,input$mdFactorData],clusts(),table))
        if (is.list(id)) { id <- do.call(cbind,id) }
        idylab <- switch(input$mdFactorRA,
                         "relative"="Proportion of cells per cluster",
                         "absolute"="Number of cells per cluster")
        if (length(levels(md[,input$mdFactorData])) <= 8) {
          idcol <- RColorBrewer::brewer.pal(length(levels(md[,input$mdFactorData])),
                              "Dark2")[1:length(levels(md[,input$mdFactorData]))]
        } else {
          idcol <- rainbow2(length(levels(md[,input$mdFactorData])))
        }
        par(mar=c(3,3,2,1),mgp=2:0)
        barplot(id,col=idcol,ylab=idylab,
                legend.text=levels(md[,input$mdFactorData]),
                args.legend=list(x="topright",horiz=T,inset=c(0,-.08),bty="n"))
        mtext(input$mdFactorData,side=3,adj=0,font=2,line=1,cex=1.2)
      } else {
        par(mar=c(3,3,2,1),mgp=2:0)
        boxplot(tapply(md[,input$mdFactorData],cl[,res()],c),
                ylab=input$mdFactorData,notch="notch" %in% input$mdFactorOpts,
                log=sub("notch","",paste(input$mdFactorOpts,collapse="")),
                border=clustCols(res()),col=alpha(clustCols(res()),0.3))
      }
    }
    
    output$mdFactor <- renderPlot({
      if (length(res()) > 0) {
        print(plot_mdFactor())
      }
    })
    
    output$mdFactorSave <- downloadHandler(
      filename="mdFactor.pdf",
      content=function(file) {
        pdf(file,width=7,height=7)
        print(plot_mdFactor())
        dev.off()
      }
    )
    
    # ^ Differentially Expressed Genes per Cluster -----------------------------------------
    
    # ^^ Dotplot of differentially expressed genes per cluster ---------------------------------------
    
    output$heatDEtype <- renderUI({
      if (!is.null(res())) {
        if (grepl("^Comp",res())) {
          temp <- list("DE vs rest"="deTissue",
                       "Set A vs Set B"="deMarker")
        } else {
          temp <- list("DE vs rest"="deTissue",
                       "Marker genes"="deMarker",
                       "DE vs neighbour"="deNeighb")
        }
        radioButtons("heatG","Dotplot Genes:",choices=temp,selected="deMarker")
      }
    })
    
    output$DEgeneSlider <- renderUI({
      if (length(res()) > 0) {
        if (input$heatG == "deTissue") {
          sliderInput("DEgeneCount",min=1,max=max(sapply(d$deTissue[[res()]],nrow)),
                      value=5,step=1,ticks=T,width="100%",
                      label=HTML(paste(
                        "Positive differential gene expression of cluster over tissue",
                        "# of genes per cluster to show",sep="<br/>")))
        } else if (input$heatG == "deMarker") {
          if (grepl("^Comp",res())) {
            temp_label <- HTML(paste(
              "Positive differential gene expression between sets",
              "# of genes per set to show",sep="<br/>"))
          } else {
            temp_label <- HTML(paste(
              "Positive differential gene expression between cluster and all other clusters",
              "# of genes per cluster to show",sep="<br/>"))
          }
          sliderInput("DEgeneCount",min=1,max=max(sapply(d$deMarker[[res()]],nrow)),
                      value=5,step=1,ticks=T,width="100%",
                      label=temp_label)
        } else if (input$heatG == "deNeighb") {
          sliderInput(
            "DEgeneCount",min=1,max=max(sapply(deNeighb[[res()]],nrow)),
            value=5,step=1,ticks=T,width="100%",
            label=HTML(paste(
              "Positive differential gene expression between cluster and nearest neighbour",
              "# of genes per cluster to show",sep="<br/>"))
            )
        }
      }
    })
    
    output$DEclustSelect <- renderUI({
      if (length(res()) > 0) {
        selectInput("DEclustNum","Cluster # for gene list",
                    choices=levels(clusts())[!levels(clusts()) == "Unselected"])
      }
    })
    
    heatGenes <- reactive({
      temp <- unique(unlist(lapply(
        switch(input$heatG,
               deTissue=d$deTissue[[res()]],
               deMarker=d$deMarker[[res()]],
               deNeighb=deNeighb[[res()]]),
        function(X) 
          if (nrow(X) == 0) { NA } else { rownames(X)[1:input$DEgeneCount] }
      )))
      temp <- temp[!is.na(temp)]
      return(temp)
    })
    
    temp_DR <- reactive({
      sapply(d$CGS[[res()]],function(X) X[heatGenes(),"DR"])
    })
    temp_MDTC <- reactive({
      sapply(d$CGS[[res()]],function(X) X[heatGenes(),"MDTC"])
    })
    
    hC <- reactive({ 
      if (exists("deDist")) {
        if (res() %in% names(deDist)) {
          return(hclust(as.dist(deDist[[res()]]),"single"))
        } else {
          return(hclust(dist(t(temp_DR())),"single"))
        }
      } else {
        return(hclust(dist(t(temp_DR())),"single"))
      }
    })
    hG <- reactive(hclust(dist(temp_DR()),"complete"))
    
    plot_dotplot <- function() {
      if (length(levels(clusts())) <= 1) {
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,paste("Heatmap cannot be computed",
                         "with less than two clusters.",sep="\n"))
      } else if (length(heatGenes()) < 1) {
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,"There are no differentially expressed genes.")
      } else {
        if ("Unselected" %in% levels(clusts())) {
          tempLabRow <- c(paste(levels(clusts())[!levels(clusts()) == "Unselected"],
                                paste(sapply(switch(input$heatG,
                                                    deTissue=d$deTissue[[res()]],
                                                    deMarker=d$deMarker[[res()]],
                                                    deNeighb=deNeighb[[res()]]),nrow),"DE"),
                                sep=": "),"Unselected")
        } else {
          tempLabRow <- paste(paste0("Cluster ",levels(clusts())),
                              paste(sapply(switch(input$heatG,
                                                  deTissue=d$deTissue[[res()]],
                                                  deMarker=d$deMarker[[res()]],
                                                  deNeighb=deNeighb[[res()]]),nrow),"DE"),
                              sep=": ")
        }
        tempLabCol <- heatGenes()
        if (exists("symbolMap")) {
          temp <- symbolMap[tempLabCol]
          tempLabCol[!is.na(temp)] <- temp[!is.na(temp)]
        }
        DR <- temp_DR()[hG()$order,hC()$order]
        temp <- range(sapply(d$CGS[[res()]],function(X) X[,"MDTC"]))
        temp <- seq(temp[1],temp[2],length.out=101)
        MDTC <- findInterval(as.vector(temp_MDTC()[hG()$order,hC()$order]),
                             vec=temp,all.inside=T)

        layout(matrix(c(0,2,3,1),2),widths=c(1,11),heights=c(1,5))
        par(mar=c(9,0,0,8))
        plot(x=NULL,y=NULL,xlim=c(0.5,nrow(DR)+.5),ylim=c(0.5,ncol(DR)+.5),
             xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n")
        abline(v=1:nrow(DR),col="grey90")
        symbols(x=rep(1:nrow(DR),ncol(DR)),
                y=as.vector(sapply(1:ncol(DR),function(X) rep(X,nrow(DR)))),
                circles=as.vector(DR)/2,inches=F,add=T,
                fg=viridis::viridis(100,d=-1)[MDTC],
                bg=viridis::viridis(100,d=-1)[MDTC])
        axis(side=1,at=1:nrow(DR),lwd=0,labels=tempLabCol[hG()$order],las=2)
        axis(side=4,at=1:ncol(DR),lwd=0,labels=tempLabRow[hC()$order],las=1)
        
        par(mar=c(9,0,0,0))
        plot(as.dendrogram(hC()),leaflab="none",horiz=T,
             ylim=c(0.5,length(hC()$order)+.5),yaxs="i",yaxt="n")
        
        par(mar=c(0,0,0,8))
        plot(as.dendrogram(hG()),leaflab="none",
             xlim=c(0.5,length(hG()$order)+.5),xaxs="i",yaxt="n")
      }
    }
    
    output$dotplot <- renderPlot({
      if (length(res()) > 0) {
        print(plot_dotplot())
      }
    })
    
    output$heatmapSave <- downloadHandler(
      filename="heatmap.pdf",
      content=function(file) {
        pdf(file,width=10,height=5)
        print(plot_dotplot())
        dev.off()
      }
    )
    
    output$deGeneSave <- downloadHandler(
      filename=function() { paste0(input$heatG,"_",input$DEclustNum,".txt") },
      content=function(file) {
        outTable <- switch(input$heatG,
                           deTissue=d$deTissue[[res()]][[input$DEclustNum]],
                           deMarker=d$deMarker[[res()]][[input$DEclustNum]],
                           deNeighb=deNeighb[[res()]][[input$DEclustNum]])
        write.table(outTable,file,quote=F,sep="\t",row.names=T,col.names=NA)
      }
    )
    
    # ^ Gene Expression Distributions per Cluster ------------------------------------------
    # ^^ Gene search box -----------------------------------------------------------------
    output$geneSearchBox <- renderUI({
      switch(input$searchType,
             comma=textInput("GOI",width="100%",
                             label=paste("Enter list of genes,",
                                         "(comma/space-separated, case-insensitive)",
                                         "and click Search")),
             regex=textInput("GOI",value="^ACTB$",width="100%",
                             label="Search for genes by regular expression and click Search"))
    })
    
    geneSearch <- function(txt,st) {
      if (exists("symbolMap")) {
        temp <- switch(st,
                       comma={
                         temp_in <- strsplit(txt,split="[\\s,]",perl=T)[[1]]
                         temp_out <- names(symbolMap)[toupper(symbolMap) %in% toupper(temp_in)]
                         names(temp_out) <- symbolMap[toupper(symbolMap) %in% toupper(temp_in)]
                         temp_out
                       },
                       regex={
                         temp_in <- grep(txt,symbolMap,ignore.case=T)
                         temp_out <- names(symbolMap)[temp_in]
                         names(temp_out) <- symbolMap[temp_in]
                         temp_out
                       })
        if (length(temp) > 0) {
          return(temp)
        } else {
          return(switch(st,
                        comma={
                          temp_in <- strsplit(txt,split="[\\s,]",perl=T)[[1]]
                          rownames(nge)[which(toupper(rownames(nge)) %in% toupper(temp_in))]
                        },
                        regex=grep(txt,rownames(nge),value=T,ignore.case=T)))
        }
        
      } else {
        return(switch(st,
                      comma={
                        temp_in <- strsplit(txt,split="[\\s,]",perl=T)[[1]]
                        rownames(nge)[which(toupper(rownames(nge)) %in% toupper(temp_in))]
                      },
                      regex=grep(txt,rownames(nge),value=T,ignore.case=T)))
      }
    }
    
    GOI <- eventReactive(input$GOIgo,geneSearch(input$GOI,input$searchType),ignoreNULL=F)
    
    # ^^ Scatterplot of gene expression in cluster ------------------------------------------------------
    output$genePlotClustSelect <- renderUI({
      if (length(res()) > 0) {
        selectInput("genePlotClust","Cluster:",selected=cSelected(),
                    choices=c("",levels(clusts())[!levels(clusts()) == "Unselected"]))
      }
    })
    
    cellMarkCols <- reactive(rainbow2(length(cellMarkers)))
    
    plot_clusterGenes <- function() {
      doubleDot <- function(col1,col2) {
        upper.half.circle <- function(col1){  
          rs <- seq(0,pi,len=100) + pi/2
          xc <- 0+cos(rs) 
          yc <- 0+sin(rs) 
          polygon(xc,yc,col=col1,border=NA)
        } 
        lower.half.circle <- function(col2){ 
          rs <- seq(0,pi,len=100) + pi/2
          xc <- 0-cos(rs) 
          yc <- 0-sin(rs) 
          polygon(xc,yc,col=col2,border=NA)
        }
        upper.half.circle(col1)
        lower.half.circle(col2)
        rs <- seq(0,2*pi,len=200)
        polygon(cos(rs),sin(rs),border="white")
      }
      singleDot <- function(col1){  
        rs <- seq(0,2*pi,len=200)
        xc <- 0+cos(rs) 
        yc <- 0+sin(rs) 
        polygon(xc,yc,col=col1,border="white")
      } 
      par(mar=c(3,3,3,20),mgp=2:0)
      if (hiC() == "") {
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,paste("Click a cell from a cluster on the tSNE plot above",
                         "or select a cluster from the drop-down list above left",
                         "to see gene expression for that cluster.",sep="\n"))
      } else {
        temp_ylab <- switch(as.character(exponent == exp(1)),
                            "TRUE"="(natural log scale)",
                            "FALSE"=paste0("(log",exponent," scale)"))
        plot(MDTC~DR,
             data=d$CGS[[res()]][[hiC()]][
               !((d$CGS[[res()]][[hiC()]]$cMu | d$CGS[[res()]][[hiC()]]$cMs) & 
                   d$CGS[[res()]][[hiC()]]$overCut),],
             col=alpha("black",0.3),
             xlab="Proportion of cells in which gene was detected",
             ylab=paste("Mean normalized gene expression where detected",temp_ylab))
        title(paste0("Cluster ", hiC(),": ",d$clusterID[[res()]][hiC()]),cex=1.2)
        mtext(paste("Cells:",sum(clusts()==hiC()),
                    "   Genes detected:",length(d$CGS[[res()]][[hiC()]]$DR)),side=3,line=0,cex=0.9)
        box(col=clustCols(res())[hiC()],lwd=2)
        
        if (input$cgLegend == "markers") {
          for (x in which(d$CGS[[res()]][[hiC()]]$cMu)) {
            TeachingDemos::my.symbols(x=d$CGS[[res()]][[hiC()]]$DR[x],
                       y=d$CGS[[res()]][[hiC()]]$MDTC[x],
                       symb=singleDot,inches=0.1,
                       MoreArgs=list(col1=cellMarkCols()[which(sapply(cellMarkersU,function(X) 
                         d$CGS[[res()]][[hiC()]]$genes[x] %in% X))]))
          }
          for (x in which(d$CGS[[res()]][[hiC()]]$cMs)) {
            temp <- unlist(strsplit(names(which(sapply(cellMarkersS,function(X) 
              d$CGS[[res()]][[hiC()]]$genes[x] %in% X))),"&"))
            TeachingDemos::my.symbols(x=d$CGS[[res()]][[hiC()]]$DR[x],
                       y=d$CGS[[res()]][[hiC()]]$MDTC[x],
                       symb=doubleDot,inches=0.1,
                       MoreArgs=list(col1=cellMarkCols()[as.integer(temp[1])],
                                     col2=cellMarkCols()[as.integer(temp[2])]))
          }
          for (x in which(d$CGS[[res()]][[hiC()]]$cMu & d$CGS[[res()]][[hiC()]]$overCut)) {
            text(x=d$CGS[[res()]][[hiC()]]$DR[x],y=d$CGS[[res()]][[hiC()]]$MDTC[x],
                 labels=d$CGS[[res()]][[hiC()]]$genes[x],srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
                 col=cellMarkCols()[which(sapply(cellMarkersU,function(X) 
                   d$CGS[[res()]][[hiC()]]$genes[x] %in% X))])
          }
          for (x in which(d$CGS[[res()]][[hiC()]]$cMs & d$CGS[[res()]][[hiC()]]$overCut)) {
            text(x=d$CGS[[res()]][[hiC()]]$DR[x],y=d$CGS[[res()]][[hiC()]]$MDTC[x],
                 labels=d$CGS[[res()]][[hiC()]]$genes[x],srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
                 col=cellMarkCols()[as.integer(temp[2])])
          }
          legend(x=1.05,y=max(d$CGS[[res()]][[hiC()]]$MDTC),xpd=NA,bty="n",ncol=1,
                 pch=19,col=cellMarkCols(),legend=names(cellMarkersU))
          
        } else if (input$cgLegend == "heatmap") {
          degl <- rownames(d$CGS[[res()]][[hiC()]]) %in% 
            rownames(switch(input$heatG,
                            deTissue=d$deTissue[[res()]],
                            deMarker=d$deMarker[[res()]],
                            deNeighb=deNeighb[[res()]])[[hiC()]])[1:input$DEgeneCount]
          if (any(degl)) {
            points(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
                   pch=16,cex=1.2,col="darkred")
            text(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
                 srt=315,cex=1.5,font=2,adj=c(1.1,-.1),col="darkred",
                 labels=d$CGS[[res()]][[hiC()]]$genes[degl])
          }
          temp_n <- nrow(switch(input$heatG,
                                deTissue=d$deTissue,
                                deMarker=d$deMarker,
                                deNeighb=deNeighb)[[res()]][[hiC()]])
          temp_lab <- switch(input$heatG,
                             deTissue=" DE genes vs rest of cells in sample",
                             deMarker=" marker genes",
                             deNeighb=" DE genes vs nearest neighbouring cluster")
          legend("top",bty="n",pch=16,col="darkred",
                 legend=paste0(temp_n,temp_lab," (showing top ",
                               min(temp_n,input$DEgeneCount),")"))
        } else if (input$cgLegend == "search" & length(GOI()) > 0) {
          degl <- which(rownames(nge) %in% GOI())
          points(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
                 pch=16,cex=1.2,col="darkred")
          text(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
               srt=315,cex=1.5,font=2,adj=c(1.1,-.1),col="darkred",
               labels=d$CGS[[res()]][[hiC()]]$genes[degl])
        }
      }
    }
    
    clickGenes <- reactiveVal() 
    observeEvent(input$cgClick,{
      t <- nearPoints(d$CGS[[res()]][[hiC()]],input$cgClick,xvar="DR",yvar="MDTC")
      temp_out <- rownames(t)
      names(temp_out) <- t$genes
      clickGenes(temp_out)
    })
    observeEvent(input$scatterClick,{
      t <- nearPoints(compDF(),input$scatterClick,xvar="x",yvar="y")
      temp_out <- rownames(t)
      names(temp_out) <- t$genes
      clickGenes(temp_out)
    })
    
    output$clusterGenes <- renderPlot({
      if (length(res()) > 0) {
        print(plot_clusterGenes())
      }
    })
    
    output$clusterGenesSave <- downloadHandler(
      filename="clusterGenes.pdf",
      content=function(file) {
        pdf(file,width=12,height=7)
        print(plot_clusterGenes())
        dev.off()
      }
    )
    
    # ^ Gene expression comparison -------------------------------------------------------
    # ^^ Gene selection by click/search --------------------------------------------------
    output$cgSelect <- renderUI({
      if (length(res()) > 0) {
          temp_choices <- switch(input$boxplotGene,
                                 click=clickGenes(),
                                 search=GOI())
          if (is.null(names(temp_choices))) {
            temp_choices <- sort(temp_choices)
          } else {
            temp_choices <- temp_choices[order(names(temp_choices))]
          }
        selectInput("cgGene",choices=temp_choices,label="Select gene from list:")
      }
    })
    
    # ^^ Boxplots for gene expression comparison ------------------------------------------------------
    plot_geneTest <- function() {
      if (input$cgGene == "") {
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,paste("Select a gene by either clicking on the plot above",
                         "or searching for genes of interest in the search bar above,",
                         "then pick the gene from the list just above this figure",
                         "to see a comparison of that gene's expression across all clusters.",
                         sep="\n"))
      } else {
        temp_ylab <- switch(as.character(exponent == exp(1)),
                            "TRUE"="(natural log scale)",
                            "FALSE"=paste0("(log",exponent," scale)"))
        temp_pos <- switch(as.character(length(levels(clusts())) > 1),"TRUE"=hC()$order,"FALSE"=1)
        layout(matrix(2:1,nrow=2),heights=c(1,4))
        par(mar=c(3,3,0,3),mgp=2:0)
        suppressWarnings(boxplot(vector("list",length(levels(clusts()))),
                                 ylim=range(nge[input$cgGene,]),
                                 ylab=paste(input$cgGene,"normalized gene expression",temp_ylab),
                                 xlab=NA,xaxt="n"))
        mtext(levels(clusts())[temp_pos],side=1,line=0,at=seq_along(temp_pos))
        mtext("Clusters, ordered by heatmap dendrogram",side=1,line=1)
        try(tempGeneName <- mapIds(annotationDB,keys=input$cgGene,keytype=rownameKeytype,
                                   column="GENENAME",multiVals="first"),silent=T)
        if (exists("tempGeneName")) { 
          mtext(paste(paste0(names(tempGeneName),": ",tempGeneName),collapse="\n"),
                side=1,line=2,font=2) 
        } else {
          mtext(paste(input$cgGene,collapse="\n"),side=1,line=2,font=2) 
        }
        if ("sct" %in% input$bxpOpts) {
          bxpCol <- alpha(clustCols(res()),.2)
        } else {
          bxpCol <- alpha(clustCols(res()),.8)
        }
        for (i in temp_pos) {
          boxplot(nge[input$cgGene,clusts() == levels(clusts())[i]],add=T,
                  at=which(temp_pos == i),notch="notch" %in% input$bxpOpts,col=bxpCol[i],outline=F)
          if ("sct" %in% input$bxpOpts) {
            points(jitter(rep(which(temp_pos == i),sum(clusts() == levels(clusts())[i])),amount=.2),
                   nge[input$cgGene,clusts() == levels(clusts())[i]],
                   pch=20,col=alpha(clustCols(res())[i],.4))
          }
        }
        if ("rnk" %in% input$bxpOpts) {
          points(x=seq_along(d$CGS[[res()]]),
                 y=sapply(d$CGS[[res()]][temp_pos],function(X) X[input$cgGene,"MTCrank"]) * 
                   max(nge[input$cgGene,]) + min(nge[input$cgGene,]),
                 pch=25,cex=1.2,col="darkred",bg="firebrick2")
          axis(side=4,at=seq(0,1,.25) * max(nge[input$cgGene,]) + min(nge[input$cgGene,]),
               labels=paste0(seq(0,1,.25) * 100,"%"),col.ticks="darkred",col.axis="darkred")
          mtext(side=4,line=2,text="Quantile of gene expression per cluster",col="darkred")
        }
        if (length(temp_pos) > 1) { 
          par(new=F,mar=c(0,3,1,3))
          plot(as.dendrogram(hC()),leaflab="none") 
        }
      }
    }
    
    output$geneTest <- renderPlot({
      if (length(res()) > 0) {
        print(plot_geneTest())
      }
    })
    
    output$geneTestSave <- downloadHandler(
      filename="geneTest.pdf",
      content=function(file) {
        pdf(file,width=12,height=7)
        print(plot_geneTest())
        dev.off()
      }
    )
    
    
    # ^ Distribution of genes of interest ------------------------------------------------
    output$geneSearchBox1 <- renderUI({
      if (input$searchType1 == "comma") {
        textInput("GOI1",width="100%",
                  label=paste("Enter list of genes"))
      } else if (input$searchType1 == "regex") {
        textInput("GOI1",value="^ACTB$",width="100%",
                  label="Enter regular expression")
      }
    })
    GOI1 <- eventReactive(input$GOI1go,geneSearch(input$GOI1,input$searchType1))
    output$GOI1select <- renderUI({ 
      temp_choices <- GOI1()
      if (is.null(names(temp_choices))) {
        temp_choices <- sort(temp_choices)
      } else {
        temp_choices <- temp_choices[order(names(temp_choices))]
      }
      selectInput("goi1",label="Select genes:",choices=temp_choices,multiple=T)
    })
    
    output$geneSearchBox2 <- renderUI({
      if (input$searchType2 == "comma") {
        textInput("GOI2",width="100%",
                  label=paste("Search by list of genes"))
      } else if (input$searchType2 == "regex") {
        textInput("GOI2",value="^ACTB$",width="100%",
                  label="Search by regular expression")
      }
    })
    GOI2 <- eventReactive(input$GOI2go,geneSearch(input$GOI2,input$searchType2))
    output$GOI2select <- renderUI({ 
      temp_choices <- GOI2()
      if (is.null(names(temp_choices))) {
        temp_choices <- sort(temp_choices)
      } else {
        temp_choices <- temp_choices[order(names(temp_choices))]
      }
      selectInput("goi2",label="Select genes:",choices=temp_choices,multiple=T)
    })
    
    plot_tsneClust <- function() {
      par(mar=c(3,3,4,1),mgp=2:0)
      plot(dr_viz,pch=21,
           col=alpha(clustCols(res())[clusts()],1),
           bg=alpha(clustCols(res())[clusts()],0.5),
           xlab="tSNE_1",ylab="tSNE_2",
           main=paste("tSNE at",res(),"using",ncol(dr_clust),"PCs"))
    }
    
    plot_goi <- function(goi) {
      if (length(goi) < 1) {
        plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
        text(.5,.5,paste("To search for your gene(s) of interest type a",
                         "list of genes or regex in the box above", 
                         "then select the gene(s) from the drop-down list",
                         "in the \"Select genes:\" box above right.",sep="\n"))
      } else {
        if (length(goi) > 5) { goiL <- 5 } else { goiL <- length(goi) }
        if (goiL > 1) {
          gv <- apply(nge[goi,],2,max)
        } else {
          gv <- nge[goi,]
        }
        cv <- cut(gv,breaks=100,labels=F)
        par(mar=c(3,3,goiL+1,1),mgp=2:0)
        plot(dr_viz,pch=21,cex=1.3,xlab="tSNE_1",ylab="tSNE_2",
             col=viridis(100,.7,d=-1)[cv],bg=viridis(100,.3,d=-1)[cv])
        temp_yrange <- max(dr_viz[,2]) - min(dr_viz[,2])
        segments(x0=seq(quantile(range(dr_viz[,1]),.7),
                        quantile(range(dr_viz[,1]),1),length.out=1000),
                 y0=max(dr_viz[,2]) + temp_yrange * .045,
                 y1=max(dr_viz[,2]) + temp_yrange * .065,
                 col=viridis(1000,d=-1),xpd=NA)
        text(x=c(quantile(range(dr_viz[,1]),.7),
                 quantile(range(dr_viz[,1]),.85),
                 quantile(range(dr_viz[,1]),1)),
             y=rep(max(dr_viz[,2]) + temp_yrange * .06,3),
             labels=c(round(min(gv),2),"Max expression per cell",round(max(gv),2)),pos=2:4,xpd=NA)
        try(tempGeneName <- mapIds(annotationDB,keys=goi,keytype=rownameKeytype,
                                   column="GENENAME",multiVals="first"),silent=T)
        if (exists("tempGeneName")) { 
          tempMissing <- is.na(tempGeneName)
          tempGeneName[tempMissing] <- names(tempGeneName)[tempMissing]
          tempGeneName[!tempMissing] <- paste0(names(tempGeneName),": ",tempGeneName)[!tempMissing]
          if (length(tempGeneName) > 4) { 
            tempGeneName[5] <- "and more..."; tempGeneName <- tempGeneName[1:5] 
          }
          title(paste(tempGeneName,collapse="\n"),line=0.25,adj=0,font.main=1)
        } else {
          temp_goi <- goi
          if (length(temp_goi) > 4) { 
            temp_goi[5] <- "and more..."; temp_goi <- temp_goi[1:5] 
          }
          title(paste(temp_goi,collapse="\n"),line=0.25,adj=0,font.main=1)
        }
      }
    }
    
    output$goiPlot1 <- renderPlot({
      if (input$plotClust1 == "clust" & length(res()) > 0) {
        print(plot_tsneClust())
        if (input$plotLabel1) { print(plot_tsne_labels()) }
      } else if (input$plotClust1 == "goi") {
        print(plot_goi(input$goi1))
        if (input$plotLabel1 & length(res()) > 0 & length(input$goi1) > 0) {
          print(plot_tsne_labels())
        }
      }
    })
    
    output$goiPlot1Save <- downloadHandler(
      filename="goi1.pdf",
      content=function(file) {
        pdf(file,width=7,height=7)
        if (input$plotClust1 == "clust" & length(res()) > 0) {
          print(plot_tsneClust())
          if (input$plotLabel1) { print(plot_tsne_labels()) }
        } else if (input$plotClust1 == "goi") {
          print(plot_goi(input$goi1))
          if (input$plotLabel1 & length(res()) > 0 & length(input$goi1) > 0) {
            print(plot_tsne_labels())
          }
        }
        dev.off()
      }
    )
    
    output$goiPlot2 <- renderPlot({
      if (input$plotClust2 == "clust" & length(res()) > 0) {
        print(plot_tsneClust())
        if (input$plotLabel2) { print(plot_tsne_labels()) }
      } else if (input$plotClust2 == "goi") {
        print(plot_goi(input$goi2))
        if (input$plotLabel2 & length(res()) > 0 & length(input$goi2) > 0) {
          print(plot_tsne_labels())
        }
      }
    })
    
    output$goiPlot2Save <- downloadHandler(
      filename="goi2.pdf",
      content=function(file) {
        pdf(file,width=7,height=7)
        if (input$plotClust2 == "clust" & length(res()) > 0) {
          print(plot_tsneClust())
          if (input$plotLabel2) { print(plot_tsne_labels()) }
        } else if (input$plotClust2 == "goi") {
          print(plot_goi(input$goi2))
          if (input$plotLabel2 & length(res()) > 0 & length(input$goi2) > 0) {
            print(plot_tsne_labels())
          }
        }
        dev.off()
      }
    )

    
    # ^ MA plot for cluster comparison ---------------------------------------------------------------
    output$resSelect2 <- renderUI({
      selectInput("res2","Resolution:",choices=clustList(),selected=res(),width="100%")
    })
    output$saveButton2 <- renderUI({
      if (grepl("^Comp",input$res2)) {
        actionButton("updateForViz2","Save this comparison to disk",icon("save"))
      } 
    })
    output$setScatterY <- renderUI({
      if ("Unselected" %in% levels(clusts())) {
        selectInput("ssY",label="Cluster A (A-B comparison)",selected="Set A",
                    choices=levels(clusts())[!levels(clusts()) == "Unselected"])
      } else {
        selectInput("ssY",label="Cluster A (A-B comparison)",choices=c("",levels(clusts())),selected=hiC())
      }
    })
    output$setScatterX <- renderUI({
      if ("Unselected" %in% levels(clusts())) {
        selectInput("ssX",label="Cluster B (A-B comparison)",selected="Set B",
                    choices=levels(clusts())[!levels(clusts()) == "Unselected"])
      } else {
        selectInput("ssX",label="Cluster B (A-B comparison)",choices=c("",levels(clusts())),
                    selected=unique(gsub(pattern="^vs\\.|\\.[A-Za-z]+?$","",
                                         colnames(deNeighb[[res()]][[hiC()]]))))
      }
    })
    output$diffLabelSelect <- renderUI({
      if (input$diffLabelType == "diff") {
        sliderInput("diffCount",min=1,max=100,value=5,step=1,width="100%",
                    label="Number of genes to label")
      } else if (input$diffLabelType == "de") {
        if (input$heatG == "deTissue") {
          sliderInput("diffCount",value=5,step=1,ticks=T,width="100%",
                      min=1,max=max(sapply(d$deTissue[[res()]][c(input$ssX,input$ssY)],nrow)),
                      label="DE vs rest: # of genes to label")
        } else if (input$heatG == "deMarker") {
          if (grepl("^Comp",res())) {
            temp_label <- "Set A vs Set B: # of genes to label"
          } else {
            temp_label <- "Marker genes: # of genes to label"
          }
          sliderInput("diffCount",
                      min=1,max=max(sapply(d$deMarker[[res()]][c(input$ssX,input$ssY)],nrow)),
                      value=5,step=1,ticks=T,width="100%",
                      label=temp_label)
        } else if (input$heatG == "deNeighb") {
          sliderInput("diffCount",value=5,step=1,ticks=T,width="100%",
                      min=1,max=max(sapply(deNeighb[[res()]][c(input$ssX,input$ssY)],nrow)),
                      label="DE vs neighbour: # of genes to label")
        }
      }
    })
    
    compDF <- reactive({
      data.frame(y=d$CGS[[res()]][[input$ssY]][,input$scatterInput] - 
                   d$CGS[[res()]][[input$ssX]][,input$scatterInput],
                 x=rowMeans(cbind(d$CGS[[res()]][[input$ssX]][,input$scatterInput],
                                  d$CGS[[res()]][[input$ssY]][,input$scatterInput])),
                 genes=d$CGS[[res()]][[input$ssX]]$genes,
                 row.names=rownames(d$CGS[[res()]][[input$ssX]]))
    })
    
    LBF <- reactive({
      lm(y~x,data=compDF())
    })
    
    diffRanked <- reactive({
      if (input$scatterLine == "sub") {
        temp <- d$CGS[[res()]][[input$ssY]][,input$scatterInput] - 
          d$CGS[[res()]][[input$ssX]][,input$scatterInput]
        names(temp) <- rownames(d$CGS[[res()]][[input$ssY]])
        return(sort(temp,decreasing=T))
      } else if (input$scatterLine == "lbf") {
        temp <- LBF()$residuals
        names(temp) <- rownames(d$CGS[[res()]][[input$ssY]])
        return(sort(temp,decreasing=T))
      }
    })
    
    plot_setScatter <- function() {
      if (!is.null(res())) {
        if (input$ssX %in% levels(clusts()) & input$ssY %in% levels(clusts())) {
          temp_exp <- switch(as.character(exponent == exp(1)),
                             "TRUE"="(natural log scale)",
                             "FALSE"=paste0("(log",exponent," scale)"))
          temp_label <- switch(input$scatterInput,
                               "MTC"=paste("mean normalized gene expression",temp_exp),
                               "MDTC"=paste("mean normalized gene expression where detected",temp_exp),
                               "DR"="proportion of cells in which gene was detected")
          par(mar=c(3,3,2,1),mgp=2:0)
          plot(y~x,data=compDF(),
               ylab=paste0("Difference in ",temp_label," (",input$ssY," - ",input$ssX,")"),
               xlab=paste0("Average of ",temp_label," between ",input$ssY," & ",input$ssX),
               main=paste0("MA plot of ",
                           switch(input$scatterInput,
                                  "MTC"="mean gene expression",
                                  "MDTC"="mean detected gene expression",
                                  "DR"="detection rate"),
                           " (",input$ssY," vs. ",input$ssX,")"),
               pch=20,col=alpha("black",0.3))
          lines(x=c(par("usr")[1],par("usr")[2]),y=c(par("usr")[3],par("usr")[3]),
                lwd=2,col=clustCols(res())[which(levels(clusts()) == input$ssX)],xpd=NA)
          lines(x=c(par("usr")[1],par("usr")[2]),y=c(par("usr")[4],par("usr")[4]),
                lwd=2,col=clustCols(res())[which(levels(clusts()) == input$ssY)],xpd=NA)
          if ("flip" %in% input$scatterLabelAngle) {
            temp_srtX <- 315
            temp_srtY <- 45
            temp_adjX <- c(-0.15,0.5)
            temp_adjY <- c(-0.15,0.5)
          } else {
            temp_srtX <- 45
            temp_srtY <- 315
            temp_adjX <- c(-0.15,0.5)
            temp_adjY <- c(-0.15,0.5)
          }
          if (input$scatterLine == "sub") {
            abline(h=0)
          } else if (input$scatterLine == "lbf") {
            abline(LBF())
          }
          if (input$diffLabelType == "diff") {
            temp_tY <- names(head(diffRanked(),input$diffCount))
            temp_tX <- names(tail(diffRanked(),input$diffCount))
            points(y~x,data=compDF()[temp_tY,],
                   pch=16,col=alpha(clustCols(res())[which(levels(clusts()) == input$ssY)],0.8))
            text(compDF()[temp_tY,"x"],compDF()[temp_tY,"y"],
                 labels=compDF()[temp_tY,"genes"],srt=temp_srtY,adj=temp_adjY,
                 col=clustCols(res())[which(levels(clusts()) == input$ssY)],cex=1.5,font=2)
            points(y~x,data=compDF()[temp_tX,],
                   pch=16,col=alpha(clustCols(res())[which(levels(clusts()) == input$ssX)],0.8))
            text(compDF()[temp_tX,"x"],compDF()[temp_tX,"y"],
                 labels=compDF()[temp_tX,"genes"],srt=temp_srtX,adj=temp_adjX,
                 col=clustCols(res())[which(levels(clusts()) == input$ssX)],cex=1.5,font=2)
          } else if (input$diffLabelType == "de") {
            degX <- rownames(switch(input$heatG,
                                    deTissue=d$deTissue[[res()]],
                                    deMarker=d$deMarker[[res()]],
                                    deNeighb=deNeighb[[res()]])[[input$ssX]])[1:input$diffCount]
            if (length(degX) > 0) {
              points(y~x,data=compDF()[degX,],
                     pch=16,col=alpha(clustCols(res())[which(levels(clusts()) == input$ssX)],0.8))
              text(compDF()[degX,"x"],compDF()[degX,"y"],
                   labels=compDF()[degX,"genes"],srt=temp_srtX,adj=temp_adjX,
                   col=clustCols(res())[which(levels(clusts()) == input$ssX)],cex=1.5,font=2)
            }
            degY <- rownames(switch(input$heatG,
                                    deTissue=d$deTissue[[res()]],
                                    deMarker=d$deMarker[[res()]],
                                    deNeighb=deNeighb[[res()]])[[input$ssY]])[1:input$diffCount]
            if (length(degY) > 0) {
              points(y~x,data=compDF()[degY,],
                     pch=16,col=alpha(clustCols(res())[which(levels(clusts()) == input$ssY)],0.8))
              text(compDF()[degY,"x"],compDF()[degY,"y"],
                   labels=compDF()[degY,"genes"],srt=temp_srtY,adj=temp_adjY,
                   col=clustCols(res())[which(levels(clusts()) == input$ssY)],cex=1.5,font=2)
            }
          } else if (input$diffLabelType == "search" & length(GOI()) > 0) {
            points(y~x,data=compDF()[GOI(),],pch=16,col=alpha("darkred",0.8))
            text(compDF()[GOI(),"x"],compDF()[GOI(),"y"],labels=compDF()[GOI(),"genes"],
                 srt=temp_srtX,adj=temp_adjX,col="darkred",cex=1.5,font=2)
          }
        }
      }
    }
    
    output$setScatter <- renderPlot(print(plot_setScatter()))
    
    output$setScatterSave <- downloadHandler(
      filename="setScatter.pdf",
      content=function(file) {
        pdf(file,width=7,height=7)
        print(plot_setScatter())
        dev.off()
      }
    )
    
    output$setComparisonSave <- downloadHandler(
      filename=function() { paste0(input$ssY,"vs",input$ssX,"_",
                                   input$scatterInput,"_",input$scatterLine,".txt") },
      content=function(file) {
        write.table(as.data.frame(diffRanked()),file,quote=F,sep="\t",row.names=T,col.names=F)
      }
    )
    
    
    # ^ Custom sets for DE ---------------------------------------------------------------
    selectedSets <- reactiveValues(a=NULL,b=NULL)
    
    plot_tsne_selDE <- function() {
      if (input$tsneSelDEcol == "") {
        id <- rep(1,nrow(md))
        idcol <- "grey20"
      } else if (is.factor(md[,input$tsneSelDEcol]) | is.character(md[,input$tsneSelDEcol])) {
        id <- as.factor(md[,input$tsneSelDEcol])
        if (length(levels(md[,input$tsneSelDEcol])) <= 8) {
          idcol <- RColorBrewer::brewer.pal(length(levels(md[,input$tsneSelDEcol])),
                              "Dark2")[1:length(levels(md[,input$tsneSelDEcol]))]
        } else {
          idcol <- rainbow2(length(levels(md[,input$tsneSelDEcol])))
        }
      } else {
        id <- cut(md[,input$tsneSelDEcol],100)
        idcol <- viridis(100,d=-1)
      }
      par(mar=c(3,3,3,1),mgp=2:0)
      plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
           xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
      points(dr_viz,pch=21,
             col=alpha(idcol,.8)[id],
             bg=alpha(idcol,0.4)[id])
      
      if (input$tsneSelDEcol == "") {
      } else if (is.factor(md[,input$tsneSelDEcol]) | is.character(md[,input$tsneSelDEcol])) {
        legend("topleft",bty="n",horiz=T,xpd=NA,inset=c(0,-.09),
               pch=21,col=idcol,pt.bg=alpha(idcol,0.5),
               title=input$tsneSelDEcol,legend=levels(md[,input$tsneSelDEcol]))
      } else {
        legend("topleft",bty="n",horiz=T,xpd=NA,inset=c(0,-.09),
               pch=21,col=viridis(3,d=-1),pt.bg=viridis(3,.5,d=-1),
               title=input$tsneSelDEcol,
               legend=c(round(min(md[,input$tsneSelDEcol]),2),
                        round((max(md[,input$tsneSelDEcol]) - 
                                 min(md[,input$tsneSelDEcol])) / 2,2),
                        round(max(md[,input$tsneSelDEcol]),2)))
      }
      
      
      points(dr_viz[selectedSets$a,],pch=19,col="#a50026")
      points(dr_viz[selectedSets$b,],pch=19,col="#313695")
      points(dr_viz[intersect(selectedSets$a,selectedSets$b),],pch=19,col="#ffffbf")
      points(dr_viz[intersect(selectedSets$a,selectedSets$b),],pch=4,col="red")
      
      legend("topright",horiz=T,bty="n",xpd=NA,inset=c(0,-.09),
             title="Selected Cells",legend=c("Set A","Set B","Both"),
             pch=c(19,19,4),col=c("#a50026","#313695","red"))
    }
    output$tsneSelDE <- renderPlot({ print(plot_tsne_selDE()) })
    
    
    currSel <- reactive(rownames(brushedPoints(as.data.frame(dr_viz),
                                               input$tsneBrush,xvar="tSNE_1",yvar="tSNE_2")))
    observeEvent(input$addCellsA,{ 
      selectedSets$a <- append(selectedSets$a,currSel()[!currSel() %in% selectedSets$a]) 
    })
    observeEvent(input$removeCellsA,{ 
      selectedSets$a <- selectedSets$a[!selectedSets$a %in% currSel()]
    })
    observeEvent(input$addCellsB,{ 
      selectedSets$b <- append(selectedSets$b,currSel()[!currSel() %in% selectedSets$b]) 
    })
    observeEvent(input$removeCellsB,{ 
      selectedSets$b <- selectedSets$b[!selectedSets$b %in% currSel()]
    })
    output$textSetA <- renderText(paste(length(selectedSets$a),"cells in Set A."))
    output$textSetB <- renderText(paste(length(selectedSets$b),"cells in Set B."))
    output$textOverlap <- renderText(
      paste(length(intersect(selectedSets$a,selectedSets$b)),"cells in both sets.",
            "Cells must be assigned to a single set prior to calculation.")
      )
    
    observeEvent(input$calcDE,{
      newRes <- paste0("Comp.",gsub("[^A-Za-z0-9]","",input$DEsetName))
      if (length(intersect(selectedSets$a,selectedSets$b)) > 0) {
        output$calcText <- renderText("Sets can't overlap (please assign red cells to only one set).")
      } else if (any(sapply(list(selectedSets$a,selectedSets$b),length) < 3)) {
        output$calcText <- renderText("Each set must contain at least 3 cells.")
      } else if (nchar(input$DEsetName) < 1) {
        output$calcText <- renderText("Please name this comparison (in text box above).")
      } else if (newRes %in% colnames(d$cl)) {
        output$calcText <- renderText("This comparison name has already been used.")
      } else {
        output$calcText <- renderText("")
        withProgress({
          temp_warn <- options("warn")
          options(warn=-1)
          
          temp <- rep("Unselected",nrow(d$cl))
          names(temp) <- rownames(d$cl)
          temp[selectedSets$a] <- "Set A"
          temp[selectedSets$b] <- "Set B"
          d$cl[[newRes]] <- factor(temp)
          
          # ^^ Gene stats per set --------------------------------------------------------
          incProgress(amount=1/6,detail="Gene detection rate per set")
          DR <- apply(nge,1,function(X) 
            tapply(X,d$cl[,newRes],function(Y) sum(Y>0)/length(Y)))
          
          incProgress(amount=1/6,detail="Mean detected gene expression per set")
          MDTC <- apply(nge,1,function(X) 
            tapply(X,d$cl[,newRes],function(Y) {
              temp <- meanLogX(Y[Y>0],ncell=ncol(nge),ex=exponent,pc=pseudocount)
              if (is.na(temp)) { temp <- 0 }
              return(temp)
            }))
          
          incProgress(amount=1/6,detail="Mean gene expression per set")
          MTC <- apply(nge,1,function(X) 
            tapply(X,d$cl[,newRes],function(Y) 
              meanLogX(Y,ncell=ncol(nge),ex=exponent,pc=pseudocount)))

          d$CGS[[newRes]] <- sapply(levels(d$cl[,newRes]),function(X) 
            data.frame(DR=DR[X,],MDTC=MDTC[X,],MTC=MTC[X,]),simplify=F)
          for (i in names(d$CGS[[newRes]])) {
            if (exists("symbolMap")) {
              d$CGS[[newRes]][[i]]$genes <- symbolMap[rownames(d$CGS[[newRes]][[i]])]
              d$CGS[[newRes]][[i]]$genes[is.na(d$CGS[[newRes]][[i]]$genes)] <- 
                rownames(d$CGS[[newRes]][[i]])[is.na(d$CGS[[newRes]][[i]]$genes)]
            } else {
              d$CGS[[newRes]][[i]]$genes <- rownames(d$CGS[[newRes]][[i]])
            }
            d$CGS[[newRes]][[i]]$MTCrank <- rank(d$CGS[[newRes]][[i]]$MTC,
                                                 ties.method="min")/nrow(d$CGS[[newRes]][[i]])
            if (i == "Unselected") { next }
            d$CGS[[newRes]][[i]]$cMu <- rownames(d$CGS[[newRes]][[i]]) %in% unlist(cellMarkersU)
            d$CGS[[newRes]][[i]]$cMs <- rownames(d$CGS[[newRes]][[i]]) %in% unlist(cellMarkersS)
            d$CGS[[newRes]][[i]]$overCut <- d$CGS[[newRes]][[i]]$MTC > mean(d$CGS[[newRes]][[i]]$MTC)
          }
          if (length(cellMarkers) < 1) {
            d$clusterID[[newRes]] <- sapply(d$CGS[[newRes]],function(Z) return(""))
          } else if (!any(unlist(cellMarkers) %in% rownames(nge))) {
            warning(paste("None of the provided cellMarkers are found in the data",
                          "(check your gene IDs against rownames in your data)."))
            d$clusterID[[newRes]] <- sapply(d$CGS[[newRes]],function(Z) return(""))
          } else {
            d$clusterID[[newRes]] <- c(names(cellMarkers)[sapply(d$CGS[[newRes]][1:2],function(Y) 
              which.max(sapply(cellMarkers,function(X) median(Y$MTC[rownames(Y) %in% X]))))],
              "Unselected")
            names(d$clusterID[[newRes]]) <- names(d$CGS[[newRes]])
          }
          
          # ^^ deTissue - DE per cluster vs all other data -------------------------------
          incProgress(amount=1/6,detail="DE vs tissue logGER calculations")
          deT_logGER <- sapply(levels(d$cl[,newRes])[1:2],function(i) 
            MTC[i,] - apply(nge[,d$cl[,newRes] != i],1,function(Y) 
              meanLogX(Y,ncell=ncol(nge),ex=exponent,pc=pseudocount)))
          deT_genesUsed <- apply(deT_logGER,2,function(X) which(X > logGERthresh))  
          if (any(sapply(deT_genesUsed,length) < 1)) {
            stop(paste0("logGERthresh should be set to less than ",
                        min(apply(deT_logGER,2,function(X) max(abs(X)))),
                        ", the largest magnitude logGER between cluster ",
                        names(which.min(apply(deT_logGER,2,function(X) max(abs(X))))),
                        " and the remaining data."))
          }
          incProgress(amount=1/6,detail="DE vs tissue Wilcoxon rank sum calculations")
          deT_pVal <- sapply(levels(d$cl[,newRes])[1:2],function(i)
            apply(nge[deT_genesUsed[[i]],],1,function(X) 
              wilcox.test(X[d$cl[,newRes] == i],X[d$cl[,newRes] != i])$p.value),simplify=F)
          d$deTissue[[newRes]] <- sapply(levels(d$cl[,newRes])[1:2],function(i) 
            data.frame(logGER=deT_logGER[deT_genesUsed[[i]],i],
                       pVal=deT_pVal[[i]])[order(deT_pVal[[i]]),],simplify=F)
          tempQval <- tapply(
            p.adjust(do.call(rbind,d$deTissue[[newRes]])$pVal,"fdr"),
            rep(names(sapply(d$deTissue[[newRes]],nrow)),sapply(d$deTissue[[newRes]],nrow)),
            c)
          for (i in names(d$deTissue[[newRes]])) { 
            d$deTissue[[newRes]][[i]] <- d$deTissue[[newRes]][[i]][tempQval[[i]] <= FDRthresh,]
            d$deTissue[[newRes]][[i]]$qVal <- tempQval[[i]][tempQval[[i]] <= FDRthresh] 
          }
          
          # ^^ deMarker - DE per cluster vs each other cluster ---------------------------
          incProgress(amount=1/6,detail="Calculating Set A vs Set B")
          
          deM_dDR <- DR["Set A",] - DR["Set B",]
          deM_logGER <- MTC["Set A",] - MTC["Set B",]
          deM_genesUsed <- switch(threshType,
                                  dDR=which(abs(deM_dDR) > dDRthresh),
                                  logGER=which(abs(deM_logGER) > logGERthresh))
          if (length(deM_genesUsed) < 1) {
            stop("Gene filtering threshold is set too high.")
          }
          
          deM_pVal <- apply(nge[deM_genesUsed,],1,function(X) 
            wilcox.test(X[d$cl[,newRes] == "Set A"],
                        X[d$cl[,newRes] == "Set B"])$p.value)
          
          temp_deVS <- data.frame(dDR=deM_dDR[deM_genesUsed],
                                  logGER=deM_logGER[deM_genesUsed],
                                  pVal=deM_pVal)[order(deM_pVal),]
          temp_deVS$qVal <- p.adjust(temp_deVS$pVal,"fdr")
          
          d$deMarker[[newRes]] <- list(
            "Set A"=temp_deVS[temp_deVS[,threshType] > 0 & temp_deVS$qVal <= FDRthresh,],
            "Set B"=temp_deVS[temp_deVS[,threshType] < 0 & temp_deVS$qVal <= FDRthresh,]
          )
          d$deMarker[[newRes]][["Set B"]]$dDR <- d$deMarker[[newRes]][["Set B"]]$dDR * -1
          d$deMarker[[newRes]][["Set B"]]$logGER <- d$deMarker[[newRes]][["Set B"]]$logGER * -1
          
          selectedSets$a <- selectedSets$b <- NULL
          options(warn=temp_warn$warn)
        },message="DE calculations:")      
        
        res(newRes) # Automatically update the view to show the calculated results.
      }
    })
    observeEvent(input$updateForViz, {
      withProgress({
        new_cl <- d$cl[input$res]
        new_CGS <- list()
        for (i in names(d$CGS[[input$res]])) {
          new_CGS[[input$res]][[i]] <- 
            d$CGS[[input$res]][[i]][colnames(d$CGS[[input$res]][[i]]) %in% c("DR","MDTC","MTC")]
        }
        new_deTissue <- d$deTissue[input$res]
        new_deMarker <- d$deMarker[input$res]
        incProgress(.5)
        save(new_cl,new_CGS,new_deTissue,new_deMarker,
             file=paste0(dataPath,dataTitle,"_selDE_",sub("Comp.","",input$res,fixed=T),".RData"))
      },message=paste0(
        "Saving ",dataTitle,"_selDE_",sub("Comp.","",input$res,fixed=T),".RData to ",dataPath))
    })
    observeEvent(input$updateForViz2, {
      withProgress({
        new_cl <- d$cl[input$res]
        new_CGS <- list()
        for (i in names(d$CGS[[input$res]])) {
          new_CGS[[input$res]][[i]] <- 
            d$CGS[[input$res]][[i]][colnames(d$CGS[[input$res]][[i]]) %in% c("DR","MDTC","MTC")]
        }
        new_deTissue <- d$deTissue[input$res]
        new_deMarker <- d$deMarker[input$res]
        incProgress(.5)
        save(new_cl,new_CGS,new_deTissue,new_deMarker,
             file=paste0(dataPath,dataTitle,"_selDE_",sub("Comp.","",input$res,fixed=T),".RData"))
      },message=paste0(
        "Saving ",dataTitle,"_selDE_",sub("Comp.","",input$res,fixed=T),".RData to ",dataPath))
    })
    

  }
  
  
  # Run the Shiny App! -------------------------------------------------------------------
  shinyApp(ui,server)
}
