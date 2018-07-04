# Note: Change egDB code from main RunVizScript to mouse only version from e1*/app.R

library(rsconnect)
options(repos = BiocInstaller::biocinstallRepos())
setwd("meCortex/e11")
deployApp(appName="e11cortex",
          appFiles=c("app.R","e11_Cortical_Only_forViz.RData","e11_Cortical_Only_savedRes.RData","intro.md"))
