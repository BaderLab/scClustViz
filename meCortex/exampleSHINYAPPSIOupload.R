# Note: Change egDB code from main RunVizScript to mouse only version from e1*/app.R

tp <- "e13"

library(rsconnect)
options(repos = BiocInstaller::biocinstallRepos())
setwd(paste0("meCortex/",tp))

deployApp(appName=paste0(tp,"cortex"),
          appFiles=c("app.R","intro.md",
                     paste0(tp,"_Cortical_Only_forViz.RData"),
                     paste0(tp,"_Cortical_Only_savedRes.RData")))
