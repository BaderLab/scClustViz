tp <- "e17"

load(paste0("meCortex/",tp,"/",tp,"_Cortical_Only_deVS.RData"))
for (res in names(deVS)) {
  for (x in names(deVS[[res]])) {
    for (y in names(deVS[[res]][[x]])) {
      print(colnames(deVS[[res]][[x]][[y]]))
      colnames(deVS[[res]][[x]][[y]])[colnames(deVS[[res]][[x]][[y]]) == "logFC"] <- "logGER"
      print(colnames(deVS[[res]][[x]][[y]]))
    }
  }
}
save(deVS,file=paste0("meCortex/",tp,"/",tp,"_Cortical_Only_deVS.RData"))

temp <- load(paste0("meCortex/",tp,"/",tp,"_Cortical_Only_forViz.RData"))
for (res in names(deTissue)) {
  for (x in names(deTissue[[res]])) {
    print(colnames(deTissue[[res]][[x]]))
    colnames(deTissue[[res]][[x]])[colnames(deTissue[[res]][[x]]) == "logFC"] <- "logGER"
    print(colnames(deTissue[[res]][[x]]))
  }
}
for (res in names(deMarker)) {
  for (x in names(deMarker[[res]])) {
    print(colnames(deMarker[[res]][[x]]))
    colnames(deMarker[[res]][[x]]) <- sub("logFC","logGER",colnames(deMarker[[res]][[x]]))
    print(colnames(deMarker[[res]][[x]]))
  }
}
for (res in names(deNeighb)) {
  for (x in names(deNeighb[[res]])) {
    print(colnames(deNeighb[[res]][[x]]))
    colnames(deNeighb[[res]][[x]]) <- sub("logFC","logGER",colnames(deNeighb[[res]][[x]]))
    print(colnames(deNeighb[[res]][[x]]))
  }
}
save(list=temp,file=paste0("meCortex/",tp,"/",tp,"_Cortical_Only_forViz.RData"))


rm(list=ls())
