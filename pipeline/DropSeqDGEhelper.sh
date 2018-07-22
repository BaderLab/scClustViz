#!/bin/bash
 
gunzip -k SORTED_out_gene_exon_tagged_dge.txt.gz

#rows
awk 'NR > 1 {print $1}' SORTED_out_gene_exon_tagged_dge.txt | awk 'END {print NR}'
awk 'NR > 1 {print $1}' SORTED_out_gene_exon_tagged_dge.txt > SORTED_out_gene_exon_tagged_dge_GENES.txt

#columns
head -n1 SORTED_out_gene_exon_tagged_dge.txt | cut -f1 --complement | awk '{print NF}'
head -n1 SORTED_out_gene_exon_tagged_dge.txt | cut -f1 --complement > SORTED_out_gene_exon_tagged_dge_BARCODES.txt

#matrix
awk 'NR > 1 {print $0}' SORTED_out_gene_exon_tagged_dge.txt | cut -f1 --complement | awk 'END {print NR, NF}' 
awk 'NR > 1 {print $0}' SORTED_out_gene_exon_tagged_dge.txt | cut -f1 --complement > SORTED_out_gene_exon_tagged_dge_MATRIX.txt

#cleanup
#rm SORTED_out_gene_exon_tagged_dge.txt




#### ADD NEW DENEIGHB METHOD TO PIPELINE ####
#### E11.5 and E17.5 swapped to new method ####
#### DropSeq loading code:

if (!file.exists(paste0(dataPath,"ebRaw.RData"))) {
  timePoints <- "e17"
  temp_path <- paste0(sys,"Dropbox/GDB/meCortex/SharedCortexData/")
  temp_preamble <- "e17/outE175_gene_exon_tagged_dge_"
  temp_cells <- scan(paste0(temp_path,temp_preamble,"BARCODES.txt"),character(),sep="\t")
  temp_genes <- scan(paste0(temp_path,temp_preamble,"GENES.txt"),character(),sep="\t")
  tempData <- Matrix(scan(paste0(temp_path,temp_preamble,"MATRIX.txt"),integer(),sep="\t"),
                     nrow=length(temp_genes),byrow=T)
  colnames(tempData) <- paste(timePoints,temp_cells,sep="_")
  rownames(tempData) <- temp_genes
  ebRaw <- data_load_processing(tempData)
  save(ebRaw,timePoints,file=paste0(dataPath,"ebRaw.RData"))
} else {
