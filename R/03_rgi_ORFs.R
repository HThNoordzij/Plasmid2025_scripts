## clusters, make prodigal .fsa file from rgi output per cluster

# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")
library(data.table)
library("Biostrings")

## functions
get_ORF <- function(x){
  x <- head(unlist(x), n=1)
  return(x)
}

# Set working directory
setwd("~/PATH/")

###############################    ORFs    #####################################
# Load clusters
file_clusters <- "summary_mobmess_genomad.csv"
df_clusters <- fread(file_clusters, 
                     sep = ",",
                     fill = TRUE)
children <- df_clusters$child %>% unique()

for(i in 1:length(children)) {
  tmp_child <- children[i]
  
  df_clusters_reference <- df_clusters %>%
    filter(child == tmp_child) %>%
    select(cluster) %>%
    unique()
  
  # Load contigs
  file_orfs <- paste0("PATH/",
                      tmp_child,
                      "_plasmids.fasta.temp.contig.fsa")
  df_orfs <- readAAStringSet(file_orfs)
  
  # rename sequences to only include ORF
  fa_given_names <- names(df_orfs)
  fa_new_names <- sapply(strsplit(fa_given_names, " "), get_ORF)
  names(df_orfs) <- fa_new_names
  
  clusters_no <- df_clusters_reference$cluster %>% unique()
  
  for (j in 1:length(clusters_no)) {
    tmp_cluster <- clusters_no[j]
    df_orfs_cluster <- df_orfs[grepl(paste0("cluster",tmp_cluster, "_"), 
                                     fa_new_names)]
    
    ## make output dir
    output_dir <- file.path(paste0("PATH/prodigal/",
                                   toupper(tmp_child)))
    
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    } 
    
    # Safe ORFs of clusters as fasta
    file_ORFs_clusters_fasta <-  paste0("PATH/prodigal/",
                                    toupper(tmp_child),
                                    "/cluster", 
                                    tmp_cluster, 
                                    "_ORFs.fasta")
    writeXStringSet(df_orfs_cluster ,file_ORFs_clusters_fasta)
  }
}

