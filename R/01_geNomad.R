# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")
library("Biostrings")

# set working directory
setwd("~/PATH/")

# functions
get_ID <- function(x) {
  x <- toupper(head(unlist(x),n=1))
  return(x)
}


###############################   geNomad  #####################################
### Load geNomad file of metaplasmidSPAdes contigs
file_genomad_spades <- "PATH/plasmid_summary.tsv"
df_genomad_spades <- read_tsv(file_genomad_spades)
df_genomad <- df_genomad_spades %>%
  mutate(child = sapply(strsplit(seq_name, "_"), get_ID))


###############################    fasta   #####################################
### Load contigs of metaplasmidSPAdes
file_contigs_spades <- "PATH/scaffolds.fasta"
df_contigs_spades <- readDNAStringSet(file_contigs_spades)

### select reference sequences
df_contigs_spades_filter <- extractList(df_contigs_spades, df_genomad$seq_name)
df_contigs_spades_filter <- unlist(df_contigs_spades_filter)

# Save potential plasmids as fasta
file_contigs_bins <- "scaffolds/scaffolds_geNomad.fasta"
writeXStringSet(df_contigs_spades_filter ,file_contigs_bins)

# Save circular file for MobMess
file_mobmess_circ <- "MobMess/circular_scaffolds_geNomad.txt"
df_mobless_circ <- df_genomad %>%
  mutate(circ = case_when(topology == "DTR" ~ 1,
                          topology == "ITR" ~ 1,
                          .default = 0)) %>%
  select(seq_name, circ)

write.table(df_mobless_circ, file = file_mobmess_circ,
            sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### split per child
for (i in 1:12) {
  tmp_child <- paste("id", i, sep = "")
  
  # split fasta
  row_filter <- extractList(df_contigs_spades_filter, 
                            names(df_contigs_spades_filter)[grepl(paste0(tmp_child, "_"), 
                                                                  names(df_contigs_spades_filter))])
  row_filter <- unlist(row_filter)
  names(row_filter) <- gsub("_", ".", names(row_filter))
  row_filter <- row_filter %>% as.character()
  
  df_contigs_child <- DNAStringSet(row_filter)
  
  file_contigs_child <- paste0("scaffolds/scaffolds_geNomad_",
                              tmp_child, 
                              ".fasta")
  writeXStringSet(df_contigs_child ,file_contigs_child)
  
  # split circ file for mobmess
  df_mobless_circ_child <- df_mobless_circ %>%
    filter(grepl(paste0(tmp_child, "_"), seq_name)) %>%
    mutate(seq_name = gsub("_", ".", seq_name))
  
  write.table(df_mobless_circ_child, file = paste0("MobMess/circular_scaffolds_geNomad_",
                                                   tmp_child,
                                                   ".txt"),
              sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
