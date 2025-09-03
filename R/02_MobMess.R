
# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")
library(data.table)
library("Biostrings")
library(RColorBrewer)
# display.brewer.all()
palette<-brewer.pal(n=12,name='Paired')

# set working directory
setwd("~/PATH/")

# functions
get_ID <- function(x) {
  x <- toupper(head(unlist(x),n=1))
  return(x)
}

#############################  MobMess   #######################################
# Load clusters
files_cluster <- list.files(path = "MobMess/output/clusters/", 
                            full.names = TRUE)
names(files_cluster) <- as.integer(gsub("ID", "", str_remove(basename(files_cluster), 
                                    "\\-mobmess_clusters.txt")))
df_cluster <- map_dfr(files_cluster, 
                      fread, .id = 'child')

df_count <- df_cluster %>% 
  group_by(child) %>% 
  mutate(child = as.integer(child),
         contigs = n()) %>% 
  select(child,cluster, contigs) %>% 
  unique() %>% 
  mutate(plasmids = n()) %>% 
  select(child, contigs,plasmids) %>% 
  unique() %>% 
  arrange(child)
names(df_count) <- c("child", "#plasmids", "#contigs")

file_count <- "count_contigs_plasmids.csv"
fwrite(df_count, 
       file = file_count,
       sep = ",", 
       row.names=FALSE, 
       col.names=FALSE)

############################  add geNomad  #####################################
## Make summary of mobmess and genomad
df_summary <- df_cluster %>% 
  select(child,cluster, systems,cluster_type, 
         representative_contig) %>%
  arrange(child,cluster) %>%
  mutate(seq_name = gsub(".", "_", representative_contig, fixed = T),
         child = paste("ID", child, sep = "")) %>%
  select(-representative_contig) %>%
  unique()

#### add genomad to the summary
file_genomad <- "geNomad/plasmid_summary.tsv"
df_genomad <- fread(file_genomad, 
                    sep = "\t",
                    fill = TRUE)
df_genomad <- df_genomad %>%
  mutate(child = sapply(strsplit(seq_name, "_"), get_ID))

df_mobmess_genomad <- left_join(df_summary, df_genomad)

## Save summary mobmess and genomad
file_summary <- "summary_mobmess_genomad.csv"
fwrite(df_mobmess_genomad, 
       file = file_summary,
       sep = ",", 
       row.names=FALSE, 
       col.names=TRUE)

########################  fasta of plasmids   ##################################
### Load contigs
file_contigs <- "scaffolds/scaffolds.fasta"
df_contigs <- readDNAStringSet(file_contigs)

### select reference sequences
df_contigs_filter <- extractList(df_contigs, df_mobmess_genomad$seq_name)
df_contigs_filter <- unlist(df_contigs_filter)

### rename to "child_plasmid"
newnames <- df_mobmess_genomad %>%
  mutate(
    newname = paste(child, cluster, sep = "_cluster")
  )

names(df_contigs_filter) <- newnames$newname

## split per child and save fasta
for (i in 1:12) {
  tmp_child <- paste("ID", i, sep = "")
  
  # split fasta
  row_filter <- extractList(df_contigs_filter, 
                            names(df_contigs_filter)[grepl(paste0(tmp_child, "_"), 
                            names(df_contigs_filter))])
  row_filter <- unlist(row_filter)
  names(row_filter) <- gsub("_", ".", names(row_filter))
  row_filter <- row_filter %>% as.character()
  
  df_contigs_child <- DNAStringSet(row_filter)
  
  file_contigs_child <- paste0("scaffolds/",
                               tmp_child, 
                               "_plasmids.fasta")
  writeXStringSet(df_contigs_child ,file_contigs_child)
}

## Make circular file for second round of MobMess
file_circular <- "MobMess/circular_scaffolds_geNomad.txt"
df_circular <- fread(file_circular,
                 sep = "\t",
                 col.names = c("seq_name", "circular"))

df_circular_new <- left_join(df_circular, df_mobmess_genomad) %>%
  filter(!is.na(child)) %>%
  mutate(ID = paste(toupper(child), ".cluster", cluster ,sep = "")) %>%
  arrange(as.integer(gsub("ID","",child)), cluster) %>%
  select(ID, circular)

file_mobmess_circ <- "MobMess/circular_all_plasmids.txt"
fwrite(df_circular_new, 
       file = file_mobmess_circ,
       sep = "\t", 
       row.names=FALSE, 
       col.names=FALSE)

############ Count plots ###################
df_count_plot <- as.matrix(df_count)
rownames(df_count_plot) <-  as.character(seq(1,12,1))
df_count_plot <- df_count_plot[,2:3]

grDevices::windows(15,15)
par(mfrow=c(1,1))
par(mar=c(4,6,2,1))

barplot(t(df_count_plot),
        col = c(palette[2], palette[1]),
        # legend = c("contigs","bins"),
        beside = TRUE,
        las = 3,
        cex.axis=1.5, cex.names=1.5,
        cex.lab=1.5,
        xlab = "Infant ID",
        ylab = "# Contigs/Plasmids",
        ylim = c(0,250))
legend("topleft", 
       legend = c("Contigs","Plasmids"), 
       fill = c(palette[2], palette[1]),
       cex = 1.5,
       inset = .03)
