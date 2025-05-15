
# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")
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

df_cluster <- data.frame()
df_count <- data.frame()

for (i in 1:length(files_cluster)) {
  tmp_child <- as.integer(gsub("ID", "",tail(unlist(
                                strsplit(head(unlist(
                                  strsplit(files_cluster[i], "-")),
                                    n=1), "/")),n=1)))
  tmp_file <- files_cluster[i]
  
  tmp_cluster <- read_tsv(tmp_file, 
                          show_col_types = FALSE)
  tmp_cluster <- tmp_cluster %>%
    mutate(child = tmp_child)
  
  df_cluster <- rbind(df_cluster,
                      tmp_cluster)
  
  ##count cluster (plasmids) and system per infant
  tmp_plasmids <- max(tmp_cluster$cluster)
  tmp_contigs <- nrow(tmp_cluster)
  tmp_row <- c(tmp_child, tmp_plasmids,tmp_contigs)
  df_count <- rbind(df_count,
                    tmp_row)
}

# save count for later plotting
names(df_count) <- c("child", "#plasmids", "#contigs")
df_count <- df_count %>%
  select(child,
         `#contigs`,
         `#plasmids`) %>%
  arrange(child)
file_count <- "count_contigs_plasmids.csv"
write_csv(df_count, file = file_count)

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
df_genomad <- read_tsv(file_genomad)
df_genomad <- df_genomad %>%
  mutate(child = sapply(strsplit(seq_name, "_"), get_ID))

df_mobmess_genomad <- left_join(df_summary, df_genomad)

## Save summary mobmess and genomad
file_summary <- "summary_mobmess_genomad.csv"
write_csv(df_mobmess_genomad, file = file_summary)

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
df_circular <- read_tsv(file_circular,
                        col_names = F)
colnames(df_circular) <- c("seq_name", "circular")

df_circular_new <- left_join(df_circular, df_mobmess_genomad) %>%
  filter(!is.na(child)) %>%
  mutate(ID = paste(toupper(child), ".cluster", cluster ,sep = "")) %>%
  arrange(as.integer(gsub("ID","",child)), cluster) %>%
  select(ID, circular)

file_mobmess_circ <- "MobMess/circular_all_plasmids.txt"
write.table(df_circular_new, file = file_mobmess_circ,
            sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

############ Count plots ###################

rownames(df_count) <- df_count$child
df_count <- df_count %>% select(-child)
df_count <- as.matrix(df_count, ncol= 2)

df_count <- df_count[ order(as.integer(rownames(df_count))), ]

grDevices::windows(15,15)
par(mfrow=c(1,1))
par(mar=c(4,6,2,1))

barplot(t(df_count),
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
