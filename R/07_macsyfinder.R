# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

# Load libraries
library("tidyverse")
library("Biostrings")

# Set working directory
setwd("~/PATH/")

## functions
get_cluster <- function(x) {
  x <- as.integer(unlist(strsplit(unlist(x)[1], "cluster"))[2])
}
# make type name shorter
adjust_type <- function(x) {
  tmp <- tail(unlist(x), n=1)
  return(tmp)
}

###############################   infants  #####################################
# Load annotation summary file
annotation_file <- "plasmids_annotation_clusters.csv"
df_annotation <- read_csv(annotation_file)

## clusters with oriT
df_annotation_ori <- df_annotation %>% filter(DB == "OriT") %>% select(cluster) %>% unique()

# empty dataframe
df_annotation_conj_child <- data.frame()
children <- df_annotation$child %>% unique()

for (i in 1:length(children)) {
  tmp_child <- children[i]
  
  ## Load macsyfinder CONJ summary files in one dataframe
  df_macsyfinder_CONJ <- data.frame()
  macsyfinder_CONJ_files <- list.files(
    path = paste0("macsyfinder/",
                  tmp_child,
                  "/best_solution_summary/"), 
    full.names = TRUE)
  for (i in 1:length(macsyfinder_CONJ_files)) {
    file_name <- tail(strsplit(macsyfinder_CONJ_files[i], "/")[[1]], n=1)
    file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
    assign(file_name, read.table(macsyfinder_CONJ_files[i], 
                                 comment.char = "#",
                                 header = T))
    
    tmp_df <- get(file_name)
    tmp_df <- tmp_df %>%
      mutate(child = tmp_child)
    
    df_macsyfinder_CONJ <- rbind(df_macsyfinder_CONJ, tmp_df)
  }
  ## add cluster column
  df_macsyfinder_CONJ <- df_macsyfinder_CONJ %>%
    mutate(cluster = sapply(strsplit(replicon, "_"), get_cluster)) %>%
    select(-replicon)
  
  ## clusters per type
  # pivot table
  df_plasmid_type <- df_macsyfinder_CONJ %>% 
    select(-child) %>%
    pivot_longer(!cluster, names_to = "type", values_to = "count") %>%
    mutate(type = sapply(strsplit(type, "[.]"), adjust_type)) 
  ## filter out zero's from df_plasmid_type
  df_plasmid_type <- df_plasmid_type %>% 
    filter(count > 0) %>%
    mutate(conjugation = type,
           child = tmp_child) %>%
    select(child,cluster, conjugation)
  
  ## add mobility type to annotation dataframe
  df_annotation_conj_child <- rbind(df_annotation_conj_child, df_plasmid_type)
}
df_annotation_conj <- left_join(df_annotation, df_annotation_conj_child)

## add oriT as conjugation when no other conjugation type
df_annotation_conj <- df_annotation_conj %>%
  mutate(conjugation = case_when(
    is.na(conjugation) & cluster %in% df_annotation_ori$cluster ~ "OriT",
    .default = conjugation
  ))

## Load taxonomy data and add plasmid mobility there
file_taxa <- "plasmids_taxonomy_clusters.csv"
df_taxonomy <- read_csv(file_taxa)

df_taxonomy_mobility <- left_join(df_taxonomy, df_annotation_conj) 
df_taxonomy_mobility <- df_taxonomy_mobility %>%
  mutate(
    "mating-pair formation class" = case_when(
      conjugation == "MOB" ~ NA,
      conjugation == "OriT" ~ NA,
      .default = conjugation
      ),
    mobility = case_when(
      grepl("MOB", conjugation) ~ "MOB",
      grepl("OriT", conjugation) ~ "MOB",
      grepl("T4SS", conjugation) ~ "CONJ",
      grepl("dCONJ", conjugation) ~ "MOB",
      .default = "Mobless"
    )) %>%
  select(
    child,
    cluster,
    length,
    contig,
    `putative host`,
    `putative broad host range (pBHR) order`,
    mobility,
    `mating-pair formation class`,
    `putative host assignment origin`,
    `PLSDB genera`,
    `PLSDB genera count`,
    `GTDB plasmid genera`,
    `GTDB genome genera`,
    `UniRef90 genera`,
    `UniRef90 protein count per genera`,
    `PLSDB species`,
    `PLSDB species count`,
    `PLSDB NCBI`,
    `GTDB plasmid species`,
    `GTDB genome species`
    ) %>%
  unique()

######## Add MobMess plasmid systems
file_mobmess_system <- "MobMess/output/all_plasmids/all-mobmess_clusters.txt"
df_mobmess <- read_tsv(file_mobmess_system) %>%
  mutate(cluster_new = cluster) %>%
  select(-cluster)

df_taxonomy_mobility_mobmess <- left_join(df_taxonomy_mobility, df_mobmess) %>%
  select(
    child,
    cluster,
    cluster_new,
    systems,
    length,
    representative_contig,
    cluster_type,
    contig,
    `putative host`,
    `putative broad host range (pBHR) order`,
    mobility,
    `mating-pair formation class`,
    `putative host assignment origin`,
    `PLSDB genera`,
    `PLSDB genera count`,
    `GTDB plasmid genera`,
    `GTDB genome genera`,
    `UniRef90 genera`,
    `UniRef90 protein count per genera`,
    `PLSDB species`,
    `PLSDB species count`,
    `PLSDB NCBI`,
    `GTDB plasmid species`,
    `GTDB genome species`)

## Rename cluster, 
##  so different infants with the same plasmid have the same cluster number
df_rename <- df_taxonomy_mobility_mobmess %>%
  mutate(contig_old = contig,
         contig_new = paste(child, ".cluster", cluster_new, sep = "")) %>%
  select(child, cluster, cluster_new,
         contig_old, contig_new)
## keep renaming dataframe
file_rename <- "rename_plasmids.csv"
write_csv(df_rename, file = file_rename)

df_taxonomy_mobility_summary <- df_taxonomy_mobility_mobmess %>%
  mutate("cluster type" = cluster_type,
         cluster = cluster_new,
         "representative contig child" = sapply(contig, function(x) {
           return(df_rename$contig_new[which(
             df_rename$contig_old == x)])}  ),
         "representative contig system" = sapply(representative_contig, function(x) {
           return(df_rename$contig_new[which(
             df_rename$contig_old == x)])}  ))  %>%
  select(
    child,
    cluster,
    systems,
    length,
    `representative contig system`,
    `cluster type`,
    `representative contig child`,
    `putative host`,
    `putative broad host range (pBHR) order`,
    mobility,
    `mating-pair formation class`,
    `putative host assignment origin`,
    `PLSDB genera`,
    `PLSDB genera count`,
    `GTDB plasmid genera`,
    `GTDB genome genera`,
    `UniRef90 genera`,
    `UniRef90 protein count per genera`,
    `PLSDB species`,
    `PLSDB species count`,
    `PLSDB NCBI`,
    `GTDB plasmid species`,
    `GTDB genome species`) %>%
  arrange(cluster, child)

## Save dataframe with conjugation info
outfile <- "plasmids_tax_mobility_clusters.csv"
write_csv(df_taxonomy_mobility_summary, file = outfile)


#####################  rename plasmid in fasta   ###############################
### Load contigs
file_contigs <- "scaffolds/all_plasmids.fasta"
df_contigs <- readDNAStringSet(file_contigs)

names(df_contigs) <- df_rename$contig_new[which(
  df_rename$contig_old == names(df_contigs))]

file_plasmids <- paste0("scaffolds/plasmids.fasta")
writeXStringSet(df_contigs ,file_plasmids)
