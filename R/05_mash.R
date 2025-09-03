## filter mash hits

# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")
library(data.table)

# Set working directory
setwd("~/PATH/")

## functions
get_child <- function(x) {
  return(head(unlist(x), n=1))
}
get_cluster <- function(x) {
  return(as.integer(gsub("cluster","",
                         tail(unlist(x), n=1)))) 
}
get_species <- function(x) {
  x <- unlist(x)
  x <- paste(head(x, n=2), collapse = " ")
  return(x)
}
get_genus <- function(x) {
  x <- unlist(x)
  x <- head(x, n=1)
  return(x)
}

###############################   infants  #####################################
children <-  paste("ID",seq(1:12), sep = "")

## Load mash files and concat in one dataframe
df_mash <- data.frame()
for (j in 1:length(children)) {
  tmp_child <- children[j]
  
  mash_files <- list.files(path = paste0("mash/output/",
                                         tmp_child,
                                         "/"),
                           full.names = TRUE)
  names(mash_files) <- str_remove(basename(mash_files), 
                                      "\\.tab")
  tmp_row <- map_dfr(mash_files, 
                         fread, .id = 'ID')
  colnames(tmp_row) <- c("ID",
                         "identity", 
                         "shared-hashes", 
                         "median-multiplicity", 
                         "p-value", 
                         "query-ID", 
                         "query-comment")
  df_mash <- rbind(df_mash,
                   tmp_row)
}

df_mash <- df_mash %>%
  mutate(child = sapply(strsplit(ID, ".", fixed = T), get_child),
         cluster = sapply(strsplit(ID, ".", fixed = T), get_cluster)) %>%
  select(-ID)

### Add Redondo-Salvo 2020 PTUs
file_ptu <- "PTUs/Redondo-Salvo_2020_Suppl_2.csv"
df_PTU <- fread(file_ptu,
                sep = ";",
                fill = TRUE,
                header = TRUE)
df_PTU <- df_PTU %>%
  mutate("query-ID" = AccessionVersion,
         PTU = ifelse(PTU == "-", NA, PTU),
         Hrange = ifelse(Hrange == "-", NA, Hrange)) %>%
  select(`query-ID`, PTU, Hrange, Filtered)

df_mash_PTU <- left_join(df_mash, df_PTU)

## summarize species
df_mash_species <- df_mash_PTU %>% 
  mutate(species = sapply(strsplit(`query-comment`," "), get_species)) %>%
  filter(!species %in% c("Uncultured bacterium",
                        "UNVERIFIED: Uncultured",
                        "Endophytic bacterium")) %>%
  unite(ID, c(child,cluster,species),sep="_",remove = F) 
count_species <- as.data.frame(table(df_mash_species$ID))
colnames(count_species) <- c("ID", "count")
df_mash_species_count <- left_join(df_mash_species, count_species) %>%
  select(child, cluster, species, count, PTU, Hrange, Filtered) %>%
  unique()

df_mash_species_summary <- df_mash_species_count %>%
  group_by(child, cluster) %>%
  mutate("PLSDB species" = paste(species, collapse = "; "),
         "PLSDB species count" = paste(count, collapse = "; "),
         "PTUs" = paste(unique(na.omit(PTU)), collapse = "; "),
         "Hranges" = paste(unique(na.omit(Hrange)), collapse = "; ")) %>%
  select(child, cluster, `PLSDB species`, `PLSDB species count`, PTUs, Hranges) %>%
  unique()

## summarize genus
df_mash_genus <- df_mash %>%
  group_by(child, cluster) %>% 
  mutate(genus = sapply(strsplit(`query-comment`," "), get_genus),
         "PLSDB NCBI" = paste(`query-ID`, collapse = "; ")) %>%
  ungroup() %>%
  filter(!genus %in% c("Uncultured",
                       "UNVERIFIED:",
                       "Endophytic")) %>%
  unite(ID, c(child,cluster,genus),sep="_",remove = F) 
count_genus <- as.data.frame(table(df_mash_genus$ID))
colnames(count_genus) <- c("ID", "count")
df_mash_genus_count <- left_join(df_mash_genus, count_genus) %>%
  select(child, cluster, genus, count, `PLSDB NCBI`) %>%
  unique()

df_mash_genus_summary <- df_mash_genus_count %>%
  group_by(child, cluster) %>%
  mutate("PLSDB genus" = paste(genus, collapse = "; "),
         "PLSDB genus count" = paste(count, collapse = "; ")) %>%
  select(child, cluster, `PLSDB genus`, `PLSDB genus count`, `PLSDB NCBI`) %>%
  unique()

df_safe <- left_join(df_mash_genus_summary,df_mash_species_summary)

## Save dataframe with PLSDB info
outfile <- "plasmids_PLSDB_taxa_clusters.csv"
fwrite(df_safe, 
       file = outfile,
       sep = ",", 
       row.names=FALSE, 
       col.names=TRUE)

