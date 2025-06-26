## filter mash hits

# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")

# Set working directory
setwd("~/PATH/")

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
df_mash <- data.frame(child = as.character(),
                      cluster=as.integer(),
                      identity=as.double(),
                      sharedhashes=as.character(),
                      medianmultiplicity=as.integer(),
                      pvalue=as.double(),
                      queryID=as.character(),
                      querycomment=as.character())
for (j in 1:length(children)) {
  tmp_child <- children[j]
  
  mash_files <- list.files(path = paste0("PATH/",
                                         tmp_child,
                                         "/"),
                           full.names = TRUE)
  for (i in 1:length(mash_files)) {
    file_name <- tail(strsplit(mash_files[i], "/")[[1]], n=1)
    file_name <- paste(head(strsplit(file_name, "[.]")[[1]], n=2), 
                       collapse = "_")
    sample_name <- tail(strsplit(file_name, "_")[[1]], n=1)
    assign(file_name, read_tsv(mash_files[i], 
                               col_names = c("identity", 
                                             "shared-hashes", 
                                             "median-multiplicity", 
                                             "p-value", 
                                             "query-ID", 
                                             "query-comment")))
    tmp_row <- get(file_name)
    tmp_row <- tmp_row %>%
      mutate(
        child = tmp_child,
        cluster = as.integer(gsub("cluster","",sample_name))
      )
    df_mash <- rbind(df_mash, tmp_row)
  }
}

## summarize species
df_mash_species <- df_mash %>% 
  mutate(species = sapply(strsplit(`query-comment`," "), get_species)) %>%
  filter(!species %in% c("Uncultured bacterium",
                        "UNVERIFIED: Uncultured",
                        "Endophytic bacterium")) %>%
  unite(ID, c(child,cluster,species),sep="_",remove = F) 
count_species <- as.data.frame(table(df_mash_species$ID))
colnames(count_species) <- c("ID", "count")
df_mash_species_count <- left_join(df_mash_species, count_species) %>%
  select(child, cluster, species, count) %>%
  unique()

df_mash_species_summary <- df_mash_species_count %>%
  group_by(child, cluster) %>%
  mutate("PLSDB species" = paste(species, collapse = "; "),
         "PLSDB species count" = paste(count, collapse = "; ")) %>%
  select(child, cluster, `PLSDB species`, `PLSDB species count`) %>%
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
write.csv(df_safe, file = outfile, row.names = FALSE, quote = FALSE)

