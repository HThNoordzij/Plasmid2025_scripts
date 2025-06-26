# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")

# set working directory
setwd("~/PATH/")

##functions
get_child <- function(x) {
  x <- head(unlist(x),n=1)
  return(toupper(x))
}
get_contig <- function(x){
  x <- head(unlist(x),n=1)
  return(x)
}
get_cluster <- function(x) {
  x <- head(unlist(x), n=1)
  x <- tail(unlist(strsplit(x, ".", fixed = T)),n=1)
  x <- as.integer(gsub("cluster", "", x))
  return(x)
}
get_taxa_dorit <- function(x) {
  x <- unlist(x)
  x <- paste(x[1:2], collapse = " ")
  return(x)
}
get_taxa_orit <- function(x) {
  x <- unlist(x)
  x <- paste(x[-1:-2], collapse = " ")
  return(x)
}
get_function_uniprot <- function(x) {
  x <- head(unlist(x),n=1)
  x <- unlist(strsplit(x, " "))[-1]
  x <- paste(x, collapse = " ")
  return(x)
}
get_taxa_uniprot <- function(x) {
  x <- tail(unlist(x),n=1)
  x <- head(unlist(strsplit(x, " TaxID")), n=1)
  return(x)
}
get_taxid <- function(x) {
  x <- tail(unlist(x),n=1)
  x <- head(unlist(strsplit(x, " ")),n=1)
  return(as.integer(x))
}

###############################################################################
##############################   infants  #####################################
###############################################################################

## Load prodigal (intermediate files from running RGI)
# create empty dataframe
df_prodigal <- data.frame(stringsAsFactors = FALSE)

## Load mobileOG files into one dataframe
prodigal_files <- list.files(path = "PATH/", 
                             full.names = TRUE)
for (i in 1:length(prodigal_files)) {
  file_name <- tail(strsplit(prodigal_files[i], "/")[[1]], n=1)
  file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
  assign(file_name, read.table(prodigal_files[i], sep="#", 
                               comment.char="", 
                               fill = TRUE))
  row <- get(file_name)
  row <- row %>% filter(!is.na(V2)) 
  colnames(row) <- c("ORF", "start", "end", "Orientation", "info")
  
  tmp_child <- tail(head(unlist(strsplit(file_name, "_")), n=3), n=1)
  row <- row %>%
    mutate(child = tmp_child,
           ORF = gsub(">", "", ORF))
  
  df_prodigal <- rbind(df_prodigal, row)
}
df_prodigal <- df_prodigal %>%
  mutate(ORF = gsub(" ", "", gsub(">", "", ORF)),
         child = sapply(strsplit(ORF, ".", fixed = T), get_child),
         contig = sapply(strsplit(ORF, "_"), get_contig),
         cluster = sapply(strsplit(ORF, "_"), get_cluster),
         Orientation = ifelse(Orientation == 1, "+", "-"),
         partial = case_when(
           grepl("partial=00", info) ~ FALSE,
           .default = TRUE
         )
  )
df_prodigal_select <- df_prodigal %>%
  select(child, cluster, contig, ORF, start, end, Orientation, info, partial) %>%
  unique()

df_prodigal_total <-  df_prodigal_select %>%
  mutate(Name = NA,
         DB = "Prodigal",
         ARO = NA,
         "Drug Class" = NA,
         "Resistance Mechanism" = NA,
         UniRef_function = NA,
         DB_taxa = NA,
         ID = contig) %>%
  select(child,
         cluster,
         contig,
         ORF,
         Name,
         DB,
         ARO,
         `Drug Class`,
         `Resistance Mechanism`,
         UniRef_function,
         DB_taxa,
         ID,
         start,
         end,
         Orientation,
         info,
         partial
  )


## load rgi
# create empty dataframe
df_rgi <- data.frame(stringsAsFactors = FALSE)

## Load rgi files into one dataframe
rgi_files <- list.files(path = "PATH/", 
                        full.names = TRUE)
for (i in 1:length(rgi_files)) {
  file_name <- tail(strsplit(rgi_files[i], "/")[[1]], n=1)
  file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
  assign(file_name, read_tsv(rgi_files[i], col_names = TRUE))
  row <- get(file_name)
  child <- tail(head(unlist(strsplit(file_name, "_")), n=3),n=1)
  row$child <- child
  df_rgi <- rbind(df_rgi, row)
}
df_rgi <- df_rgi %>%
  mutate(child = sapply(strsplit(Contig, ".", fixed = T), get_child),
         cluster = sapply(strsplit(Contig, "_"), get_cluster))

# filter
df_rgi_filter <- df_rgi %>% 
  filter(Best_Identities >= 80 &
           `Percentage Length of Reference Sequence` >= 80) 

## select right columns rgi
df_rgi_total <- df_rgi_filter%>%
  mutate(ORF = Contig,
         contig = sapply(strsplit(Contig, "_"), get_contig),
         # cluster = as.integer(gsub("cluster", "",sapply(strsplit(Contig, "_"), get_cluster))),
         ID = contig,
         "Name" = Best_Hit_ARO,
         "DB" = "CARD",
         "UniRef_function" = NA,
         "DB_taxa" = NA) %>%
  select(child,
         cluster,
         contig,
         ORF,
         Name,
         DB,
         ARO,
         `Drug Class`,
         `Resistance Mechanism`,
         UniRef_function,
         DB_taxa,
         ID
  )


## DoriC
doric_file <- "PATH/all_plasmids_DoriC.txt"
df_doric <- read_tsv(doric_file, col_names = FALSE)
names(df_doric) <- c("query", "subject", "qlength", "slength", "qcovs", "qcovhsp",
                     "identity", "alignment length", "mismatches", "gap opens", 
                     "qstart", "qend", "sstart", "send", "Evalue", 
                     "bit score", "STitel")

df_doric <- df_doric %>%
  mutate(child = sapply(strsplit(query, ".", fixed = T), get_child),
         contig = query)

# Add scovs column, percentage of subject/ARG covered
df_doric$scovs <- df_doric$`alignment length`/df_doric$slength * 100  

# make dataframe with only DoriC hits, filtered on 95% cov and 90% ident.
df_doric_filter <- df_doric %>% filter(scovs >= 95 & identity >= 90) 

## Some DoriC hits are cut in half due to the linearization of the plasmid
## calculate the approximation of the coverage of these:
# Max coverage of query/subject pair at the start of the contig
df_doric_begin <- df_doric %>% filter(qstart == 1) %>% 
  ungroup() %>%
  group_by(query, subject) %>%
  summarise(maxcov = max(scovs), .groups = "keep")
# Max coverage of the query/contig pair at the end of the contig
df_doric_end <- df_doric %>% filter(qend == qlength)  %>% 
  group_by(query, subject) %>%
  summarise(maxcov = max(scovs), .groups = "keep")
# add together the max coverage of the end and start of the query/subject pair
df_doric_sum <- rbind(df_doric_begin, df_doric_end) %>% 
  group_by(query, subject) %>% 
  summarise(scovs = sum(maxcov), .groups = "keep") %>% ungroup() %>% group_by(query)
# take the subject with highest scovs per query 
df_doric_sum <- top_n(df_doric_sum, 1, scovs) %>% top_n(1, subject)

# filter to only contain contigs that have a ori with >= 90 % coverage
df_doric_Sum_filter <- df_doric_sum %>% 
  filter(scovs >= 95) %>% 
  filter(!query %in% df_doric_filter$query) %>%
  select(query, subject) %>% 
  unique() 
df_doric_Sum_filter <- left_join(df_doric_Sum_filter, df_doric)

# Add all the ori hits into one dataframe
df_doric_filter_total <- rbind(df_doric_Sum_filter, df_doric_filter) 
## deduplicate
df_doric_dedup <- df_doric_filter_total %>% 
  group_by(query, qstart, child) %>% 
  top_n(-1, Evalue) %>%
  top_n(1, identity) %>%
  top_n(1, `bit score`) %>%
  top_n(1, subject) %>%
  ungroup() %>%
  group_by(query, qend, child) %>%
  top_n(-1, Evalue) %>%
  top_n(1, identity) %>%
  top_n(1, `bit score`) %>%
  top_n(1, subject) %>%
  ungroup()

## Load metadata DoriC
doric_meta_file <- "PATH/doric10/Plasmids.csv"
df_doric_metadata <- read_delim(doric_meta_file, delim = ";")
df_doric_metadata <- df_doric_metadata %>%
  mutate(subject = paste(doricAC,Refseq, sep = "_"))

df_doric_taxa <- left_join(df_doric_dedup, df_doric_metadata)

df_doric_total <- df_doric_taxa %>% 
  mutate(cluster = sapply(strsplit(query, "_"), get_cluster),
         "ORF" = NA,
         "Name" = subject,
         "DB" = "DoriC",
         "ARO" = NA,
         "Drug Class" = NA,
         "Resistance Mechanism" = NA,
         UniRef_function = NA,
         DB_taxa = sapply(strsplit(Organism, " "), get_taxa_dorit),
         ID = contig,
         start = qstart,
         end = qend,
         Orientation = ifelse(send > sstart, "+", "-"),
         info = NA,
         partial = FALSE) %>%
  select(child,
         cluster,
         contig,
         ORF,
         Name,
         DB,
         ARO,
         `Drug Class`,
         `Resistance Mechanism`,
         UniRef_function,
         DB_taxa,
         ID,
         start,
         end,
         Orientation,
         info,
         partial)


## OriT
orit_file <- "PATH/all_plasmids_oriT.txt"
df_orit <- read_tsv(orit_file, col_names = FALSE)
names(df_orit) <- c("query", "subject", "qlength", "slength", "qcovs", "qcovhsp",
                     "identity", "alignment length", "mismatches", "gap opens", 
                     "qstart", "qend", "sstart", "send", "Evalue", 
                     "bit score", "STitel")

df_orit <- df_orit %>%
  mutate(child = sapply(strsplit(query, ".", fixed = T), get_child),
         contig = query)

# Add scovs column, percentage of subject/ARG covered
df_orit$scovs <- df_orit$`alignment length`/df_orit$slength * 100  

# make dataframe with only orit hits, filtered on 95% cov and 90% ident.
df_orit_filter <- df_orit %>% filter(scovs >= 95 & identity >= 90) 

df_orit_total <- df_orit_filter %>% 
  mutate(cluster = sapply(strsplit(query, "_"), get_cluster),
         "ORF" = NA,
         "Name" = subject,
         "DB" = "OriT",
         "ARO" = NA,
         "Drug Class" = NA,
         "Resistance Mechanism" = NA,
         UniRef_function = NA,
         DB_taxa = sapply(strsplit(subject, "_"), get_taxa_orit),
         ID = contig,
         start = qstart,
         end = qend,
         Orientation = ifelse(send > sstart, "+", "-"),
         info = NA,
         partial = FALSE) %>%
  select(child,
         cluster,
         contig,
         ORF,
         Name,
         DB,
         ARO,
         `Drug Class`,
         `Resistance Mechanism`,
         UniRef_function,
         DB_taxa,
         ID,
         start,
         end,
         Orientation,
         info,
         partial)


## uniprot90
# create empty dataframe
df_uniprot <- data.frame(stringsAsFactors = FALSE)
uniprot_files <- list.files(path = "PATH/diamond/uniref90/", 
                            full.names = TRUE)
for (i in 1:length(uniprot_files)) {
  file_name <- tail(strsplit(uniprot_files[i], "/")[[1]], n=1)
  file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
  assign(file_name, read_tsv(uniprot_files[i], col_names = FALSE))
  row <- get(file_name)
  tmp_child <- head(unlist(strsplit(file_name, "_")), n=1)
  row <- row %>%
    mutate(
      child = toupper(tmp_child)
    )
  df_uniprot <- rbind(df_uniprot, row)
}
names(df_uniprot) <- c("query", "subject", "qlength", "slength", "qcovhsp",
                       "identity", "alignment length", "mismatches", 
                       "gap opens", "qstart", "qend", "sstart", "send", 
                       "Evalue","bit score", "STitel", "child")

df_uniprot_dedup <- df_uniprot %>%
  group_by(query, child) %>%
  top_n(-1, Evalue) %>%
  top_n(1, `alignment length`) %>%
  top_n(1, identity) %>%
  # top_n(1, subject) %>%
  ungroup()

df_uniprot_filter <- df_uniprot_dedup %>%
  filter(qcovhsp >= 80) %>%
  # filter(scovs >= 80) %>%
  # filter(identity >= 80)
  filter(Evalue <= 0.001)

## Load taxid translation
file_taxid <- "PATH/kraken2/taxonomy_all_2025_04_23.tsv"
df_taxid <- read_tsv(file_taxid)

df_taxid_genus <- df_taxid %>%
  filter(Rank %in% c("species", 
                     "genus",
                     "strain",
                     "subspecies")) %>%
  ## filter out Rank = species that are actually unclassified 
  ## at genus/species level
  filter(!grepl(" bacterium", `Scientific name`)) %>%
  filter(!grepl("uncultured", `Scientific name`, ignore.case = T)) %>%
  filter(!grepl("unculturable", `Scientific name`, ignore.case = T)) %>%
  filter(!grepl("uncultivated", `Scientific name`, ignore.case = T)) %>%
  ## remove brackets []
  mutate(`Scientific name` = gsub("[", "", `Scientific name`, fixed = T),
         `Scientific name` = gsub("]", "", `Scientific name`, fixed = T))


df_uniprot_tax <- df_uniprot_filter %>%
  mutate("Taxon Id" = sapply(strsplit(STitel, "TaxID="), get_taxid))
df_uniprot_tax <- left_join(df_uniprot_tax, df_taxid_genus)


df_uniprot_total <- df_uniprot_tax %>%
  mutate(UniRef_function = str_replace_all(sapply(strsplit(STitel, " n="), get_function_uniprot), ",", ";"),
         DB_taxa = `Scientific name`,
         cluster = sapply(strsplit(query, "_"), get_cluster),
         ORF = query,
         contig = sapply(strsplit(ORF, "_"), get_contig),
         Name = subject,
         DB = "UniRef90",
         ARO = subject,
         "Drug Class" = NA,
         "Resistance Mechanism" = NA,
         ID = contig) %>%
  select(child,
         cluster, 
         contig, 
         ORF, 
         Name,
         DB,
         ARO,
         `Drug Class`,
         `Resistance Mechanism`,
         UniRef_function,
         DB_taxa,
         ID)

## gather
df_total <- rbind(df_rgi_total, df_uniprot_total)
## add start, end, orientation etc from prodigal
df_total <- left_join(df_total, df_prodigal_select)
## fill in "empty" ORFs
df_prodigal_filter <- df_prodigal_total %>%
  filter(!ORF %in% df_total$ORF)
df_total <- rbind(df_total, df_prodigal_filter)
## add DoriC and OriT
df_total <- rbind(df_total, df_orit_total)
df_total <- rbind(df_total, df_doric_total)

df_total <- df_total %>%
  unique() %>%
  arrange(child, cluster, start) %>%
  ## remove ,
  mutate(UniRef_function = gsub(",", ";", UniRef_function, fixed = TRUE))

## save annotated
outfile <- "plasmids_annotation_clusters.csv"
write.csv(df_total, file = outfile, row.names = FALSE, quote = FALSE)


