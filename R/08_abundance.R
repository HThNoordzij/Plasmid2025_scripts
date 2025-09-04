
# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

# Load libraries
library("tidyverse")
library(data.table)

# Set working directory
setwd("~/PATH/")

# functions
fread_with_arguments <- function(file) {
  return(fread(file, 
               sep = " ",
               skip = 1,
               col.names = c("reads", 
                             "day",
                             "sample"),
               fill = TRUE))
}
get_child <- function(x) {
  x <- head(unlist(x), n=1)
  return(x)
}
get_day <- function(x) {
  x <- tail(unlist(x), n=1)
  x <- as.integer(gsub("d", "", x))
  return(x)
}
get_cluster <- function(x) {
  x <- tail(unlist(x), n=1)
  x <- as.integer(gsub("cluster", "", x))
  return(x)
}

#################################   infants   ##################################
## Load reads file
file_name_reads <- list.files(path = "PATH/", full.names = TRUE)
names(file_name_reads) <- tolower(str_remove(basename(file_name_reads),
                                             "\\_depth.txt"))
df_reads <- map_dfr(file_name_reads, 
                    fread_with_arguments, .id = 'child')

df_reads <- df_reads %>%
  mutate(day = as.integer(gsub("d", "", day))) %>%
  select(child,day, reads)

## Load coverage files
coverage_files <- list.files(path = "bwa/coverage/", full.names = TRUE)
names(coverage_files) <- str_remove(basename(coverage_files),
                                     "\\_plasmids_coverage.txt")
df_coverage <- map_dfr(coverage_files, 
                    fread, .id = 'ID')
df_coverage <- df_coverage %>%
  mutate(child = sapply(strsplit(ID, "_", fixed = T), get_child),
         day = sapply(strsplit(ID, "_", fixed = T), get_day),
         rname = `#rname`,
         cluster = sapply(strsplit(rname, ".", fixed = T), get_cluster)) %>%
  select(-`#rname`)

## add child, cluster, day
## Rename cluster, see 07_macsyfinder.R
##  so different infants with the same plasmid have the same cluster number
file_rename <- "rename_plasmids.csv"
df_rename <- fread(file_rename,
                   sep = ",",
                   fill = T)

df_coverage <- df_coverage %>%
  mutate(contig = sapply(rname, function(x) {
    return(df_rename$contig_new[which(
      df_rename$contig_old == x)])}),
    cluster = sapply(strsplit(contig, ".", fixed = TRUE), get_cluster)
  )
df_coverage_reads <- left_join(df_coverage, df_reads)

## add RPKM
df_coverage_reads <- df_coverage_reads %>%
  mutate(RPKM = (numreads / 
                   (as.double(endpos) * as.double(reads))
  ) * 1000000000 )

## filter on 80% coverage and meandepth 10
df_coverage_filter <- df_coverage_reads %>%
  filter(coverage >= 80) %>% filter(meandepth >= 10)

## Safe clusters abundance summary
outfile <- "abundance_summary_clusters.csv"
fwrite(df_coverage_filter,
      file = outfile,
      sep = ",",
      row.names = F,
      col.names = T)
