
# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

# Load libraries
library("tidyverse")

# Set working directory
setwd("~/PATH/")

# functions
get_cluster <- function(x) {
  x <- tail(unlist(x), n=1)
  x <- as.integer(gsub("cluster", "", x))
  return(x)
}

#################################   infants   ##################################

## Load reads file
# empty dataframe
df_reads <- data.frame()
file_name_reads <- list.files(path = "PATH/", full.names = TRUE)
for (i in 1:length(file_name_reads)) {
  file_name <- tail(strsplit(file_name_reads[i], "/")[[1]], n=1)
  file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
  assign(file_name, read.table(file_name_reads[i], sep = "",header = TRUE, 
                               col.names = c("reads", "day", "sample"),
                               fill = TRUE))
  
  tmp_child <- tolower(head(unlist(strsplit(file_name, "_")), n=1))
  tmp_df <- get(file_name)
  tmp_df <- tmp_df %>%
    mutate(day = as.integer(gsub("d", "", day)),
           child = tmp_child) %>%
    select(child,day, reads)
  df_reads <- rbind(df_reads, tmp_df)
}

## Load coverage files
# empty dataframe
df_coverage <- data.frame()
coverage_files <- list.files(path = "bwa/coverage/", full.names = TRUE)
for (i in 1:length(coverage_files)) {
  file_name <- tail(strsplit(coverage_files[i], "/")[[1]], n=1)
  file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
  assign(file_name, read.table(coverage_files[i], 
                               col.names = c("rname","startpos","endpos",
                                             "numreads","covbases","coverage",
                                             "meandepth","meanbaseq","meanmapq")))
  tmp_child <- head(unlist(strsplit(file_name, "_")), n=1)
  tmp_day <-as.integer(gsub("d", "", tail(head(unlist(strsplit(file_name, "_")),n=2),n=1)))  
  tmp_df <- get(file_name)
  tmp_df <- tmp_df %>%
    mutate(child = tmp_child,
           day = tmp_day)
  df_coverage <- rbind(df_coverage, tmp_df)
}

## add child, cluster, day
## Rename cluster, see 07_macsyfinder.R
##  so different infants with the same plasmid have the same cluster number
file_rename <- "rename_plasmids.csv"
df_rename <- read_csv(file_rename)

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
write.csv(df_coverage_filter,
          file = outfile, row.names = FALSE, quote = FALSE)
