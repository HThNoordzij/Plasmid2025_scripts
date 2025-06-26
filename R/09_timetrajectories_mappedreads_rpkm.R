# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")
library(RColorBrewer)
# display.brewer.all()
palette<-brewer.pal(n=12,name='Paired')

# set working directory
setwd("~/PATH/")

## functions
get_child <- function(x) {
  x <- unlist(x)[2]
  x <- as.integer(gsub("id", "", x, ignore.case = T))
  return(x)
}
get_cluster <- function(x) {
  x <- unlist(x)[1]
  x <- as.integer(gsub("cluster", "", x, ignore.case = T))
  return(x)
}

## Load abundance files
abundance_file <- "abundance_summary_clusters.csv"
df_abundance <- read_csv(abundance_file)

## load read counts
read_files <- list.files(
  path = "PATH/",
  full.names = T
)

# create empty dataframe
df_reads <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(read_files)) {
  file_name <- tail(strsplit(read_files[i], "/")[[1]], n=1)
  file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
  
  assign(file_name, read.table(read_files[i],
                               header = TRUE))
  
  row <- get(file_name)
  childno <-as.integer(gsub("ID", "", 
                          head(unlist(strsplit(file_name, "_")), n=1), 
                          ignore.case = T))  
  row$child <- childno

  row <- row %>%
    mutate(day = as.integer(gsub("d", "", day))) %>%
    select(reads, day, child)
  
  df_reads <- rbind(df_reads, row)
}

df_reads_millions <- df_reads %>%
  mutate(million = reads / 1000000)

## Load flagstat files
flagstat_files <- list.files(
  path = "PATH/flagstat/", 
  full.names = TRUE)

# create empty dataframe
df_flagstat <- data.frame(stringsAsFactors = FALSE)

## since I only keep mapped reads in bam file
## need to caluclate % mapped with read counts

for (i in 1:length(flagstat_files)) {
  file_name <- tail(strsplit(flagstat_files[i], "/")[[1]], n=1)
  file_name <- head(strsplit(file_name, "[.]")[[1]], n=1)
  assign(file_name, read.table(flagstat_files[i], 
                               header = FALSE,
                               sep = " "))
  tmp_df <- get(file_name)
  row <- list()
  row$mapped <- tmp_df[1,1]
  
  child <- head(unlist(strsplit(file_name, "_")), n=1)
  row$child <- child
  
  childno <- as.integer(gsub("id", "", child, ignore.case = T))
  row$child <- childno
  
  day <- tail(head(unlist(strsplit(file_name, "_")), n=2), n=1)
  day <- gsub("d", "", day)
  row$day <- as.integer(day)
  
  df_flagstat <- rbind(df_flagstat, row)
}

df_flagstat_reads <- left_join(df_flagstat, df_reads)

df_flagstat_reads <- df_flagstat_reads %>%
  filter(!is.na(reads)) %>%
  mutate(perc = as.integer(mapped) / as.integer(reads) * 100)

df_flagstat_reads <- df_flagstat_reads %>% 
  select(perc, child,day) %>%
  arrange(day) %>%
  group_by(child)

df_coverage <- df_abundance %>%
  group_by(child, day) %>%
  summarise(TotalRPKM = sum(RPKM)) %>%
  ungroup() %>%
  mutate(child = as.integer(gsub("id", "", child))) %>%
  arrange(day) %>%
  group_by(child)


############# plot
colors <- palette
grDevices::windows(10,10)
par(mfrow=c(3,4))
par(oma=c(2,2,0,2))
par(mar=c(2,2,2,3))

for (i in 1:12) {
  df_plot_flagstat <- df_flagstat_reads %>%
    filter(child == i) %>%
    as.data.frame()
  
  df_plot_RPKM <- df_coverage %>%
    filter(child == i) %>%
    as.data.frame()
  
  tmp_days <- df_reads %>%
    filter(child == i) %>%
    select(day) %>%
    unique() %>%
    unlist() %>%
    as.character() %>%
    as.integer()
  
  tmp_title <- paste0("ID",i)
  
  ## add missing days
  for (j in 1:length(tmp_days)) {
    tmp_day <- tmp_days[[j]]
    if(nrow(df_plot_flagstat %>% filter(day == tmp_day)) == 0){
      tmp_row <- c(as.double(0), i, tmp_day)
      df_plot_flagstat <- rbind(df_plot_flagstat, tmp_row)
    }
    if(nrow(df_plot_RPKM %>% filter(day == tmp_day)) == 0){
      tmp_row <- c(as.double(0), i, tmp_day)
      df_plot_RPKM <- rbind(df_plot_RPKM, tmp_row)
    }
  }
  df_plot_RPKM <- df_plot_RPKM %>% 
    arrange(day) 
  df_plot_flagstat <- df_plot_flagstat %>% 
    arrange(day) 
  
  # RPKM plot
  plot(df_plot_RPKM$day,
       df_plot_RPKM$TotalRPKM,
       pch = 20,
       # col=colors[df_plot_RPKM$child],
       col = "grey",
       lwd = 2,
       main = tmp_title,
       axes=FALSE,
       xlab="", 
       ylab="")
  # RPKM line
  lines(df_plot_RPKM$day,
        df_plot_RPKM$TotalRPKM,
        # col=colors[df_plot_RPKM$child]
        col = "grey",
        lwd = 2,)
  
  axis(4,
       pretty(range(0,max(df_plot_RPKM$TotalRPKM))),
       col="black")
  box()
  
  ## legend
  if(i == 1) {
    legend('topright',
           legend = c("Mapping", "RPKM"),
           fill = c("black", "grey"))
  }
  
  ## add second plot
  par(new = TRUE)
  
  #flagstat plot with dots
  plot(df_plot_flagstat$day,
       df_plot_flagstat$perc,
       pch = 20,
       # col=colors[df_plot_flagstat$child],
       col = "black",
       lwd = 2,
       axes=FALSE, 
       ylim=c(0,max(df_plot_flagstat$perc)),
       xlab="",
       ylab="")
  # flagstat line
  lines(df_plot_flagstat$day,
        df_plot_flagstat$perc,
        # col=colors[df_plot_flagstat$child]
        col = "black",
        lwd = 2,
  )
  axis(2, 
       pretty(range(df_plot_flagstat$perc)),
       col="black")
  axis(1, 
       c(0, 100, 200, 300, 400),
       col="black")
  
  if(i == 12) {
    mtext("% mapped",
          side = 2, 
          line = 0.8,
          outer=TRUE)
    mtext("RPKM",
          side = 4, 
          line = 0.8,
          outer=TRUE, 
          las=0)
    mtext("Days",
          side = 1, 
          line = 0.8,
          outer=TRUE)
  }
}

