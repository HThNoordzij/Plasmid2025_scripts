
# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

# Load libraries
library("tidyverse")

# Set working directory
setwd("~/PATH/")

df_annotation_clusters <- data.frame()

for (i in 1:12) {
  child_no <- paste("id", i, sep = "")
  
  ### Open metaphlan
  file_metaphlan <- paste("metaphlan/",
                          child_no, 
                          "_metaphlan_merged.txt",sep = "") 
  df_metaphlan <- read_tsv(file_metaphlan, comment = "#")
  df_species <- df_metaphlan %>%
    filter(grepl("s__", clade_name)) %>%
    mutate(species = sapply(strsplit(clade_name, "s__"), function(x)x[2])) %>%
    select(-clade_name, -NCBI_tax_id) %>%
    as.matrix()
  species_names <- df_species[,ncol(df_species)]
  df_species <- df_species[,1:ncol(df_species)-1]
  sample_names <- sapply(colnames(df_species), 
                         function(x) {
                           gsub(paste(child_no,"_d", sep=""), "",
                                head(unlist(strsplit(x, "_s")),n=1))
                         } )
  df_species_matrix <- matrix(as.double(df_species), 
                              nrow = nrow(df_species), 
                              ncol = ncol(df_species))
  colnames(df_species_matrix) <- sample_names
  rownames(df_species_matrix) <- species_names
  df_species_matrix<-df_species_matrix[,order(as.numeric(colnames(df_species_matrix)))]
  df_species_matrix <- df_species_matrix[apply(df_species_matrix,1,function(x)sum(as.double(x)>0))>=5,]
    
  ## Load summary file for abundance
  file_abundance_summary <- "abundance_summary_clusters.csv"
  df_abundance_summary <- read.table(file_abundance_summary, sep = ",", header = T)
  df_abundance_summary <- df_abundance_summary %>% 
    filter(child == child_no) %>%
    select(cluster, day, RPKM) 
  
  days <- df_abundance_summary %>% select(day) %>% unlist() %>% unique() %>% as.numeric()
  
  df_abundance_child <- df_abundance_summary %>%
    pivot_wider(names_from = day, values_from = RPKM, values_fill = 0)
  
  df_abundance_child <- as.matrix(df_abundance_child)
  child_no_names <- df_abundance_child[,1]
  df_abundance_child <- df_abundance_child[,2:ncol(df_abundance_child)]
  rownames(df_abundance_child) <- child_no_names
  child_samples <- colnames(df_abundance_child)
  df_abundance_child<-df_abundance_child[,order(as.numeric(child_samples))]
  child_samples <- colnames(df_abundance_child)
  
  df_abundance_child <- df_abundance_child[apply(df_abundance_child,1,function(x)sum(x>0))>=5,]
  child_allcor<-cor(t(df_abundance_child),t(df_species_matrix))
  
  dim(child_allcor)
  hist(apply(child_allcor,1,max),breaks=10)
  child_sets<-list()
  for(i in 1:ncol(child_allcor)){
    temp<-which(child_allcor[,i]>.95)
    child_sets[[i]]<-temp
  }
    
  names(child_sets)<-colnames(child_allcor)
  child_sets<-child_sets[lapply(child_sets,length)>0]
  length(child_sets)
  names(child_sets)
  
  if(length(child_sets) < 3) {
    tiff(filename = paste("Figures/correlations_species_",
                          child_no, 
                          ".tif", sep = ""),
         width = 20,
         height = 10,
         units = "in",
         res = 100
    )
  } else {
    tiff(filename = paste("Figures/correlations_species_",
                          child_no, 
                          ".tif", sep = ""),
         width = 20,
         height = 20,
         units = "in",
         res = 300
    )
  }
  
  ## Plot all
  for (i in 1:length(child_sets)) {
    # set margins etc
    if(i == 1) {
      plots <- length(child_sets)
      rows <- ceiling(plots / 2)
      par(mfrow=c(rows,2))
      par(oma=c(2,2,0,2))
      par(mar=c(2,2,2,3))
    }
    
    tmp_species <- names(child_sets[i])
    tmp_species_abundance <- df_species_matrix[tmp_species,]
    tmp_cluster_abundance <- df_abundance_child[child_sets[[i]],, drop = FALSE]
    rownames(tmp_cluster_abundance) <- gsub("tmp_","", rownames(tmp_cluster_abundance))
    ifelse(is.null(rownames(tmp_cluster_abundance)),
           tmp_clusters <- "23",
           tmp_clusters <- rownames(tmp_cluster_abundance))
    
    plot(as.integer(names(tmp_species_abundance)),
         tmp_species_abundance,
         type = "l",
         main = toupper(child_no),
         col = "black",
         lwd = 5,
         axes=FALSE, 
         xlab="",
         ylab="",
         cex.main = 1.5,
         cex.axis = 1.5)
    box()
    axis(2, 
         pretty(range(tmp_species_abundance)),
         col="black",
         cex.axis = 1.5)
    axis(1, 
         c(0, 100, 200, 300, 400),
         col="black",
         cex.axis = 1.5)
    
    ## extra plot and lines per cluster
    for (j in 1:ifelse(is.null(nrow(tmp_cluster_abundance)),1,nrow(tmp_cluster_abundance))) {
      if(j == 1) {
        
        ifelse(is.null(nrow(tmp_cluster_abundance)),
               tmp_cluster_abundance_first <- tmp_cluster_abundance,
               tmp_cluster_abundance_first <- tmp_cluster_abundance[j,])
        
        par(new = TRUE)
        plot(as.integer(names(tmp_cluster_abundance_first)),
             tmp_cluster_abundance_first,
             type = "l",
             col = 2,
             lwd = 5,
             lty=2,
             axes=FALSE, 
             ylim = c(0,max(tmp_cluster_abundance)),
             xlab="",
             ylab="")
        axis(4,
             pretty(range(0,max(tmp_cluster_abundance))),
             col="black",
             cex.axis = 1.5)
      } else {
        tmp_cluster_abundance_later <- tmp_cluster_abundance[j,]
        lines(as.integer(names(tmp_cluster_abundance_later)),
              tmp_cluster_abundance_later,
              col = 1+j,
              lwd = 5,
              lty=2)
      }
      
      if(j == ifelse(is.null(nrow(tmp_cluster_abundance)), 
                     1,
                     nrow(tmp_cluster_abundance)))
        legend('topleft',
               legend = c(gsub("_", " ", 
                               gsub("complex", "", tmp_species)),
                          paste0(tmp_clusters,
                                 " (",
                                 round(child_allcor[tmp_clusters,tmp_species], digits = 3),
                                 ")")),
               fill = c(1,2:(ifelse(is.null(nrow(tmp_cluster_abundance)), 
                                    1,
                                    nrow(tmp_cluster_abundance))+1)),
               title = "cluster (Pearson):",
               cex = 1.5)
    }
    
    # if(i == length(child_sets)) {
    mtext("Species relative abundance",
          side = 2, 
          line = 0.8,
          outer=TRUE,
          cex = 1.5)
    mtext("Plasmid RPKM",
          side = 4, 
          line = 0.8,
          outer=TRUE,
          cex = 1.5, 
          las=0)
    mtext("Days",
          side = 1, 
          line = 0.8,
          cex = 1.5,
          outer=TRUE)
    # }
    
  }
  dev.off()
  
  ###### end ######
  ############ Add species to annotation summary ##################
  
  df_species <- data.frame(child = as.character(),
                           contig = as.character(),
                           species = as.character(),
                           species_corr_number = as.character(),
                           stringsAsFactors = FALSE)
  species <- names(child_sets)  
  
  for (i in 1:length(species)) {
    s <- species[i]
    contigs <- names(child_sets[[i]])
    for (j in 1:length(contigs)) {
      contig <- contigs[j]
      correlate_number <- round(child_allcor[rownames(child_allcor)==contig,colnames(child_allcor)==s], digits = 4)
      
      if(contig %in% df_species[[2]]){
        add_species <- paste0(df_species[df_species[[2]] == contig, ][[3]], ";", s)
        add_species_corr_numb <- paste0(df_species[df_species[[2]] == contig, ][[4]], ";", correlate_number)
        df_species[df_species[[2]] == contig, ][[3]] <- add_species
        df_species[df_species[[2]] == contig, ][[4]] <- add_species_corr_numb
      } else {
        new_row <- c(child_no,contig, s, correlate_number)
        df_species <- rbind(df_species, new_row)
      }
    }
  }
  colnames(df_species) <- c("child", "cluster", "species", "species_corr_number")
  df_species$cluster <- as.integer(df_species$cluster)
  df_species$species <- gsub("_", " ", df_species$species)
  df_species
  
  # load annotaion/summary so far
  file_mash <- "plasmids_tax_mobility_GC_clusters.csv"
  df_mash <- read_csv(file_mash)
  
  df_annotation_child <- df_mash %>% 
    filter(child == toupper(child_no)) %>%
    mutate(cluster = as.integer(cluster),
           child = tolower(child)
    )
  
  df_annotation_species <- left_join(df_annotation_child, df_species)
  
  ### Correlate clusters with clusters ####
  cluster_corr <- cor(t(df_abundance_child))
  dim(cluster_corr)
  hist(apply(cluster_corr,1,max),breaks=10)
  child_sets<-list()
  for(i in 1:ncol(cluster_corr)){
    temp <- cluster_corr[,i] 
    temp<-which(temp>.95)
    child_sets[[i]]<-temp
  }
  names(child_sets)<-colnames(cluster_corr)
  child_sets<-child_sets[lapply(child_sets,length)>1]
  length(child_sets)
  names(child_sets)
  
  ###### add correlating clusters to summary dataframe
  df_clusters_corr <- data.frame(child = as.character(),
                                 contig = as.character(),
                                 clusters_corr = as.character(),
                                 clusters_corr_number = as.character(),
                                 stringsAsFactors = FALSE)
  
  clusters_corr <- names(child_sets)
  child_sets
  
  df_clusters_corr <- data.frame(child = as.character(),
                                 cluster = as.integer(),
                                 clusters_corr = as.character(),
                                 clusters_corr_number = as.character(),
                                 stringsAsFactors = F)
  
  if(length(clusters_corr) > 0) {
    for (i in 1:length(clusters_corr)) {
      s <- clusters_corr[i]
      contigs <- names(child_sets[[i]])
      
      for (j in 1:length(contigs)) {
        contig <- contigs[j]
        if(contig != s) {
          correlate_number <- round(cluster_corr[rownames(cluster_corr)==s,colnames(cluster_corr)==contig], digits = 4)
          
          if(contig %in% df_clusters_corr[[2]] &
             nchar(contig) > 0){
            add_clusters_corr <- paste0(df_clusters_corr[df_clusters_corr[[2]] == contig, ][[3]], ";", s)
            add_clusters_corr_numb <- paste0(df_clusters_corr[df_clusters_corr[[2]] == contig, ][[4]], ";", correlate_number)
            
            df_clusters_corr[df_clusters_corr[[2]] == contig, ][[3]] <- add_clusters_corr
            df_clusters_corr[df_clusters_corr[[2]] == contig, ][[4]] <- add_clusters_corr_numb
          } else if(nchar(contig) > 0) {
            new_row <- c(child_no,contig, s, correlate_number)
            df_clusters_corr <- rbind(df_clusters_corr, new_row)
          }
        }
      }
    }
    colnames(df_clusters_corr) <- c("child", "cluster", "clusters_corr", "clusters_corr_number")
    df_clusters_corr$cluster <- as.integer(df_clusters_corr$cluster)
  }
    
  df_annotation_species <- left_join(df_annotation_child, df_species)
  df_annotation_clusters_child <- left_join(df_annotation_species, df_clusters_corr)
  
  df_annotation_clusters <- rbind(df_annotation_clusters,
                                  df_annotation_clusters_child)
  
}

##### Save new summary
df_annotation_clusters <- df_annotation_clusters %>%
  mutate(child = toupper(child)) %>%
  arrange(child, cluster)

summary_outfile <- "cluster_species_correlation.csv"
write_csv(df_annotation_clusters, file = summary_outfile)

