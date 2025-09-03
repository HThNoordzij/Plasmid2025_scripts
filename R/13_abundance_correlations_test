
# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

# Load libraries
library("tidyverse")
library(data.table)

# Set working directory
setwd("~/PATH/")

## Load summary file for abundance
file_abundance_summary <- "abundance_summary_clusters.csv"
df_abundance_summary <- fread(file_abundance_summary, 
                              sep = ",",
                              fill = TRUE)

## Loop through infants and gather correlations
df_correlations <- data.frame()
for (j in 1:12) {
  
  child_abundance <- paste("id", j, sep = "")
  
  df_abundance_child <- df_abundance_summary %>% 
    filter(child == child_abundance) %>%
    select(cluster, day, RPKM) 
  
  days <- df_abundance_child %>% 
    select(day) %>% 
    unlist() %>% 
    unique() %>% 
    as.numeric()
  
  df_abundance_pivot <- df_abundance_child %>%
    pivot_wider(names_from = day, values_from = RPKM, values_fill = 0)
  
  df_abundance_child <- as.matrix(df_abundance_pivot)
  child_no_names <- df_abundance_child[,1]
  df_abundance_child <- df_abundance_child[,2:ncol(df_abundance_child)]
  rownames(df_abundance_child) <- child_no_names
  child_samples <- colnames(df_abundance_child)
  df_abundance_child<-df_abundance_child[,order(as.numeric(child_samples))]
  child_samples <- colnames(df_abundance_child)
  
  df_abundance_child <- df_abundance_child[apply(df_abundance_child,
                                                 1,
                                                 function(x)sum(x>0))>=5,]
  
  ### make correlations with other infants
    for (i in 1:12) {
    
    child_no <- paste("id", i, sep = "")
    
    ### Load metaphlan
    file_metaphlan <- paste("metaphlan/",
                            child_no, 
                            "_metaphlan_merged.txt",sep = "") 
    df_metaphlan <- fread(file_metaphlan, 
                          skip = 1, 
                          fill = TRUE)
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
    
    ## select equal number of days
    df_abundance_child_equal <- df_abundance_child
    df_species_matrix_equal <- df_species_matrix
    if(ncol(df_abundance_child) > ncol(df_species_matrix)) {
      df_abundance_child_equal <- df_abundance_child[,1:ncol(df_species_matrix)]
    } else if (ncol(df_abundance_child) < ncol(df_species_matrix)) {
      df_species_matrix_equal <- df_species_matrix[,1:ncol(df_abundance_child)]
    }
    
    ## Correlate abundance plasmids with species abundance
    child_allcor<-cor(t(df_abundance_child_equal),t(df_species_matrix_equal))
    
    child_sets<-list()
    for(i in 1:ncol(child_allcor)){
      temp<-which(child_allcor[,i]>.95)
      child_sets[[i]]<-temp
    }
    
    names(child_sets)<-colnames(child_allcor)
    child_sets<-child_sets[lapply(child_sets,length)>0]
    
    tmp_row <- c(child_abundance, 
                 child_no, 
                 length(child_sets),
                 ncol(child_allcor) * nrow(child_allcor))
    
    df_correlations <- rbind(df_correlations,
                             tmp_row)
  }
}

colnames(df_correlations) <- c("Plasmid", 
                               "Species", 
                               "correlations",
                               "total_comparisons")
df_correlations <- df_correlations %>%
  mutate(not_correlated = as.integer(total_comparisons) - 
           as.integer(correlations)) %>%
  select(-total_comparisons)


df_correlations_within <- df_correlations %>%
  filter(Plasmid == Species) %>%
  mutate(observerd = as.integer(correlations)) 

df_correlations_between <- df_correlations %>%
  filter(Plasmid != Species) 

df_correlations_exp <- df_correlations_between %>%
  select(-Species) %>%
  group_by(Plasmid) %>%
  mutate(expected = mean(as.integer(correlations))) %>%
  select(Plasmid, expected) %>%
  unique()

df_correlations_obs_ex <- left_join(df_correlations_within,
                                    df_correlations_exp) %>%
  select(-Plasmid)

## test chi squared
df_chi <- df_correlations_obs_ex %>%
  mutate(o_min_e_square = (observerd - expected)^2,
         div_e = o_min_e_square / expected)

file_chi <- "correlations_chisquared.csv"
write_csv(df_chi, file = file_chi)
