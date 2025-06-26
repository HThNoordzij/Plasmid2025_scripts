## Get GC content per plasmid

# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

#Load package
library("tidyverse")
library("Biostrings")

library(RColorBrewer)
# display.brewer.all()
palette<-brewer.pal(n=8, name='Set1')

# Set working directory
setwd("~/PATH/")

## functions
get_contig <- function(x) {
  x <- unlist(x)
  x <- paste(x[2], x[3], x[4], x[5], sep = "_")
  return(x)
}
get_length <- function(x) {
  x <- as.integer(tail(unlist(x), n=1))
  return(x/1000)
}
get_child <- function(x) {
  x <- tail(unlist(x), n=1)
  x <- head(unlist(strsplit(x, "_")),n=1)
  return(x)
}


# Load cluster summary
file_clusters <- "plasmids_tax_mobility_clusters.csv"
df_clusters <- read_csv(file_clusters)

file_fasta <- "plasmids.fasta"
df_fasta <- readAAStringSet(file_fasta)

df_GC <- data.frame(
  "contig" = names(df_fasta),
  "GC" = letterFrequency(df_fasta, letters = "GC") / letterFrequency(df_fasta, letters = "ATGC") * 100
) %>%
  mutate("representative contig child" = contig,
         GC = G.C,
         ) %>%
  select(`representative contig child`, GC)

df_clusters_GC <- left_join(df_clusters, df_GC)

df_clusters_summary <- df_clusters_GC %>%
  select(c(colnames(df_clusters)[1:4],
           GC,
           colnames(df_clusters)[5:22]))

## Save summary
file_out <- "plasmids_tax_mobility_GC_clusters.csv"
write.csv(df_clusters_summary, file_out,row.names = F)

df_bacteroides <- df_clusters_summary %>% 
  filter(`putative host` == "Bacteroides") %>%
  mutate(genus = "Bacteroides")

df_bifido <- df_clusters_summary %>% 
  filter(`putative host` == "Bifidobacterium") %>%
  mutate(genus = "Bifidobacterium")

df_clostridium <- df_clusters_summary %>% 
  filter(`putative host` == "Clostridium") %>%
  mutate(genus = "Clostridium")

df_collinsella <- df_clusters_summary %>% 
  filter(`putative host` == "Collinsella") %>%
  mutate(genus = "Collinsella")

df_enterococcus <- df_clusters_summary %>% 
  filter(`putative host` == "Enterococcus") %>%
  mutate(genus = "Enterococcus")

df_escherichia <- df_clusters_summary %>% 
  filter(`putative host` == "Escherichia") %>%
  mutate(genus = "Escherichia")

df_klebsiella <- df_clusters_summary %>% 
  filter(`putative host` == "Klebsiella") %>%
  mutate(genus = "Klebsiella")

df_staphylococcus <- df_clusters_summary %>% 
  filter(`putative host` == "Staphylococcus") %>%
  mutate(genus = "Staphylococcus")

df_streptococcus <- df_clusters_summary %>% 
  filter(`putative host` == "Streptococcus") %>%
  mutate(genus = "Streptococcus")

df_Bacteroidales <- df_clusters_summary %>% 
  filter(`putative broad host range (pBHR) order` == "pBHR_Bacteroidales") %>%
  mutate(genus = "pBHR_Bacteroidales")

df_enterobacterales <- df_clusters_summary %>% 
  filter(`putative broad host range (pBHR) order` == "pBHR_Enterobacterales") %>%
  mutate(genus = "pBHR_Enterobacterales")

df_lactobacillales <- df_clusters_summary %>% 
  filter(`putative broad host range (pBHR) order` == "pBHR_Lactobacillales") %>%
  mutate(genus = "pBHR_Lactobacillales")

df_pbhrother <- df_clusters_summary %>% 
  filter(`putative broad host range (pBHR) order` == "pBHR_other") %>%
  mutate(genus = "pBHR_other")


df_bacteroides_sum <- data.frame(
  genus = "Bacteroides",
  length = mean(df_bacteroides$length),
  length_sd = sd(df_bacteroides$length),
  GC = mean(df_bacteroides$GC),
  GC_sd = sd(df_bacteroides$GC),
  n = nrow(df_bacteroides)
)
df_bifido_sum <- data.frame(
  genus = "Bifidobacterium",
  length = mean(df_bifido$length),
  length_sd = sd(df_bifido$length),
  GC = mean(df_bifido$GC),
  GC_sd = sd(df_bifido$GC),
  n = nrow(df_bifido)
)
df_clostridium_sum <- data.frame(
  genus = "Clostridium",
  length = mean(df_clostridium$length),
  length_sd = sd(df_clostridium$length),
  GC = mean(df_clostridium$GC),
  GC_sd = sd(df_clostridium$GC),
  n = nrow(df_clostridium)
)
df_collinsella_sum <- data.frame(
  genus = "Collinsella",
  length = mean(df_collinsella$length),
  length_sd = sd(df_collinsella$length),
  GC = mean(df_collinsella$GC),
  GC_sd = sd(df_collinsella$GC),
  n = nrow(df_collinsella)
)
df_enterococcus_sum <- data.frame(
  genus = "Enterococcus",
  length = mean(df_enterococcus$length),
  length_sd = sd(df_enterococcus$length),
  GC = mean(df_enterococcus$GC),
  GC_sd = sd(df_enterococcus$GC),
  n = nrow(df_enterococcus)
)
df_escherichia_sum <- data.frame(
  genus = "Escherichia",
  length = mean(df_escherichia$length),
  length_sd = sd(df_escherichia$length),
  GC = mean(df_escherichia$GC),
  GC_sd = sd(df_escherichia$GC),
  n = nrow(df_escherichia)
)
df_klebsiella_sum <- data.frame(
  genus = "Klebsiella",
  length = mean(df_klebsiella$length),
  length_sd = sd(df_klebsiella$length),
  GC = mean(df_klebsiella$GC),
  GC_sd = sd(df_klebsiella$GC),
  n = nrow(df_klebsiella)
)
df_staphylococcus_sum <- data.frame(
  genus = "Staphylococcus",
  length = mean(df_staphylococcus$length),
  length_sd = sd(df_staphylococcus$length),
  GC = mean(df_staphylococcus$GC),
  GC_sd = sd(df_staphylococcus$GC),
  n = nrow(df_staphylococcus)
)
df_streptococcus_sum <- data.frame(
  genus = "Streptococcus",
  length = mean(df_streptococcus$length),
  length_sd = sd(df_streptococcus$length),
  GC = mean(df_streptococcus$GC),
  GC_sd = sd(df_streptococcus$GC),
  n = nrow(df_streptococcus)
)
df_Bacteroidales_sum <- data.frame(
  genus = "pBHR_Bacteroidales",
  length = mean(df_Bacteroidales$length),
  length_sd = sd(df_Bacteroidales$length),
  GC = mean(df_Bacteroidales$GC),
  GC_sd = sd(df_Bacteroidales$GC),
  n = nrow(df_Bacteroidales)
)
df_enterobacterales_sum <- data.frame(
  genus = "pBHR_Enterobacterales",
  length = mean(df_enterobacterales$length),
  length_sd = sd(df_enterobacterales$length),
  GC = mean(df_enterobacterales$GC),
  GC_sd = sd(df_enterobacterales$GC),
  n = nrow(df_enterobacterales)
)
df_lactobacillales_sum <- data.frame(
  genus = "pBHR_Lactobacillales",
  length = mean(df_lactobacillales$length),
  length_sd = sd(df_lactobacillales$length),
  GC = mean(df_lactobacillales$GC),
  GC_sd = sd(df_lactobacillales$GC),
  n = nrow(df_lactobacillales)
)
df_pbhrother_sum <- data.frame(
  genus = "pBHR_other",
  length = mean(df_pbhrother$length),
  length_sd = sd(df_pbhrother$length),
  GC = mean(df_pbhrother$GC),
  GC_sd = sd(df_pbhrother$GC),
  n = nrow(df_pbhrother)
)

df_sum <- rbind(df_bacteroides_sum, df_bifido_sum)
df_sum <- rbind(df_sum, df_clostridium_sum)
df_sum <- rbind(df_sum, df_collinsella_sum)
df_sum <- rbind(df_sum, df_enterococcus_sum)
df_sum <- rbind(df_sum, df_escherichia_sum)
df_sum <- rbind(df_sum, df_klebsiella_sum)
df_sum <- rbind(df_sum, df_staphylococcus_sum)
df_sum <- rbind(df_sum, df_streptococcus_sum)
df_sum <- rbind(df_sum, df_Bacteroidales_sum)
df_sum <- rbind(df_sum, df_enterobacterales_sum)
df_sum <- rbind(df_sum, df_lactobacillales_sum)
df_sum <- rbind(df_sum, df_pbhrother_sum)

################################################################################
### GC of MAGs
file_gtdb <- list.files(path = "PATH/concoct/gtdbtk/",
                        full.names = TRUE)
df_gtdb <- read_tsv(file_gtdb, id = file_gtdb)
df_gtdb <- df_gtdb %>%
  mutate("Bin Id" = user_genome,
         id = `PATH/id10_gtdbtk.bac120.summary.tsv`,
         child = sapply(strsplit(id, "/", fixed = T), get_child)) %>%
  select(-id)

file_checkm_gc <- list.files(path = "PATH/concoct/bins_taxa/",
                             full.names = TRUE)
df_checkm_gc <- read_tsv(file_checkm_gc,
                         id = file_checkm_gc)
df_checkm_gc <- df_checkm_gc %>%
  mutate(id = `PATH/concoct/bins_taxa/id10_checkmTax.txt`,
         child = sapply(strsplit(id, "/", fixed = T), get_child)) %>%
  select(`Bin Id`, child, GC)

file_checkm_report <- list.files(path = "PATH/concoct/bins_report/",
                                 full.names = TRUE)
df_checkm_report <- read_tsv(file_checkm_report,
                             id = file_checkm_report)
df_checkm_report_filter <- df_checkm_report %>%
  mutate(id = `PATH/concoct/bins_report/id10_checkm_report.txt`,
         child = sapply(strsplit(id, "/", fixed = T), get_child)) %>%
  select(-id) %>%
  filter(Completeness >= 50 &
           (Contamination >= 10 |
              "Strain heterogeneity" >= 50))
df_checkm_report_filter <- left_join(df_checkm_report_filter, df_checkm_gc)
df_checkm_report_filter <- left_join(df_checkm_report_filter, df_gtdb)

df_bacteroides_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Bacteroides", classification, ignore.case = T) |
           grepl("g__Phocaeicola", classification, ignore.case = T)) %>%
  mutate(genus = "Bacteroides") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_bifido_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Bifidobacterium", classification, ignore.case = T)) %>%
  mutate(genus = "Bifidobacterium") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_clostridium_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Clostridium", classification, ignore.case = T)) %>%
  mutate(genus = "Clostridium") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_collinsella_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Collinsella", classification, ignore.case = T)) %>%
  mutate(genus = "Collinsella") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_enterococcus_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Enterococcus", classification, ignore.case = T)) %>%
  mutate(genus = "Enterococcus") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_escherichia_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Escherichia", classification, ignore.case = T)) %>%
  mutate(genus = "Escherichia") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_klebsiella_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Klebsiella", classification, ignore.case = T)) %>%
  mutate(genus = "Klebsiella") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_staphylococcus_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__staphylococcus", classification, ignore.case = T)) %>%
  mutate(genus = "Staphylococcus") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_streptococcus_MAG <- df_checkm_report_filter %>% 
  filter(grepl("g__Streptococcus", classification, ignore.case = T)) %>%
  mutate(genus = "Streptococcus") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_Bacteroidales_MAG <- df_checkm_report_filter %>% 
  filter(grepl("o__Bacteroidales", classification, ignore.case = T)) %>%
  mutate(genus = "pBHR_Bacteroidales") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_enterobacterales_MAG <- df_checkm_report_filter %>% 
  filter(grepl("o__Enterobacterales", classification, ignore.case = T)) %>%
  mutate(genus = "pBHR_Enterobacterales") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

df_lactobacillales_MAG <- df_checkm_report_filter %>% 
  filter(grepl("o__Lactobacillales", classification, ignore.case = T)) %>%
  mutate(genus = "pBHR_Lactobacillales") %>%
  select(`Bin Id`,child, GC, Completeness, Contamination, 
         `Strain heterogeneity`, classification)

################################################################################
## PLOTS
grDevices::windows(10,10)
par(mfrow=c(1,1))
par(oma=c(0,0,0,0))
par(mar=c(4,16,2,1))

boxplot(df_bacteroides$length,
        df_bifido$length,
        df_clostridium$length,
        df_collinsella$length,
        df_enterococcus$length,
        df_escherichia$length,
        df_klebsiella$length,
        df_staphylococcus$length,
        df_streptococcus$length,
        df_Bacteroidales$length,
        df_enterobacterales$length,
        df_lactobacillales$length,
        df_pbhrother$length,
        names = c("Bacteroides",
                  "Bifidobacterium",
                  "Clostridium",
                  "Collinsella",
                  "Enterococcus",
                  "Escherichia",
                  "Klebsiella",
                  "Staphylococcus",
                  "Streptococcus",
                  "pBHR Bacteroidales",
                  "pBHR Enterobacterales",
                  "pBHR Lactobacillales",
                  "pBHR other"),
        las = 2,
        ylab = "Length (kbp)")

boxplot(df_lactobacillales$GC,
        df_lactobacillales_MAG$GC,
        df_enterobacterales$GC,
        df_enterobacterales_MAG$GC,
        df_Bacteroidales$GC,
        df_Bacteroidales_MAG$GC,
        df_streptococcus$GC,
        df_streptococcus_MAG$GC,
        df_staphylococcus$GC,
        df_staphylococcus_MAG$GC,
        df_klebsiella$GC,
        df_klebsiella_MAG$GC,
        df_escherichia$GC,
        df_escherichia_MAG$GC,
        df_enterococcus$GC,
        df_enterococcus_MAG$GC,
        df_collinsella$GC,
        df_collinsella_MAG$GC,
        df_clostridium$GC,
        df_clostridium_MAG$GC,
        df_bifido$GC,
        df_bifido_MAG$GC,
        df_bacteroides$GC,
        df_bacteroides_MAG$GC,
        horizontal=TRUE,
        yaxt="n",
        las = 1,
        col = palette[5:4],
        xlab = "GC content",
        cex.axis = 1.5,
        cex.lab = 1.5)
axis(2, 
     at=c(1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5,17.5,19.5,21.5,23.5),
     labels = c("pBHR Lactobacillales",
                "pBHR Enterobacterales",
                "pBHR Bacteroidales",
                "Streptococcus",
                "Staphylococcus",
                "Klebsiella",
                "Escherichia",
                "Enterococcus",
                "Collinsella",
                "Clostridium",
                "Bifidobacterium",
                "Bacteroides"),
     las = 2,
     cex.axis = 1.5)
abline(h=c(2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5), lwd=1, col=1)
legend("topright",
       legend = c("MAG", "Plasmid"),
       fill = palette[4:5],
       cex = 1.5)

bacteroides_test <- wilcox.test(df_bacteroides$GC,
                                df_bacteroides_MAG$GC)
bifidobacter_test <- wilcox.test(df_bifido$GC,
                                 df_bifido_MAG$GC)
clostridium_test <- wilcox.test(df_clostridium$GC,
                                df_clostridium_MAG$GC)
collinsella_test <- wilcox.test(df_collinsella$GC,
                                df_collinsella_MAG$GC)
entereococcus_test <- wilcox.test(df_enterococcus$GC,
                                  df_enterococcus_MAG$GC)
escherichia_test <- wilcox.test(df_escherichia$GC,
                                df_escherichia_MAG$GC)
klebsiella_test <- wilcox.test(df_klebsiella$GC,
                               df_klebsiella_MAG$GC)
staphylococcus_test <- wilcox.test(df_staphylococcus$GC,
                                   df_staphylococcus_MAG$GC)
streptococcus_test <- wilcox.test(df_streptococcus$GC,
                                  df_streptococcus_MAG$GC)
Bacteroidales_test <- wilcox.test(df_Bacteroidales$GC,
                                 df_Bacteroidales_MAG$GC)
enterobacterales_test <- wilcox.test(df_enterobacterales$GC,
                                     df_enterobacterales_MAG$GC)
lactobacillales_test <- wilcox.test(df_lactobacillales$GC,
                                    df_lactobacillales_MAG$GC)

p_values <- c(bacteroides_test$p.value, 
              bifidobacter_test$p.value,
              clostridium_test$p.value,
              collinsella_test$p.value,
              entereococcus_test$p.value,
              escherichia_test$p.value,
              klebsiella_test$p.value,
              staphylococcus_test$p.value,
              streptococcus_test$p.value,
              Bacteroidales_test$p.value,
              enterobacterales_test$p.value,
              lactobacillales_test$p.value)
p_adj <- p.adjust(p_values, method = "BH")


case_when(is.na(p_adj) ~ "",
          p_adj > 0.05 ~ "",
          p_adj > 0.01 ~ "*",
          p_adj > 0.001 ~ "**",
          p_adj > 0.0001 ~ "***",
          .default = "****")
