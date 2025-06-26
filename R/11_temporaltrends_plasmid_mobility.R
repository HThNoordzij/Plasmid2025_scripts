# Clean up the memory before we start. 
rm(list=ls(all=TRUE)) 

# Load libraries
library("tidyverse")
library(mgcv)

library(RColorBrewer)
# display.brewer.all()
palette<-brewer.pal(n=8,name='Set1')
palette_children<-brewer.pal(n=12,name='Paired')

# Set working directory
setwd("~/PATH/")

# Load macSyFinder summary file
macSyFinder_file <- "plasmids_tax_mobility_clusters.csv"
df_macSyFinder <- read_csv(macSyFinder_file)

children <- df_macSyFinder %>% 
  select(child) %>% 
  unique() %>%
  unlist() 

df_conj <- df_macSyFinder %>%
  filter(mobility == "CONJ") %>%
  select(cluster, child) %>%
  unique()
df_mob <- df_macSyFinder %>%
  filter(mobility == "MOB") %>%
  select(cluster, child) %>%
  unique()
df_mobless <- df_macSyFinder %>%
  filter(mobility == "Mobless") %>%
  select(cluster, child) %>%
  unique()
df_all <- df_macSyFinder %>%
  select(cluster, child) %>%
  unique()

## Load abundance files
df_RPKM <- data.frame()
RPKM_file <- "abundance_summary_clusters.csv"
df_RPKM <- read_csv(RPKM_file, col_names = T)
df_RPKM <- df_RPKM %>%
  mutate(child = toupper(child)) %>%
  select(cluster, day, RPKM, child)

## separate per mobility type (and total)
df_RPKM_conj <- left_join(df_conj, df_RPKM)
df_RPKM_mob <- left_join(df_mob, df_RPKM)
df_RPKM_mobless <- left_join(df_mobless, df_RPKM)
df_RPKM_all <- left_join(df_all, df_RPKM)

## remove NAs
df_RPKM_conj <- df_RPKM_conj %>%
  filter(!is.na(RPKM))
df_RPKM_mob <- df_RPKM_mob %>%
  filter(!is.na(RPKM))
df_RPKM_mobless <- df_RPKM_mobless %>%
  filter(!is.na(RPKM))
df_RPKM_all <- df_RPKM_all %>%
  filter(!is.na(RPKM))

## sum total abundances per child 
df_RPKM_conj_sum <- data.frame()
df_RPKM_mob_sum <- data.frame()
df_RPKM_mobless_sum <- data.frame()
df_RPKM_all_sum <- data.frame()

for(i in 1:length(children)) {
  tmp_child <- children[i]
  
  ## pivot, have days as columns
  df_RPKM_conj_pivot <- df_RPKM_conj %>% 
    filter(child == tmp_child) %>%
    select(cluster,day,RPKM) %>%
    pivot_wider(names_from = day,
                values_from = RPKM,
                names_sort = TRUE,
                values_fill = 0)
  df_RPKM_mob_pivot <- df_RPKM_mob %>% 
    filter(child == tmp_child) %>%
    select(cluster,day,RPKM) %>%
    pivot_wider(names_from = day,
                values_from = RPKM,
                names_sort = TRUE,
                values_fill = 0)
  df_RPKM_mobless_pivot <- df_RPKM_mobless %>% 
    filter(child == tmp_child) %>%
    select(cluster,day,RPKM) %>%
    pivot_wider(names_from = day,
                values_from = RPKM,
                names_sort = TRUE,
                values_fill = 0)
  df_RPKM_all_pivot <- df_RPKM_all %>% 
    filter(child == tmp_child) %>%
    select(cluster,day,RPKM) %>%
    pivot_wider(names_from = day,
                values_from = RPKM,
                names_sort = TRUE,
                values_fill = 0)
  
  df_RPKM_conj_pivot <- as.data.frame(df_RPKM_conj_pivot)
  df_RPKM_mob_pivot <- as.data.frame(df_RPKM_mob_pivot)
  df_RPKM_mobless_pivot <- as.data.frame(df_RPKM_mobless_pivot)
  df_RPKM_all_pivot <- as.data.frame(df_RPKM_all_pivot)
  
  if(nrow(df_RPKM_conj_pivot) > 0) {
    cluster_names <- df_RPKM_conj_pivot$cluster
    rownames(df_RPKM_conj_pivot) <- cluster_names
    df_RPKM_conj_pivot <- df_RPKM_conj_pivot %>%
      select(-cluster)
  }
  
  if(nrow(df_RPKM_mob_pivot) > 0) {
    cluster_names <- df_RPKM_mob_pivot$cluster
    rownames(df_RPKM_mob_pivot) <- cluster_names
    df_RPKM_mob_pivot <- df_RPKM_mob_pivot %>%
      select(-cluster)
  }
  
  if(nrow(df_RPKM_mobless_pivot) > 0) {
    cluster_names <- df_RPKM_mobless_pivot$cluster
    rownames(df_RPKM_mobless_pivot) <- cluster_names
    df_RPKM_mobless_pivot <- df_RPKM_mobless_pivot %>%
      select(-cluster)
  }
  
  if(nrow(df_RPKM_all_pivot) > 0) {
    cluster_names <- df_RPKM_all_pivot$cluster
    rownames(df_RPKM_all_pivot) <- cluster_names
    df_RPKM_all_pivot <- df_RPKM_all_pivot %>%
      select(-cluster)
  }
  
  tmp_days <- colnames(df_RPKM_conj_pivot)
  if(nrow(df_RPKM_conj_pivot) > 0) {
    for (j in 1:length(tmp_days)) {
      tmp_row <- c(as.integer(tail(unlist(strsplit(tmp_child, "ID")),n=1)), 
                   tmp_days[j],
                   sum(df_RPKM_conj_pivot[,j], na.rm = TRUE))
      
      df_RPKM_conj_sum <- rbind(df_RPKM_conj_sum, tmp_row)
    }
  }
  
  tmp_days <- colnames(df_RPKM_mob_pivot)
  if(nrow(df_RPKM_mob_pivot) > 0) {
    for (j in 1:length(tmp_days)) {
      tmp_row <- c(as.integer(tail(unlist(strsplit(tmp_child, "ID")),n=1)), 
                   tmp_days[j],
                   sum(df_RPKM_mob_pivot[,j], na.rm = TRUE))
      
      df_RPKM_mob_sum <- rbind(df_RPKM_mob_sum, tmp_row)
    }
  }
  
  tmp_days <- colnames(df_RPKM_mobless_pivot)
  if(nrow(df_RPKM_mobless_pivot) > 0) {
    for (j in 1:length(tmp_days)) {
      tmp_row <- c(as.integer(tail(unlist(strsplit(tmp_child, "ID")),n=1)), 
                   tmp_days[j],
                   sum(df_RPKM_mobless_pivot[,j], na.rm = TRUE))
      
      df_RPKM_mobless_sum <- rbind(df_RPKM_mobless_sum, tmp_row)
    }
  }
  
  
  tmp_days <- colnames(df_RPKM_all_pivot)
  if(nrow(df_RPKM_all_pivot) > 0) {
    for (j in 1:length(tmp_days)) {
      tmp_row <- c(as.integer(tail(unlist(strsplit(tmp_child, "ID")),n=1)), 
                   tmp_days[j],
                   sum(df_RPKM_all_pivot[,j], na.rm = TRUE))
      
      df_RPKM_all_sum <- rbind(df_RPKM_all_sum, tmp_row)
    }
  }
  
  
}
colnames(df_RPKM_conj_sum) <- c("child", "day", "sumdepth")
colnames(df_RPKM_mob_sum) <- c("child", "day", "sumdepth")
colnames(df_RPKM_mobless_sum) <- c("child", "day", "sumdepth")
colnames(df_RPKM_all_sum) <- c("child", "day", "sumdepth")

### remove NAs
df_RPKM_conj_sum <- df_RPKM_conj_sum %>%
  filter(!is.na(sumdepth))
df_RPKM_mob_sum <- df_RPKM_mob_sum %>%
  filter(!is.na(sumdepth))
df_RPKM_mobless_sum <- df_RPKM_mobless_sum %>%
  filter(!is.na(sumdepth))
df_RPKM_all_sum <- df_RPKM_all_sum %>%
  filter(!is.na(sumdepth))

### Add zero's for infants with >1 plasmid for days of no RPKM
## df_RPKM_all_sum
for (i in 1:12) {
  tmp_child <- as.character(i)
  
  
  tmp_all_days <- df_RPKM %>%
    filter(child == paste0("ID",tmp_child)) %>%
    select(day) %>%
    unique() %>%
    unlist()
  
  df_tmp <- df_RPKM_all_sum %>%
    filter(child == tmp_child)
  
  if(tmp_child %in% df_RPKM_all_sum$child){
    for (j in 1:length(tmp_all_days)) {
      tmp_day <- as.integer(tmp_all_days[j])
      if(!tmp_day %in% df_tmp$day) {
        tmp_row <- c(tmp_child, tmp_day, 0)
        df_RPKM_all_sum <- rbind(df_RPKM_all_sum, tmp_row)
      }
    }
  }
}

## df_RPKM_conj_sum
for (i in 1:12) {
  tmp_child <- as.character(i)
  
  
  tmp_all_days <- df_RPKM %>%
    filter(child == paste0("ID",tmp_child)) %>%
    select(day) %>%
    unique() %>%
    unlist()
  
  df_tmp <- df_RPKM_conj_sum %>%
    filter(child == tmp_child)
  
  if(tmp_child %in% df_RPKM_conj_sum$child) {
    for (j in 1:length(tmp_all_days)) {
      tmp_day <- as.integer(tmp_all_days[j])
      if(!tmp_day %in% df_tmp$day) {
        tmp_row <- c(tmp_child, tmp_day, as.double(0))
        df_RPKM_conj_sum <- rbind(df_RPKM_conj_sum, tmp_row)
      }
    }
  }
}

## df_RPKM_mob_sum
for (i in 1:12) {
  tmp_child <- as.character(i)
  
  
  tmp_all_days <- df_RPKM %>%
    filter(child == paste0("ID",tmp_child)) %>%
    select(day) %>%
    unique() %>%
    unlist()
  
  df_tmp <- df_RPKM_mob_sum %>%
    filter(child == tmp_child)
  
  if(tmp_child %in% df_RPKM_mob_sum$child) {
    for (j in 1:length(tmp_all_days)) {
      tmp_day <- as.integer(tmp_all_days[j])
      if(!tmp_day %in% df_tmp$day) {
        tmp_row <- c(tmp_child, tmp_day, as.double(0))
        df_RPKM_mob_sum <- rbind(df_RPKM_mob_sum, tmp_row)
      }
    }
  }
}

## df_RPKM_mobless_sum
for (i in 1:12) {
  tmp_child <- as.character(i)
  
  
  tmp_all_days <- df_RPKM %>%
    filter(child == paste0("ID",tmp_child)) %>%
    select(day) %>%
    unique() %>%
    unlist()
  
  df_tmp <- df_RPKM_mobless_sum %>%
    filter(child == tmp_child)
  
  if(tmp_child %in% df_RPKM_mobless_sum$child) {
    for (j in 1:length(tmp_all_days)) {
      tmp_day <- as.integer(tmp_all_days[j])
      if(!tmp_day %in% df_tmp$day) {
        tmp_row <- c(tmp_child, tmp_day, as.double(0))
        df_RPKM_mobless_sum <- rbind(df_RPKM_mobless_sum, tmp_row)
      }
    }
  }
}

#####
df_RPKM_conj_sum <- df_RPKM_conj_sum %>%
  mutate(rowname = paste0(child, "_d", day))
df_RPKM_mob_sum <- df_RPKM_mob_sum %>%
  mutate(rowname = paste0(child, "_d", day))
df_RPKM_mobless_sum <- df_RPKM_mobless_sum %>%
  mutate(rowname = paste0(child, "_d", day))
df_RPKM_all_sum <- df_RPKM_all_sum %>%
  mutate(rowname = paste0(child, "_d", day))

conj_mapped_plasmid <- as.double(df_RPKM_conj_sum$sumdepth)
names(conj_mapped_plasmid) <- df_RPKM_conj_sum$rowname
conj_profiles_filt_norm_plasmid <- as.data.frame(df_RPKM_conj_sum)
conj_days_plasmid<-unlist(lapply(strsplit(names(conj_mapped_plasmid),split='_d'),function(x)x[2]))

mob_mapped_plasmid <- as.double(df_RPKM_mob_sum$sumdepth)
names(mob_mapped_plasmid) <- df_RPKM_mob_sum$rowname
mob_profiles_filt_norm_plasmid <- as.data.frame(df_RPKM_mob_sum)
mob_days_plasmid<-unlist(lapply(strsplit(names(mob_mapped_plasmid),split='_d'),function(x)x[2]))

mobless_mapped_plasmid <- as.double(df_RPKM_mobless_sum$sumdepth)
names(mobless_mapped_plasmid) <- df_RPKM_mobless_sum$rowname
mobless_profiles_filt_norm_plasmid <- as.data.frame(df_RPKM_mobless_sum)
mobless_days_plasmid<-unlist(lapply(strsplit(names(mobless_mapped_plasmid),split='_d'),function(x)x[2]))

all_mapped_plasmid <- as.double(df_RPKM_all_sum$sumdepth)
names(all_mapped_plasmid) <- df_RPKM_all_sum$rowname
all_profiles_filt_norm_plasmid <- as.data.frame(df_RPKM_all_sum)
all_days_plasmid<-unlist(lapply(strsplit(names(all_mapped_plasmid),split='_d'),function(x)x[2]))


mob_id_plasmid<-unlist(lapply(strsplit(names(mob_mapped_plasmid),split='_'),function(x)x[1]))
colors_mob_plasmid<-vector()
colors_mob_plasmid[mob_id_plasmid=='1']<-1
colors_mob_plasmid[mob_id_plasmid=='2']<-2
colors_mob_plasmid[mob_id_plasmid=='3']<-3
colors_mob_plasmid[mob_id_plasmid=='4']<-4
colors_mob_plasmid[mob_id_plasmid=='5']<-5
colors_mob_plasmid[mob_id_plasmid=='6']<-6
colors_mob_plasmid[mob_id_plasmid=='7']<-7
colors_mob_plasmid[mob_id_plasmid=='8']<-8
colors_mob_plasmid[mob_id_plasmid=='9']<-'brown'
colors_mob_plasmid[mob_id_plasmid=='10']<-'purple'
colors_mob_plasmid[mob_id_plasmid=='11']<-'darkblue'
colors_mob_plasmid[mob_id_plasmid=='12']<-'orange'

mobless_id_plasmid<-unlist(lapply(strsplit(names(mobless_mapped_plasmid),split='_'),function(x)x[1]))
colors_mobless_plasmid<-vector()
colors_mobless_plasmid[mobless_id_plasmid=='1']<-1
colors_mobless_plasmid[mobless_id_plasmid=='2']<-2
colors_mobless_plasmid[mobless_id_plasmid=='3']<-3
colors_mobless_plasmid[mobless_id_plasmid=='4']<-4
colors_mobless_plasmid[mobless_id_plasmid=='5']<-5
colors_mobless_plasmid[mobless_id_plasmid=='6']<-6
colors_mobless_plasmid[mobless_id_plasmid=='7']<-7
colors_mobless_plasmid[mobless_id_plasmid=='8']<-8
colors_mobless_plasmid[mobless_id_plasmid=='9']<-'brown'
colors_mobless_plasmid[mobless_id_plasmid=='10']<-'purple'
colors_mobless_plasmid[mobless_id_plasmid=='11']<-'darkblue'
colors_mobless_plasmid[mobless_id_plasmid=='12']<-'orange'

all_id_plasmid<-unlist(lapply(strsplit(names(all_mapped_plasmid),split='_'),function(x)x[1]))
colors_all_plasmid<-vector()
colors_all_plasmid[all_id_plasmid=='1']<-1
colors_all_plasmid[all_id_plasmid=='2']<-2
colors_all_plasmid[all_id_plasmid=='3']<-3
colors_all_plasmid[all_id_plasmid=='4']<-4
colors_all_plasmid[all_id_plasmid=='5']<-5
colors_all_plasmid[all_id_plasmid=='6']<-6
colors_all_plasmid[all_id_plasmid=='7']<-7
colors_all_plasmid[all_id_plasmid=='8']<-8
colors_all_plasmid[all_id_plasmid=='9']<-'brown'
colors_all_plasmid[all_id_plasmid=='10']<-'purple'
colors_all_plasmid[all_id_plasmid=='11']<-'darkblue'
colors_all_plasmid[all_id_plasmid=='12']<-'orange'

conj_id_plasmid<-unlist(lapply(strsplit(names(conj_mapped_plasmid),split='_'),function(x)x[1]))
colors_conj_plasmid<-vector()
colors_conj_plasmid[conj_id_plasmid=='1']<-1
colors_conj_plasmid[conj_id_plasmid=='2']<-2
colors_conj_plasmid[conj_id_plasmid=='3']<-3
colors_conj_plasmid[conj_id_plasmid=='4']<-4
colors_conj_plasmid[conj_id_plasmid=='5']<-5
colors_conj_plasmid[conj_id_plasmid=='6']<-6
colors_conj_plasmid[conj_id_plasmid=='7']<-7
colors_conj_plasmid[conj_id_plasmid=='8']<-8
colors_conj_plasmid[conj_id_plasmid=='9']<-'brown'
colors_conj_plasmid[conj_id_plasmid=='10']<-'purple'
colors_conj_plasmid[conj_id_plasmid=='11']<-'darkblue'
colors_conj_plasmid[conj_id_plasmid=='12']<-'orange'
##
conj_days_plasmid<-as.numeric(conj_days_plasmid)
mod_conj_mapped_plasmid <- gam(conj_mapped_plasmid~s(conj_days_plasmid,k=5))
summary(mod_conj_mapped_plasmid)

mob_days_plasmid<-as.numeric(mob_days_plasmid)
mod_mob_mapped_plasmid <- gam(mob_mapped_plasmid~s(mob_days_plasmid,k=5))
summary(mod_mob_mapped_plasmid)

mobless_days_plasmid<-as.numeric(mobless_days_plasmid)
mod_mobless_mapped_plasmid <- gam(mobless_mapped_plasmid~s(mobless_days_plasmid,k=5))
summary(mod_mobless_mapped_plasmid)

all_days_plasmid<-as.numeric(all_days_plasmid)
mod_all_mapped_plasmid <- gam(all_mapped_plasmid~s(all_days_plasmid,k=5))
summary(mod_all_mapped_plasmid)


############ plotting
grDevices::windows(10,10)
par(mfrow=c(2,2))
par(oma=c(2,2,0,0))
par(mar=c(2,2,2,3))

plot(mod_all_mapped_plasmid,
     residuals=T,
     shift=mean(all_mapped_plasmid),
     # pch=16,
     col= 'black',
     rug=F,
     lwd=2,
     shade=T,
     shade.col='grey90',
     xlab='',
     ylab='',
     cex=1,
     main='Plasmid (all)') 
points(all_mapped_plasmid~all_days_plasmid,
       pch=16,
       cex=1,
       col=colors_all_plasmid)
legend('topleft',
       legend=c('ID1','ID2','ID3','ID4','ID5','ID6','ID7','ID8','ID9','ID10','ID11','ID12'),
       fill=c(1:8, 'brown','purple','darkblue','orange'),cex=.7)
Rlabel_plasmid <- bquote(italic(R)^2 == .(format(round(summary(mod_all_mapped_plasmid)$r.sq, digits = 3))))
Plabel_plasmid <- bquote(italic(P) == .(format(round(summary(mod_all_mapped_plasmid)$s.table[,4], digits = 1))))
Nlabel_plasmid <- bquote(n == .(format(nrow(df_all))))
mylabel_plasmid <- c(Rlabel_plasmid, Plabel_plasmid, Nlabel_plasmid)
legend('top', legend = mylabel_plasmid, bty = 'n')


plot(mod_conj_mapped_plasmid,
     residuals=T,
     shift=mean(conj_mapped_plasmid),
     # pch=16,
     col='black',
     rug=F,
     lwd=2,
     shade=T,
     shade.col='grey90',
     xlab='',ylab='',
     cex=1,
     main='Plasmid (CONJ)') 
points(conj_mapped_plasmid~conj_days_plasmid,
       pch=16,
       cex=1,
       col=colors_conj_plasmid)
Rlabel_plasmid <- bquote(italic(R)^2 == .(format(round(summary(mod_conj_mapped_plasmid)$r.sq, digits = 3))))
Plabel_plasmid <- bquote(italic(P) == .(format(round(summary(mod_conj_mapped_plasmid)$s.table[,4], digits = 2))))
Nlabel_plasmid <- bquote(n == .(format(nrow(df_conj))))
mylabel_plasmid <- c(Rlabel_plasmid, Plabel_plasmid, Nlabel_plasmid)
legend('topright', legend = mylabel_plasmid, bty = 'n')

plot(mod_mob_mapped_plasmid,
     residuals=T,
     shift=mean(mob_mapped_plasmid),
     # pch=16,
     col='black',
     rug=F,
     lwd=2,
     shade=T,
     shade.col='grey90',
     xlab='',ylab='',
     cex=1,
     main='Plasmid (MOB)') 
points(mob_mapped_plasmid~mob_days_plasmid,
       pch=16,
       cex=1,
       col=colors_mob_plasmid)
Rlabel_plasmid <- bquote(italic(R)^2 == .(format(round(summary(mod_mob_mapped_plasmid)$r.sq, digits = 3))))
Plabel_plasmid <- bquote(italic(P) == .(format(round(summary(mod_mob_mapped_plasmid)$s.table[,4], digits = 2))))
Nlabel_plasmid <- bquote(n == .(format(nrow(df_mob))))
mylabel_plasmid <- c(Rlabel_plasmid, Plabel_plasmid, Nlabel_plasmid)
legend('top', legend = mylabel_plasmid, bty = 'n')

plot(mod_mobless_mapped_plasmid,
     residuals=T,
     shift=mean(mobless_mapped_plasmid),
     # pch=16,
     col='black',
     rug=F,
     lwd=2,
     shade=T,
     shade.col='grey90',
     xlab='',ylab='',
     cex=1,
     main='Plasmid (MOBless)') 
points(mobless_mapped_plasmid~mobless_days_plasmid,
       pch=16,
       cex=1,
       col=colors_mobless_plasmid)
Rlabel_plasmid <- bquote(italic(R)^2 == .(format(round(summary(mod_mobless_mapped_plasmid)$r.sq, digits = 3))))
Plabel_plasmid <- bquote(italic(P) == .(format(round(summary(mod_mobless_mapped_plasmid)$s.table[,4], digits = 2))))
Nlabel_plasmid <- bquote(n == .(format(nrow(df_mobless))))
mylabel_plasmid <- c(Rlabel_plasmid, Plabel_plasmid, Nlabel_plasmid)
legend('top', legend = mylabel_plasmid, bty = 'n')

mtext("RPKM",
      side = 2, 
      line = 0.8,
      outer=TRUE)
mtext("Days",
      side = 1, 
      line = 0.8,
      outer=TRUE)


###############################################################################
## Plots per child
grDevices::windows(10,10)
par(mfrow=c(3,4))
par(oma=c(2,2,0,0))
par(mar=c(2,2,2,1))

for (i in 1:12) {
  tmp_child <- i
  df_tmp_all <- df_RPKM_all_sum %>%
    filter(child == tmp_child) %>%
    mutate(day = as.integer(day)) %>%
    arrange(day)
  tmp_days <- df_RPKM %>%
    filter(child == paste0("ID",tmp_child)) %>%
    select(day) %>%
    unique() %>% unlist()
  
  df_tmp_CONJ <- df_RPKM_conj_sum %>%
    filter(child == tmp_child)
  df_tmp_MOB <- df_RPKM_mob_sum %>%
    filter(child == tmp_child)
  df_tmp_MOBless <- df_RPKM_mobless_sum %>%
    filter(child == tmp_child)
  
  ## add missing days
  for (j in 1:length(tmp_days)) {
    tmp_day <- tmp_days[[j]]
    if(nrow(df_tmp_CONJ %>% filter(day == tmp_day)) == 0){
      tmp_row <- c(tmp_child, tmp_day, 0, NA)
      df_tmp_CONJ <- rbind(df_tmp_CONJ, tmp_row)
      colnames(df_tmp_CONJ) <- colnames(df_RPKM_conj_sum)
    }
    if(nrow(df_tmp_MOB %>% filter(day == tmp_day)) == 0){
      tmp_row <- c(tmp_child, tmp_day, 0, NA)
      df_tmp_MOB <- rbind(df_tmp_MOB, tmp_row)
      colnames(df_tmp_MOB) <- colnames(df_RPKM_mob_sum)
    }
    if(nrow(df_tmp_MOBless %>% filter(day == tmp_day)) == 0){
      tmp_row <- c(tmp_child, tmp_day, 0, NA)
      df_tmp_MOBless <- rbind(df_tmp_MOBless, tmp_row)
      colnames(df_tmp_MOBless) <- colnames(df_RPKM_mobless_sum)
    }
  }
  
  df_tmp_CONJ <- df_tmp_CONJ %>% 
    mutate(day = as.integer(day)) %>%
    arrange(day)
  df_tmp_MOB <- df_tmp_MOB %>% 
    mutate(day = as.integer(day)) %>%
    arrange(day)
  df_tmp_MOBless <- df_tmp_MOBless %>% 
    mutate(day = as.integer(day)) %>%
    arrange(day)
  
  tmp_title <- paste0("ID", tmp_child)
  
  if(nrow(df_tmp_all) > 0) {
    plot(df_tmp_all$day,
         df_tmp_all$sumdepth,
         xlim = c(0,366),
         xlab = "",
         ylab = "",
         ylim = c(0,max(as.integer(df_tmp_all$sumdepth)) + 50),
         main = tmp_title,
         type = "l",
         lwd = 2)
  } else {
    plot(1, 
         type="n",
         xlab="",
         ylab="", 
         xlim=c(0, 366), 
         ylim=c(0, 200),
         main = tmp_title)
    abline(h = 0,
           lwd = 2)
  }
  
  if(sum(as.double(df_tmp_CONJ$sumdepth)) > 0) {
    lines(df_tmp_CONJ$day,
          df_tmp_CONJ$sumdepth,
          col = palette[1],
          lwd = 2)
  }
  if(sum(as.double(df_tmp_MOB$sumdepth)) > 0) {
    lines(df_tmp_MOB$day,
          df_tmp_MOB$sumdepth,
          col = palette[2],
          lwd = 2)
  }
  if(sum(as.double(df_tmp_MOBless$sumdepth)) > 0)  {
    lines(df_tmp_MOBless$day,
          df_tmp_MOBless$sumdepth,
          col = palette[3],
          lwd = 2)
  }
  
  if(tmp_child == 1) {
    legend('topright',
           legend = c("All", "CONJ", "MOB", "MOBless"),
           fill = c("black", palette[1:3]))
  }
}

mtext("Days",
      side = 1, 
      line = 0.8,
      outer=TRUE)
mtext("RPKM",
      side = 2, 
      line = 0.8,
      outer=TRUE, 
      las=0)

###############################################################################

### correlation tests
corr_all <- cor.test(all_days_plasmid, all_mapped_plasmid, method = 'spearman')
corr_conj <- cor.test(conj_days_plasmid, conj_mapped_plasmid, method = 'spearman')
corr_mob <- cor.test(mob_days_plasmid, mob_mapped_plasmid, method = 'spearman')
corr_mobless <- cor.test(mobless_days_plasmid, mobless_mapped_plasmid, method = 'spearman')

df_corr <- data.frame("Class" = c("All", "CONJ", "MOB", "MOBless"),
                      "Rho" = c(corr_all$estimate,
                                corr_conj$estimate,
                                corr_mob$estimate,
                                corr_mobless$estimate),
                      "Pinit" = c(corr_all$p.value,
                                  corr_conj$p.value,
                                  corr_mob$p.value,
                                  corr_mobless$p.value))

df_corr <- df_corr %>%
  mutate(P = p.adjust(Pinit, method="BH")) %>%
  select(Class,
         Rho,
         P)

## Save Pvalues/rho of spearman
file_spearman_out <- "spearman_classification.csv"
write.csv(df_corr,
          file = file_spearman_out, row.names = FALSE, quote = FALSE)
