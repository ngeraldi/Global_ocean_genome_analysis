library(dplyr)
library(tidyr)
library(metacoder) # 

################################################################################
#  get otu table and standardise
source("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/R/projects/Global_ocean_genome_analysis/amplicon_stats_mal_source.R")
########################################################################################################

### get unique lineages  !!!!!
lin_col<-c(42,43,26,28:38)
lin_otu<-dat  %>% 
    select_at(lin_col) %>% 
    filter(!duplicated(otu_col)) %>% 
    mutate_all(as.character)
lin_sp <- lin_otu %>% 
    filter(!duplicated(sp_col))
#  for MAL
########################################################################################################
#####   use taxa to parse taxonomy     head(lin)
### input two sheets, 
###        otu - samples in columns and otu in rows
###         samples- samples in each row
##### spread based on samples  #######   names(dat_otu)  summary(dat_otu_sam)   names(dat_otu_sam$read_low)
#   use stat for mal.
#  1- "raw.otu.DNA_read_low"  2- "raw.species.DNA_read_low" 3- "rare.otu.DNA_read_low"   4- "rare.species.DNA_read_low"
sdat<-stat[["rare.species.DNA_read_low"]]  #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  names(sdat)      names(dat$read_low)
tr<-c(18:205)   ## rows of taxa
dr<-c(1:17,207:222)   ## rows of data for sam
sdat<-sdat[!duplicated(sdat$Sample.ID),]  ## check deep_65_DNA is duplicated

#####!!!!!     limit DATA !!!!!!!!!  !!!!!!!  unique(sdat$filter_size_cat)

sdat <- sdat   %>% 
  filter(!(grepl("deep",Sample.ID) | grepl("M",Sample.ID)))

otu<-sdat[,tr]
otu1<-data.frame(t(otu))  ### transpose so samples are in columns, taxa in rows
names(otu1)<-sdat$Sample.ID

  otu1 <- otu1  %>%  # join to add lineage data     names(otu1)
    mutate(sp_col=as.character(row.names(otu1))) %>% 
    left_join(lin_sp)
    #mutate(sp_num=as.character(sp_num))

#
######   further filter taxa, remove sk and k and if not metazoan     names(otu1)
tr<-c(1:81)  # names(otu1)
dr<-c(82:length(names(otu1)))
otu2<-otu1 %>% 
  filter(Kingdom == "Metazoa") 

otu3<-otu2[which(rowSums(otu2[,tr])>0),]  # remove taxa (rows) with no reads
otu3_dr<-otu3[,dr]
otu3_tr<-otu3[,tr]
mes<-otu3[,which(colSums(otu3_tr)>0)]  # remove samples (columns) with no reads
# !!! double check !!!!
otu4<-cbind(mes,otu3_dr)  # add back in lineage     names(otu4)
## make column of ocean basin
### make ocean basins
geo<-read.table("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/Global_layers_oct18.csv", sep=",",header=T)
geo1<-geo[,c(1,22)]# get lohg and ocean   names(geo)

### match samples from sam
dr<-c(1:17,205:222)   ## rows of data for sam    names(sdat)
nam_mes<-as.character(names(mes))
##   for sampls   names(sdat)        sdat$Sample.ID
sam<- sdat %>% 
  select_at(dr) %>% 
  filter(Sample.ID %in% names(mes)) %>% 
  left_join(geo1)  %>%  
  mutate(cruise_leg="profile") %>% # make column with different cruise legs
  mutate(cruise_leg=replace(cruise_leg,grepl("deep",Sample.ID),"deep")) %>% 
  mutate(cruise_leg=replace(cruise_leg,grepl("M",Sample.ID),"surface")) %>% 
  mutate(pelagic_zone=replace(pelagic_zone,grepl("M",Sample.ID),"E")) %>%   # fix pelagic zone catagories (E-<200, M-<1000 B)
  mutate(pelagic_zone=replace(pelagic_zone,grepl("deep",Sample.ID),"B")) %>%
  mutate(pp_mean_cat=cut(Present.Surface.Primary.productivity.Mean, breaks=c(-Inf,.002,.004, Inf), 
                             labels=c("<0.002",">0.002 and <0.004",">.004"))) %>% #  pp catagory
  # hist(sdat$Present.Surface.Primary.productivity.Mean)
  mutate(land_dist_cat=cut(land_dist, breaks=c(-Inf,5,10, Inf), 
                         labels=c("<500",">500 and <1000",">1000"))) %>%   # land_dist cat
  
  mutate(sst_mean_cat=cut(Present.Surface.Temperature.Mean, breaks=c(-Inf,20,25, Inf), 
                          labels=c("<20",">20 and <25",">25"))) %>% # hist(sdat$Present.Surface.Temperature.Mean)
  
  mutate(Depth_zone_cat="Epipelagic")  %>% 
  mutate(Depth_zone_cat=replace(Depth_zone_cat,pelagic_zone=="B","Bathypelagic")) %>%
  mutate(Depth_zone_cat=replace(Depth_zone_cat,pelagic_zone=="M","Mesopelagic")) 
#  3 dpeth cat, furface ppmean dist_land    hist(sdat$Depth)

####     a few checks
unique(otu2$Class)
mam<-otu1 %>% 
  filter(Class == "Aves")

otu5<-otu4
otu5$median <- apply(otu5[,1:221], 1, median)
mes<-otu5 %>% # names(otu4)
  mutate(sumread=rowSums(.[1:221])) %>% # rowSums(.[1:5])
  mutate(readmean=rowMeans(.[1:221])) %>% 
  select(sp_col:readmean) %>% 
  arrange(desc(sumread))  # median


########################################################################################################

# make tree data ----------------------------------------------------------
########################################################################################################

### all metazoa removed superkingdom
tax<-taxa::parse_tax_data(otu4, class_cols =c("Kingdom","Phylum","Class","Order",
                                              "Family","Genus","Species"))

#  tax2$data
#  sum per taxon?????
tax$data$tax_abund_per_sam <- calc_taxon_abund(tax, "tax_data",
                                       cols = sam$Sample.ID)
##  ra
# number of samples tha have reads for each taxon:  for depth catagories or DNA  names(sam)
  tax$data$tax_occ_depth <- calc_n_samples(tax, "tax_abund_per_sam", groups = sam$pelagic_zone)

#################################
names(otu4)
tr<-c(2:76)  # rows of samples taxon col 1
#  To plot read depth, you first need to add up the number of reads per taxon.
#  The function `calc_taxon_abund` is good for this. 
tax$data$taxon_counts <- calc_taxon_abund(tax, data = "tax_data")
tax$data$taxon_counts$total <- rowSums(tax$data$taxon_counts[,tr]) # 

#  print(tax)
#  get_data(tax)
#  tax$all_names()    tax$taxon_indexes  tax$classifications 

# plot --------------------------------------------------------------------

heat_tree(tax, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs, 
          node_size_axis_label = "OTU count",
          initial_layout = "large-graph", layout = "davidson-harel"
          )
# n_obs- number of otus
# total-number of reads
# n_sample - number of samples     edge_size = n_samples
# also use edge_size

#metazoa_plot <- tax %>%
#filter_taxa(name == "Chordata", subtaxa = TRUE)  %>% 
heat_tree(tax, 
          node_label = taxon_names, 
          node_size = n_obs,  # number of samples - right
          node_size_axis_label = "Number of samples",
          #node_size = "log10",
          node_size_range = c(0.01,0.05),
          node_color = total,  #total is number of reads-left of legend
          node_color_axis_label = "Number of reads",
          node_color_trans = "log10",
          node_label_size_range = c(0.02, 0.026),
          #node_label_size_trans = "area",
          node_label_max = 65,
          overlap_avoidance = 0.8,
          initial_layout = "large-graph", layout = "davidson-harel")
## initial_layout = "fruchterman-reingold", layout = "davidson-harel",
## initial_layout = "davidson-harel", layout = "reingold-tilford" # ok better labels
## large-graph , gem-BAD longtime, mds,reingold-tilford



# map comparisons ---------------------------------------------------------

sam<- sam %>% 
  arrange(factor(pp_mean_cat, levels=c("<0.002",">0.002 and <0.004",">.004"))) 


###########################    names(sam)

###   get differences     cruise_leg     Depth_zone_cat    ocean
tax$data$diff_table <- compare_groups(tax, data = "taxon_counts",
                                            cols = sam$Sample.ID,
                                            groups = sam$pp_mean_cat)
#   print(tax$data$diff_table_depth)  # keep only sig p_value (wilcox_p_value)

##  change p_vlaue to adjusted
tax$data$diff_table$wilcox_p_value <- p.adjust(tax$data$diff_table$wilcox_p_value,
                                               method = "fdr")
## get log of mean difference
tax$data$diff_table$log_mean_diff<-log(abs(tax$data$diff_table$mean_diff+1))
tax$data$diff_table$log_mean_diff[tax$data$diff_table$mean_diff<1]<-tax$data$diff_table$log_mean_diff[tax$data$diff_table$mean_diff<1]*(-1)

## set anything > 0.05 to 0 in log2_median_ratio
tax$data$diff_table$log2_median_ratio[tax$data$diff_table$wilcox_p_value >0.05]<-0
### plot
mdcol<-c(-3,3)
heat_tree_matrix(tax,
                 data = "diff_table",
                 key_size=0.6 , #0.5 means half the width/height of the graph
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log_mean_diff,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = mdcol,
                 edge_color_interval = mdcol,
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log mean difference",
                 node_label_size_range = c(0.030, 0.035),
                 node_size_range = c(0.010, 0.035),
                 node_label_max = 55,
                 overlap_avoidance = 0.75,
                 initial_layout = "large-graph", layout = "davidson-harel")
## grean or pos, more in treat 1 top ??



## only significant ?
#  hist(tax$data$diff_table$wilcox_p_value) 


#     save.image("~/Dropbox/Documents/KAUST/eDNA/DMAP/R/Enviros/metacoder_tara.RData")

########################################################################################################