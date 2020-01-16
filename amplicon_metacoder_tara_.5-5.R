library(dplyr)
library(tidyr)
library(metacoder) # 

################################################################################
#  get otu table and standardise
source("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/R/projects/Global_ocean_genome_analysis/amplicon_stats_tara_source.R")
#  source("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/Scripts/amplicon_stats_mal_source.R")
########################################################################################################
###  shan is my calc, number_of and shannon is from obis
#  for TARA : see ind stats to get list info
summary(dat_otu_reads) 
##   use dat !!!!!!!!!!      need lineage   summary(dat)   names(dat_otu_rare)   
##### combine if same species  but unlike source keep lineage  names(read_low)  names(dat)
#  for pid>90
sdat<- dat_otu_rare$read_low_sp  # 247160
sdat<-sdat[!duplicated(sdat$Sample.ID),]  ## check deep_65_DNA is duplicated
# names(sdat)  

#####!!!!!     limit DATA !!!!!!!!!  !!!!!!!  unique(sdat$filter_size_cat)
#sdat<-sdat[sdat$Depth_region=="SRF",]
sdat<-sdat[sdat$filter_size_cat==".5-5",]
################################################


tr<-c(23:length(names(sdat)))
### get unique lineage and otu tables
### get unique lineages  !!!!!      names(dat2)
lin_col<-c(2,36,37,39:45)
lin_otu<-dat2  %>% 
  select_at(lin_col) %>% 
  filter(!duplicated(cid)) %>% 
  mutate_all(as.character)
lin_sp <- lin_otu %>% 
  filter(!duplicated(sp_num))
########################################################################################################
#####   prep data for taxa to parse taxonomy     head(lin)
### input two sheets!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, 
###        otu - samples in columns and otu in rows
###         samples- samples in each row
##### prep Otu  #######   names(dat_otu)  summary(dat_otu_sam)   names(dat_otu_sam$read_low)
tr<-c(23:length(names(sdat)))  ## rows of taxa
dr<-c(1:22)   ## rows of data
otu<-sdat[,tr]
otu1<-data.frame(t(otu))  ### transpose so samples are in columns, taxa in rows
names(otu1)<-sdat$Sample.ID

otu1 <- otu1  %>%  # join to add lineage data     names(otu1)
  mutate(sp_num=as.character(row.names(otu1))) %>% 
  left_join(lin_sp)
  # names(sam$Depth)    hist(sam$Present.Surface.Temperature.Range)
# hist(sam$Land_Dist)
######   further filter taxa, remove sk and k and if not metazoan     names(otu)
otu2<-otu1 %>% 
  filter(Kingdom=="Metazoa")

####  make sure all taxa and samples have info
tr<-c(1:72)  # names(otu2) specify columns of samples
dr<-c(73:length(names(otu2)))  # data columns
#

otu3<-otu2[which(rowSums(otu2[,tr])>0),]  # remove taxa (rows) with no reads
otu3_dr<-otu3[,dr]
otu3_tr<-otu3[,tr]
mes<-otu3[,which(colSums(otu3_tr)>0)]  # remove samples (columns) with no reads
# !!! double check !!!!
otu4<-cbind(mes,otu3_dr)  # add back in lineage     names(otu4)

#   unique(otu4$Phylum)
# get most abund
### match samples from sam
#
#####  get sam info from stat       summary(stat_raw)
stat_sam <-stat_raw$read_low_sp %>%   #  names(stat_sam)
  dplyr::select(Sample.ID,shan:Pop_in_1000km)
#
dr<-c(1:22)   ## rows of data for sam     names(sdat)  
nam_mes<-as.character(names(mes))

sam<- sdat %>% ##  names(sdat)   names(sam)
  select_at(dr) %>% 
  filter(Sample.ID %in% names(mes)) %>% 
  left_join(stat_sam) %>% 
  #left_join(ob_cat) %>% 
  mutate(sst_mean_cat=cut(Present.Surface.Temperature.Mean, breaks=c(-Inf,20,25, Inf), 
                           labels=c("<20",">20 and <25",">25"))) %>% #   catagory  hist(sam$Present.Surface.Temperature.Mean)
  mutate(sst_range_cat=cut(Present.Surface.Temperature.Range, breaks=c(-Inf,4,12, Inf), 
                         labels=c("<4ºC",">4 and <12ºC",">12ºC"))) %>% #   catagory
  mutate(pp_mean_cat=cut(Present.Surface.Primary.productivity.Mean, breaks=c(-Inf,0.0025,0.005, Inf), 
                          labels=c("<0.0025",">0.0025 and <0.005",">.005"))) %>% #   catagory  hist(sam$Present.Surface.Primary.productivity.Mean)
  mutate(pp_range_cat=cut(Present.Surface.Primary.productivity.Range, breaks=c(-Inf,.005,.01, Inf), 
                           labels=c("<0.005",">0.005 and <0.01",">0.01"))) %>% #   catagory  hist(sam$Present.Surface.Primary.productivity.Range)
  mutate(land_dist_cat=cut(land_dist, breaks=c(-Inf,4,12, Inf), 
                           labels=c("<500",">500 and <1000",">1000"))) %>%  #land_dist cat
  mutate(OHI_cat=cut(OHI_2013, breaks=c(-Inf,3,4, Inf), 
                           labels=c("<3",">3 and <4",">4")))  %>%  #land_dist cat  hist(sam$OHI_2013)
    mutate(Latitude_abs=abs(Latitude)) %>% 
  mutate(Latitude_cat=cut(Latitude_abs, breaks=c(-Inf,23.5,35, Inf), 
                         labels=c("Tropical","Subtropical","Temperate")))  #   catagory

#################################################################################################
#################################################################################################
####   begin parsing taxa

tax<-taxa::parse_tax_data(otu4, class_cols =c("Kingdom","Phylum","Class","Order",
                                              "Family","Genus","Species"))

#  To plot read depth, you first need to add up the number of reads per taxon.
#  The function `calc_taxon_abund` is good for this. 
tax$data$taxon_counts <- calc_taxon_abund(tax, data = "tax_data")
#    names(otu4)
tr<-c(2:73)  # rows of samples taxon col 1   ^^^^  otu4 +1    ^^^^^^^^
tax$data$taxon_counts$total <- rowSums(tax$data$taxon_counts[,tr]) # -1 = 
#  names(tax$data$taxon_counts)

#print(tax)
#get_data(tax)
#tax$taxon_names


#metazoa_plot <- tax %>%
#filter_taxa(name == "Chordata", subtaxa = TRUE)  %>% 

heat_tree(tax, 
          node_label = taxon_names, 
          node_size = n_obs, 
          node_color = total,
          node_color_trans = "log10",
          node_color_range = rev(diverging_palette()),
          node_size_range = c(0.005, 0.026),
          node_label_size_range = c(0.020, 0.030),
          #node_label_size_trans = "area",  
          node_label_max = 70,
          overlap_avoidance = 0.7, # >per less overlap
          node_size_axis_label = "Number of species",
          node_color_axis_label = "Number of reads",
          initial_layout = "large-graph", layout = "davidson-harel")
## initial_layout = "fruchterman-reingold", layout = "davidson-harel",
## initial_layout = "davidson-harel", layout = "reingold-tilford" # ok better labels
## large-graph , gem-BAD longtime, mds,reingold-tilford

## for focus on richness - both number of species
heat_tree(tax, 
          node_label = taxon_names, 
          node_size = n_obs, 
          node_color = n_obs,
          node_color_trans = "log10",
          node_size_trans = "log10",
          node_color_range = diverging_palette(),
          node_label_size_range = c(0.020, 0.026),
          #node_label_size_trans = "area",
          node_label_max = 70,
          overlap_avoidance = 0.8, # >per less overlap
          node_size_axis_label = "Number of species",
          node_color_axis_label = "Number of species",
          initial_layout = "large-graph", layout = "davidson-harel")

#     save.image("~/Dropbox/Documents/KAUST/eDNA/DMAP/R/Enviros/metacoder_tara.RData")  # save R environment

########################################################################################################

# Compare categories ------------------------------------------------------
#   set order of catagories
# in order for tree matrix to plot corrreclty the sample data (sam) needs to be arranged for each separate run,
# each one below needs to be run separately before then running next section

sam<- sam %>% 
  arrange(factor(sst_range_cat, levels=c("<4ºC",">4 and <12ºC",">12ºC"))) 

sam<- sam %>% 
  arrange(factor(Latitude_cat, levels=c("Tropical","Subtropical","Temperate")))

sam<- sam %>% 
  arrange(factor(land_dist_cat, levels=c("<500",">500 and <1000",">1000"))) %>% 
  mutate(land_dist_cat2 = forcats::fct_recode(land_dist_cat, 
             "<500km"="<500",">500 and <1000km"=">500 and <1000",">1000km"=">1000"))

sam<- sam %>% 
  arrange(factor(filter_size_cat, levels=c(".5-5","5-20","20-200","180-2000")))  %>% 
  mutate(filter_size_cat2 = forcats::fct_recode(filter_size_cat, 
                          "0.5-5µm"=".5-5","5-20µm"="5-20","20-200µm"="20-200","180-2000µm"="180-2000"))
  #  levels(sam$filter_size_cat)
sam<- sam %>% 
  arrange(factor(Depth_region, levels=c("SRF","DCM","MEM"))) %>% #  unique(sam$Depth_region)
mutate(Depth_region2 = forcats::fct_recode(Depth_region, 
                          "Surface"="SRF","Deep chlorophyll maximum"="DCM","Epipelagic mixed layer"="MEM"))

sam <- sam %>% 
  arrange(factor(pp_mean_cat, levels=c("<0.0025",">0.0025 and <0.005",">.005"))) #  

sam<- sam %>% 
  arrange(factor(pp_range_cat, levels=c("<0.005",">0.005 and <0.01",">0.01"))) # 


###########################
###   get differences     cruise_leg     Depth_zone_cat  names(sam) 
##  run using "n_obs"   better for richness??
## was   "taxon_counts"
tax$data$diff_table <- compare_groups(tax, data = "taxon_counts",
                                      cols = sam$Sample.ID,
                                      groups = sam$pp_range_cat)
# tax$data$diff_table$log_mean_diff

###     use for above   -----  filter_size_cat  Depth_region2 sst_range_cat   Latitude_cat  land_dist_cat2
#    ocean_basin
#                    hist(tax$data$diff_table$log_mean_diff) 
##  change p_vlaue to adjusted
tax$data$diff_table$wilcox_p_value <- p.adjust(tax$data$diff_table$wilcox_p_value,
                                               method = "fdr")
## get log of mean difference
tax$data$diff_table$log_mean_diff<-log(abs(tax$data$diff_table$mean_diff+1))
tax$data$diff_table$log_mean_diff[tax$data$diff_table$mean_diff<1]<-tax$data$diff_table$log_mean_diff[tax$data$diff_table$mean_diff<1]*(-1)
## set anything > 0.05 to 0 in log2_median_ratio
#tax$data$diff_table$log_mean_diff[tax$data$diff_table$wilcox_p_value >0.05]<-0
#  mes<-tax$data$diff_table[tax$data$diff_table$wilcox_p_value <= 0.05 & complete.cases(tax$data$diff_table$wilcox_p_value),]

#
### plot
mdcol<-c(-5,5)
heat_tree_matrix(tax,
                 seed = 100,
                 data = "diff_table",
                 key_size=0.65 , #0.5 means half the width/height of the graph
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log_mean_diff, #  mean_diff
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = mdcol,
                 edge_color_interval = mdcol,
                 node_size_axis_label = "Number of taxa",
                 node_color_axis_label = "Log mean difference",
                 node_label_size_range = c(0.025, 0.035),
                 node_label_max = 55,
                 overlap_avoidance = 0.65,
                 initial_layout = "large-graph", layout = "davidson-harel")


###########################################################################################
# comparison stat sheet ---------------------------------------------------

names(sam)
tr<-c(4,22,39:46)  # pick relevent columns with variables that want results from
sam2<-sam %>% 
    select_at(tr)
pdat<-names(sam2)
# names of columns in tr
shnam<-c("dep.","filt.","sstm.","sstr.","ppm.","ppr.","ld.","OHI.")

# start lapply   i<-4
j<-1
for (i in tr) {
# get diff
  mes<-as.data.frame(cbind(sam[,1],sam[,i]))
  names(mes) [2]<-"var"
  mes2 <- compare_groups(tax, dataset = "taxon_counts",
                                      cols = mes$Sample.ID,
                                      groups = mes$var)
##  change p_vlaue to adjusted
  mes2$p_value_adj <- p.adjust(mes2$wilcox_p_value,
                                               method = "fdr")
  mes2$variable<-shnam[j]
  
  if (i==4){ tab<- mes2 
  } else { tab<-rbind(tab,mes2)
  }
  
  j<-j+1
}

### get taxa
ttt<-data.frame(cbind(tax$taxon_ids(),tax$taxon_names()))
names(ttt)<-c("taxon_id","taxon_names")
tab<-tab %>% 
  left_join(ttt)
  
# export
#     data.table::fwrite(tab,"/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/final_info/tara_metacoder2.csv",row.names=F, sep=",")


