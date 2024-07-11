
library(dplyr)
library(tidyr)
library(vegan)
# library(edgeR)

############################################################################################################
geo<-read.table("/Users/nathangeraldi/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/Global_layers_April18.csv", sep=",", header = TRUE)
##  get only uniqiue lat and long    names(geo)
geo1<-geo %>% 
  distinct(Latitude, Longitude, .keep_all = TRUE) %>% 
  mutate(Depth= replace(Depth, which(is.na(Depth) & Depth_region=="SRF"), 5)) %>% # not good when join later
  mutate(land_dist= replace(land_dist, which(is.na(land_dist) & Station=="TARA_011"), 0.01)) %>% 
  mutate(OHI_2013= replace(OHI_2013, which(is.na(OHI_2013) & Station==81), 3.77)) %>% 
  dplyr::rename(Depth_geo=Depth)
########################################################################################################
######    tara vargas data
dat<-data.table::fread(file = '/Users/nathangeraldi/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/groups/amplicon_Tara_metazoans.csv', sep = ',', header = TRUE)
names(dat)   # unique(dat$Cruise.x)
Dat<-dat[,1:47] # added 2 columns, dont want mess ups
pat<-c("Annelida","Arthropoda","Chordata","Cnidaria","Ctenophora","Echinodermata","Gastrotricha","Mollusca",
       "Nematoda","Orthonectida","Platyhelminthes","Porifera","Rotifera","Sipuncula","Tardigrada","Urochordata","Chaetognatha",
       "Craniata","Nemertea","Bryozoa","Hemichordata","Brachiopoda","Entoprocta","Cephalochordata","Entoprocta") 
pat<-paste(pat, collapse = "|")
#######   get basics   ######################################################
mes2<-dat %>% 
    group_by(filter_size,filter_size_cat) %>% 
    summarise(mean=mean(reads),sum=sum(reads),n(), sum_rar=sum(reads_rarified))

######   tdat   !!! add phylum column !!!!!!
dat2<- dat %>%    ### head(tdat)    unique(tdat$phylum)  4,871,092  unique(dat2$Sample.ID)
  filter(Cruise.x== "TARA") %>%
  filter(filter_size_cat != "0.8-20", filter_size_cat != "<0.8") %>% 
  filter(DNAorWGA=="DNA")

#### test for pid   #       hist(tdat$pid)   names(dat2)
pid_high<-97  # strict percent ID
pid_low<-90   # less strict  percent ID
###   column names should be reads and reads_rarified  !!!!!!!!!  change above if not
##  raw reads
rm(dat) ####  remove big dataframe
dat<-list()
read_low<- dat2 %>%    ### unique(read_low$otu)
  filter(pid>pid_low) %>% 
  filter(reads>0)  %>%
  mutate(otu=cid)
read_high<- dat2 %>%    ### 16133 records
  filter(pid>pid_high)  %>% 
  filter(reads>0) %>%
  mutate(otu=cid)
dat<-list(read_low=read_low,read_high=read_high)

##### combine if same species    names(read_low)  names(dat)   summary(dat)
#  for pid>90
dat[["read_low_sp"]]<- read_low %>%   # 247160
  group_by_at(vars(Sample.ID,Station:filter_size,OTU.LEVEL.SHANNON.DIVERSITY_from.Vargas:Locality.4,filter_size_cat,sp_num)) %>% 
  summarize(reads=sum(reads),reads_rarified=sum(reads_rarified), reads_percent=sum(reads_percent)) %>% 
  mutate(otu=sp_num) %>% 
  filter(complete.cases(otu))
dat[["read_high_sp"]]<- read_high %>%   # 247160
  group_by_at(vars(Sample.ID,Station:filter_size,OTU.LEVEL.SHANNON.DIVERSITY_from.Vargas:Locality.4,filter_size_cat,sp_num)) %>% 
  summarize(reads=sum(reads),reads_rarified=sum(reads_rarified), reads_percent=sum(reads_percent)) %>% 
  mutate(otu=sp_num) %>% 
  filter(complete.cases(otu))
##########################################
##### spread   #######   names(dat_otu)  summary(dat_otu_reads)   names(dat_otu_reads$read_low)   names(dat$read_low_sp)
dat_otu_reads<-lapply(dat,function(x) {
  x %>% 
    ungroup() %>% 
    dplyr::select(Sample.ID,TStation:NPP,OTU.LEVEL.SHANNON.DIVERSITY_from.Vargas:Locality.4,filter_size_cat,otu,reads) %>% #,OTU.LEVEL.SHANNON.DIVERSITY_from.Vargas:Locality.4
    spread(otu, reads, fill=0)
    })
### for rarefied      names(dat_otu_rare) summary(dat_otu_rare)  names(dat_otu_rare$read_low_sp)
dat_otu_rare<-lapply(dat,function(x) {
  x %>% 
    ungroup() %>% 
    filter(reads_rarified>0) %>%
    dplyr::select(Sample.ID,TStation:NPP,OTU.LEVEL.SHANNON.DIVERSITY_from.Vargas:Locality.4,filter_size_cat,otu,reads_rarified) %>% 
    spread(otu, reads_rarified, fill=0) 
  })
#   #   names(dat_otu_reads$read_low) 
#lapply(mylist, function(x) mean(x))
####################################################################################################
### calculate  riches and div - simpsons
#####    for all otus  --  otu     lapply(dat_otu, names)
nn<-23 ## first column of species
##### raw sequences       summary(dat_otu_reads)  names(dat_otu_reads$read_low[9530:9536])
#                                                 names(dat_otu_reads$read_low_sp[525:534])  
stat_raw<-lapply(dat_otu_reads,function(x) {
  shan <- as.numeric(diversity(x[,nn:ncol(x)]),index="simpson") # less biased to different n-Gihrig 2011
  rich <- as.numeric(rowSums(x[,nn:ncol(x)] > 0))  # plot(sample_rich)
  reads <- as.numeric(rowSums(x[,nn:ncol(x)]))
  # clean and join to geo   names(geo)  names(tdat90_otu[,1:30])    names(tdat90_otu)
  x %>% 
    bind_cols(shan=shan,rich=rich,reads=reads) %>% 
    mutate(Latitude_abs=abs(Latitude)) %>% 
    dplyr::select(Sample.ID:filter_size_cat,shan:Latitude_abs,OTU.LEVEL.SHANNON.DIVERSITY_from.Vargas)  %>% 
    mutate(filter_size_cat= factor(filter_size_cat, levels=c(".5-5","5-20","20-200","180-2000"))) %>% 
    left_join(geo1[,c(6:19)])  %>% 
    mutate(Depth= replace(Depth, which(is.na(Depth) & Depth_region=="SRF"), 5)) 
})   ###  summary(stat_raw)    names(stat_raw$read_low)  mes<-stat_raw$read_low_sp
#####rarefied
stat_rare<-lapply(dat_otu_rare,function(x) {
  shan <- as.numeric(diversity(x[,nn:ncol(x)]),index="simpson") # less biased to different n-Gihrig 2011
  rich <- as.numeric(rowSums(x[,nn:ncol(x)] > 0))  # plot(sample_rich)
  reads <- as.numeric(rowSums(x[,nn:ncol(x)]))
  # clean and join to geo   names(geo)  names(tdat90_otu[,1:30])    names(tdat90_otu)
  x %>% 
    bind_cols(shan=shan,rich=rich,reads=reads) %>% 
    mutate(Latitude_abs=abs(Latitude)) %>% 
    dplyr::select(Sample.ID:filter_size_cat,shan:Latitude_abs,OTU.LEVEL.SHANNON.DIVERSITY_from.Vargas)  %>% 
    mutate(filter_size_cat= factor(filter_size_cat, levels=c(".5-5","5-20","20-200","180-2000"))) %>% 
    left_join(geo1[,c(6:19)]) %>% 
    mutate(Depth= replace(Depth, which(is.na(Depth) & Depth_region=="SRF"), 5)) 
})
#   
########################################################################################################
########################################################################################################
