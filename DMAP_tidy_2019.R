
setwd("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/work")
##      save.image("DMAP4.rdata")    ###   to save to above directory  -  getting info from web takes time, so save it once you have it
##       load("DMAP4.rdata")         #####   load from above directory

#  library(phyloseq)
library(taxize)  #   install.packages("taxize")
#library(ape)  #   install.packages("ape")
#library(gplots)  #   install.packages("gplots")
library(RColorBrewer)
#library(repmis)  #   install.packages("repmis")
library(tidyverse)  #   install.packages("vegan")

#########################################################################################
###### "biomass"  data    april 2019  ####################################################################
##  all data
tax1<-read.csv("/Users/geraldn/Dropbox/Global_databases/DMAP/DMAP_for_biomass/all_seq.csv", header = T)
n2<-read.csv("/Users/geraldn/Dropbox/Global_databases/DMAP/DMAP_for_biomass/sp_ids.csv", header = T)

########################################################################################################################################
###### TARA data     ####################################################################
#########################################################################################################
## import data
Ttax<- tax1 %>%   # names(tax1)    unique(tax1$cruise)
  filter(cruise=="tar")
Tsam<-read.csv("/Users/geraldn/Dropbox/Global_databases/DMAP/TARA/tara.stations.good1.csv", header = T)
Tenv<-read.csv("/Users/geraldn/Dropbox/Global_databases/DMAP/TARA/TARA_env_var.csv", header = T)

#  tidy environmental data
## to do split column to match sample data at lat lon depth
 #  add meta data neead sample ID       TARA_036_DCM_0.22.1.6_IO.d17
Env.feature<-levels(Tenv$Env.feature)
Depth_region<-c("DCM","DCM","MEM","MES","MES","SRF")
dep<-data.frame(cbind(Env.feature,Depth_region))      ### names(Tenv1)
Tenv1<- Tenv %>%
  separate(Station, c("Cruise","Station"), sep = "_", remove=T)  %>%
  left_join(dep,by="Env.feature") %>%
  mutate(env_ID = paste(Station, Depth_region, sep = '.')) %>%
  dplyr::select(Sample.ID,Cruise,Station,Depth_region,env_ID,Depth.water..m.:Tpot...C.,OXYGEN..µmol.kg...calculated.from.in.situ.senso....:Chl.a..mg.m..3...calculated.from.in.situ.senso....,NPP.C..mg.m..2.day...at.the.sampling.location.for.....)  %>%
  group_by(env_ID) %>%
  summarise_at(c(5:17), mean, na.rm = TRUE) %>%
  dplyr::rename(Depth.env=Depth.water..m.,Temp.=Tpot...C.,Oxygen=OXYGEN..µmol.kg...calculated.from.in.situ.senso....,Nitrate=X.NO3....µmol.l.,Nitrite=X.NO2....µmol.l.,Phosphate=PO4..µmol.l.,Silicate=Si.OH.4..µmol.l.,Chlorophyll=Chl.a..mg.m..3...calculated.from.in.situ.senso....,NPP=NPP.C..mg.m..2.day...at.the.sampling.location.for.....) %>%
  dplyr::select(env_ID,Depth.env,Sal:Oxygen,Nitrate, Nitrite,Phosphate,Silicate:NPP)
    # names(Tenv1)     names(Ttax1)    names(Ttax1[1000:1141])
# tidy sam then tax data
Tsam1<- Tsam  %>% 
     separate(Station, c("Cruise","Station"), sep = "_", remove=T) %>% 
    select(-cat)
Ttax1<- Ttax %>%   ## names(Ttax)  head(Ttax)
   select(-cruise,-tax_catagory ) %>%   # for 2019
   separate(sample_ID, c("Cruise","Station","Depth_region","filter_size","filt2"), sep="_", remove=F) %>% 
   unite(filter_size, c("filter_size","filt2"), sep = "-", remove=T)  %>%   # unique(Ttax1$Depth)
   mutate(env_ID = paste(Station, Depth_region, sep = '.')) %>%
   left_join(Tsam1, by=c("Cruise","Station")) %>% 
   dplyr::rename(ID = sample_ID)  %>%
   left_join(Tenv1, by="env_ID")

#  add unique sample IDs
q<-data.frame(unique(Ttax1$ID))    #  names(q)
q$Sample<- c(1:length(q$unique.Ttax1.ID.))
names(q)[1]<-paste("ID")
# join and simplify     names(Ttax1)   names(Ttax1[220:length(names(Ttax1))])
Ttax1<- Ttax1 %>%  
  left_join(q, by=c("ID")) %>% 
   select(ID,Sample,Station,Cruise,area,filter_size,gene_domain,pid,Date.Time:NPP,Hydra.vulgaris:Millepora.sp..EK.2011)# names(Ttax1)  names(Ttax1[999:1131])   

########################################################################################################################################
###### Malispina data     ####################################################################
#########################################################################################################

## mala overall

#  get station !!
Msam_amp<-data.table::fread("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/amplicons_Mal_sample_info.csv",header = TRUE, sep=",")
Msam_amp <- Msam_amp %>%     # names(Msam_amp)
  select(Station, Latitude, Longitude) %>% 
  filter(!duplicated(Station)) 
# Msam_amp[,c(2:3)] <- lapply(Msam_amp[,c(2:3)], round, 4)

Mtax<- tax1 %>%   # names(tax1)    unique(tax1$cruise)
  filter(cruise=="mal") %>% 
  dplyr::rename(sample_ID_long=sample_ID) %>% 
  mutate(sample_ID=gsub("\\..*$","",sample_ID_long)) %>% 
  mutate(sample_ID=as.character(sample_ID), sample_ID_long=as.character(sample_ID_long)) %>% 
  mutate(sample_ID=if_else(area=="prof", sample_ID_long, sample_ID ))
  # get correct sampl_ID    length(unique(Mtax$sample_ID))   length(unique(Mtax$sample_ID_long))
  mes<-Mtax[duplicated(Mtax$sample_ID),]
Menv<-read.csv("/Users/geraldn/Dropbox/Global_databases/DMAP//DMAP_meta_land.dist.csv", header = T)
Msam<-read.csv("/Users/geraldn/Dropbox/Global_databases/DMAP/Mal/Mal_reads_region.csv", header = T)
#  add meta data neead sample ID     names(Mtax)        names(msam)          names(Menv1)     names(Menv)
Menv1<- Menv %>%
  filter(grepl("MALASPINA",Study))  %>% 
  select(Sample.label,Study,Date,Sampling.depth.m,Size.fraction.lower.threshold.micrometre,Size.fraction.upper.threshold.micrometre,
        Latitude,Longitude,Ocean.and.sea.regions,Mean_Temperature.deg.C.,Mean_Salinity.PSU.,Mean_Oxygen.umol.kg.,Land_Dist)  %>%
  dplyr::rename(sample_ID=Sample.label,Depth.env=Sampling.depth.m,Temp.=Mean_Temperature.deg.C.,Oxygen=Mean_Oxygen.umol.kg.,Salinity=Mean_Salinity.PSU.,Locality=Ocean.and.sea.regions) %>%
  mutate(filter_size=paste(Size.fraction.lower.threshold.micrometre,".",Size.fraction.upper.threshold.micrometre, sep="")) %>%
  select(sample_ID,Study,Date,Depth.env,filter_size,Latitude,Longitude,Locality,Temp.:Land_Dist)

Menv1<-fuzzyjoin::geo_left_join(Menv1, Msam_amp, method = "haversine", max_dist = 10)

Menv1<-Menv1 %>% 
  dplyr::rename(Latitude=Latitude.x ,Longitude=Longitude.x ) %>% 
  select(sample_ID,Study,Station, Date,Depth.env,filter_size,Latitude,Longitude,Locality,Temp.:Land_Dist)


# names(Menv1)     names(Msam1)    names(Msam)
Msam1<- Msam

Mtax1<- Mtax %>%   
  select(-cruise,-tax_catagory ) %>% 
  left_join(Menv1, by="sample_ID")  %>% 
  left_join(Msam1, by="sample_ID")

mess<-length(Mtax1$sample_ID)
# simplify     names(Mtax1)   names(Ttax1[280:length(names(Ttax1))])
Mtax1<- Mtax1 %>%  
  dplyr::rename(ID=sample_ID,Date.Time=Date,Depth=Depth.env,Paired.reads=num_reads) %>%  
  mutate(Cruise="Malaspina") %>%  ## moved to combined
  mutate(Sample=c(1:mess)) %>% 
  mutate(Paired.reads=as.numeric(gsub(",","",Paired.reads))) %>% 
  select(ID,Sample,Cruise,Station,area,filter_size,gene_domain,pid,Depth,Basin,Date.Time,Latitude:Paired.reads,Hydra.vulgaris:Millepora.sp..EK.2011) %>%     # names(Mtax1)  names(Ttax1[999:1131])   
  mutate(Station=as.character(Station))

#########################################################################################################
###### combine data     ####################################################################
#########################################################################################################

dat<-bind_rows(Ttax1,Mtax1)    ## names(dat)      names(dat[999:1041,])  
dat<- dat %>% 
  select(ID:NPP,Depth:Paired.reads,Hydra.vulgaris:Millepora.sp..EK.2011)# 
### makes sure second dat has same column check 1140


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
#               isolate species and fix names      ########################################################
############################################################################################################
## species start - 20
#spID<-names(dat[29:length(names(dat))])   ### make separate species ID dataframe
sp<-data.frame(names(dat[33:length(names(dat))]) )    ## make separate species dataframe
names(sp)[1]<-paste("species.name")
#######################################
###   fix bad NCBI species names   ####
sp1<-splitstackshape::cSplit(sp, "species.name", sep=".", drop = F)  # names(sp1)
#sp2<-sp1[complete.cases(sp1$species.name_3),]   #
sp_g<-sp1[,1:3]   ### colnames(sp_g)
sp_g$species.name_fix<-paste(sp_g$species.name_1,sp_g$species.name_2)   ### fixed names
names(sp_g)[3]<-paste("species.fixed")
sp_g<-sp_g[,c(1,3,4)]
sp_g$species.name_fix[sp_g$species.name_fix=="Nymphon unguiculatum-charcoti"]<-"Nymphon unguiculatum"
names(sp_g)[2]<-paste("species.only")
names(sp_g)[1]<-paste("species.original.colname")
####   get ncbi taxon ID  --  n2
n2<- n2 %>% 
  mutate(species.original.colname=gsub(" ", "\\.",Category)) %>% 
  mutate(species.original.colname=gsub("-", "\\.",species.original.colname))

sp_g<-sp_g %>% # names(sp_g)
  left_join(n2) 
  
sp_g$ID[sp_g$species.original.colname=="Cerapachys.biroi"]<-2015173
sp_g$species.name_fix[sp_g$species.original.colname=="Cerapachys.biroi"]<-"Ooceraea biroi"
############################################################################################################
### get taxonomy from NCBI############################################################################################
################################################    names(cl)  names(sp_g)
setwd("/Users/geraldn/Dropbox/eDNA_db/reference_data/NCBI_tax_taxonomizr")
tax<-taxonomizr::getTaxonomy(sp_g$ID, "accessionTaxa.sql", desiredTaxa = c("superkingdom","kingdom",
                                                                  "phylum", "class", "order", "family", "genus", "species"), mc.cores = 3,debug = FALSE)
tax<-data.frame(tax,stringsAsFactors = F)
sp_good<-cbind(sp_g,tax)
sp_good<-sp_good %>% 
    filter(species.only != "reads")

######################################################################################################
######################################################################################################
############################################################################################################
##############################################################
############################################################################################################
#######  #####################################################################################################
#####get info from worms  taxonomy and if marine     and if needed put before function
#############
library("worrms")
l<-length(sp_g$Category)   ## cange to name??     1-99,100-199, 200-299
i<-1  #  names(sp_good)
cworm <- wm_records_names(name = sp_good$species.name_fix[i:180] , fuzzy = F, marine_only = F)   #st i<-1
cwormm <- wm_records_names(name = sp_good$species.name_fix[181:194] , fuzzy = F, marine_only = F)  
cworm1 <-bind_rows(cworm)  
cworm2 <-bind_rows(cwormm) 
  cworm2<-rbind(cworm2,cworm1)

cwa<-cworm2[!duplicated(cworm2$valid_name),]  # names(cwa1)
cwa<-cwa[,c(9,3,11:16,19:22,8)]#  simplify
names(cwa)[2]<-paste("species.ncbi_fix") 
names(cwa)[1]<-paste("species_worrms")
names(sp_good)[3]<-paste("species.ncbi_fix") 
sp_good<-left_join(sp_good,cwa,by="species.ncbi_fix")
# make marine 1 and 0          rm(sp_good)
sp_good$Marine<-sp_good$isMarine
sp_good$Marine[sp_good$Marine!=1]<-0
sp_good$Marine[is.na(sp_good$Marine)]<-0
sp_good$Marine[sp_good$isBrackish==1]<-1
######################################################################################################
######################################################################################################
##### make good taxonomy, filll in missing worrms with ncbi    names(sp_good)
sp_good$kingdom<-sp_good$kingdom.y
ind <- is.na(sp_good$kingdom.y)
sp_good$kingdom[ind] <- sp_good$kingdom.x[ind] 
sp_good$phylum<-sp_good$phylum.y
ind <- is.na(sp_good$phylum.y)
sp_good$phylum[ind] <- sp_good$phylum.x[ind] 
sp_good$class<-sp_good$class.y
ind <- is.na(sp_good$class.y)
sp_good$class[ind] <- sp_good$class.x[ind] 
sp_good$order<-sp_good$order.y
ind <- is.na(sp_good$order.y)
sp_good$order[ind] <- sp_good$order.x[ind] 
sp_good$family<-sp_good$family.y
ind <- is.na(sp_good$family.y)
sp_good$family[ind] <- sp_good$family.x[ind] 
sp_good$genus<-sp_good$genus.y
ind <- is.na(sp_good$genus.y)
sp_good$genus[ind] <- sp_good$genus.x[ind] 
# clean up    names(sp_good)
sp_good<-sp_good[,c(1:6,13,20:31)]
##  add others to marine
sp_good$Marine[sp_good$class=="Anthozoa"]<-1

############################################################################################################
############################################################################################################
############################################################################################################
##########################################################################################################
######   make one big dataframe !!!!!!     names(dat) #######################################################################
dat1<- dat %>%
  gather(species.original.colname,reads,Hydra.vulgaris:Millepora.sp..EK.2011)
dat1$reads[is.na(dat1$reads)] <- 0
mess<-paste(dat1$ID,"_",dat1$Sample,"_",dat1$filter_crit)
dat1<- data.frame(cbind(mess,dat1))
names(dat1)[1]<-paste("unique_ID") 
dat2<- dat1 %>%
  left_join(sp_good,by="species.original.colname")
#  names(dat1)
### export
#      
data.table::fwrite(dat2,"/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/DMAP_biomass_apr19.csv",row.names=F, sep=",")
############################################################################################################
##########################################################################################################
############################################################################################################


