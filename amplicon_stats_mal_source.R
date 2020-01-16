
library(dplyr)
library(tidyr)
library(vegan)
# library(edgeR)

## function to flatten list and keep names
flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}
############################################################################################################
##   get data
dat<-data.table::fread(file = '/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/groups/amplicon_Mal_metazoans.csv', sep = ',', header = TRUE)
names(dat)   # unique(dat$Cruise.x)
dat<-dat[,1:38]  # added metazoan per sample, but messed things up, removed
### get global geo data
geo<-read.table("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/Global_layers_April18.csv", sep=",", header = TRUE)
##  get only uniqiue lat and long    names(geo)
geo1<-geo %>% 
  distinct(Latitude, Longitude, .keep_all = TRUE) %>% 
  mutate(Depth= replace(Depth, which(is.na(Depth) & Depth_region=="SRF"), 5)) %>% # not good when join later
  mutate(land_dist= replace(land_dist, which(is.na(land_dist) & Station=="TARA_011"), 0.01)) %>% 
  mutate(OHI_2013= replace(OHI_2013, which(is.na(OHI_2013) & Station==81), 3.77)) %>% 
  dplyr::rename(Depth_geo=Depth)
########################################################################################################
######  data
pat<-c("Annelida","Arthropoda","Chordata","Cnidaria","Ctenophora","Echinodermata","Gastrotricha","Mollusca",
       "Nematoda","Orthonectida","Platyhelminthes","Porifera","Rotifera","Sipuncula","Tardigrada","Urochordata","Chaetognatha",
       "Craniata","Nemertea","Bryozoa","Hemichordata","Brachiopoda","Entoprocta","Cephalochordata","Entoprocta") 
pat<-paste(pat, collapse = "|")
#######   get basics   ######################################################
mes<-dat %>% 
    group_by(Station) %>% 
    summarise(mean=mean(reads),sum=sum(reads),n(), sum_rar=sum(reads_rarified))
##  make new consistant column names
dat$read_col<-dat$reads
dat$rare_col<-dat$reads_rarified
dat$perc_col<-dat$reads_percent
dat$otu_col<-dat$taxID
dat$sp_col<-as.character(dat$sp_num)
dat$split1<-dat$DNAID
#
dat$pid_col<-dat$PI_match ## set pid filter criteria
pid_high<-97  # strict percent ID
pid_low<-90   # less strict  percent ID
#  names(dat)
##    !!!!!!!!!     neeed to check if change columns    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
samples_col<-c(2,3,8:22)  ## unique to samples keep these
data_cols<-c(39:45)  ## unique to taxa

######  initial filtering and select needed collumns
dat2<- dat    %>%        ### head(tdat)    names(dat)
  select_at(c(samples_col,data_cols))
  #filter(Cruise.x== "TARA") %>%
  #filter(filter_size_cat != "0.8-20", filter_size_cat != "<0.8") %>% 
  ##  filter(DNAID=="DNA")   ## include as variables
## split if necessary
dat3 <- split( dat2 , f = dat2$split1 )    ## names(dat2)
  ###   column names should be reads and reads_rarified  !!!!!!!!!  change above if not
##  raw reads
#  rm(dat) ####  remove big dataframe   names(dat3)
datl<-list()
datl[["DNA_read_low"]]<- dat3[[1]] %>%    ### 69212
  filter(pid_col>pid_low) %>% 
  filter(read_col>0) 
datl[["DNA_read_high"]]<- dat3[[1]] %>%    ### 16133 records
  filter(pid_col>pid_high)  %>% 
  filter(read_col>0)
datl[["RNA_read_low"]]<- dat3[[2]] %>%    ### 69212
  filter(pid_col>pid_low) %>% 
  filter(read_col>0)  
datl[["RNA_read_high"]]<- dat3[[2]] %>%    ### 16133 records
  filter(pid_col>pid_high)  %>% 
  filter(read_col>0)
###  summary(datl)    mes<-datl[["DNA_read_low"]]
##### !!!!!!!!!!!!!!!!!!!!!!! make new with species  !!!!!!!!!!!!!!!!!!!!!!!! names(datl[[1]])
samples_col<-c(1:17,22)  ## unique to samples keep these
data_cols<-c(18,19,20)  ## unique to taxa
#
datl_sp<-lapply(datl,function(x) {
  x %>% 
  group_by_at(c(samples_col)) %>% 
  summarize(reads=sum(read_col),reads_rarified=sum(rare_col), reads_percent=sum(perc_col)) %>% 
  mutate(otu=sp_col)
})
### clean up otu based table !!!!!!!!!!!!!!!!! names(datl[[1]])
samples_col<-c(1:17,21)  ## unique to samples keep these
data_cols<-c(18,19,20)  ## unique to taxa
datl_otu<-lapply(datl,function(x) {
  x %>% 
    group_by_at(c(samples_col)) %>% 
    summarize(reads=sum(read_col),reads_rarified=sum(rare_col), reads_percent=sum(perc_col)) %>% 
    mutate(otu=otu_col)
})
##
datl2<-list(datl_otu,datl_sp)  ## names(datl2[[1]])
names(datl2)<-c("otu","species")
datl3<-flattenlist(datl2) ## names(datl3)         names(datl3[[1]])
#   otu or species, DNA or RNA, 90 or 97 pid
##### spread   #######   names(dat_otu)  summary(dat_otu)   names(dat$read_low)
sel_raw<-c(1:17,22,19)
sel_rar<-c(1:17,22,20)
dat_otu_reads<-lapply(datl3,function(x) {
  x %>% 
    ungroup() %>% 
    select_at(sel_raw) %>% 
    spread(otu, reads, fill=0)})
### for rarefied      names(dat_otu_rare) summary(dat_otu_rare)    summary(dat_otu_reads)
#mes2<-datl3[[8]]
mes1<-dat_otu_reads[[8]]
#   mes<-dat_otu_rare[[8]]
dat_otu_rare<-lapply(datl3,function(x) {
  x %>% 
    ungroup() %>% 
    filter(reads_rarified>0) %>%
    select_at(sel_rar) %>% 
    spread(otu, reads_rarified, fill=0) 
  })
#  combine and split
datl2<-list(dat_otu_reads,dat_otu_rare)  ## names(datl2[[1]])
names(datl2)<-c("raw","rare")
datl3<-flattenlist(datl2) ## names(datl3)    ## names(datl3[[1]])
### calculate  riches and div - simpsons
#####    for all otus  --  otu     lapply(dat_otu, names)
nn<-18 ## first column of species!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##### raw
stat<-lapply(datl3,function(x) {
  shan <- as.numeric(diversity(x[,nn:ncol(x)]),index="simpson") # less biased to different n-Gihrig 2011
  rich <- as.numeric(rowSums(x[,nn:ncol(x)] > 0))  # plot(sample_rich)
  reads <- as.numeric(rowSums(x[,nn:ncol(x)]))
  # clean and join to geo   names(geo)  names(tdat90_otu[,1:30])    names(tdat90_otu)
  x %>% 
    bind_cols(shan=shan,rich=rich,reads=reads) %>% 
    mutate(Latitude_abs=abs(Latitude)) %>% 
    left_join(geo1[,c(6:19)]) # names(geo1)
})   ###  summary(stat_raw)
summary(stat)   ## 231 samples   mes<-stat$raw.otu.DNA_read_low
# m<-lapply(stat, dim)
#  do.call(rbind,m)
########################################################################################################
########################################################################################################
