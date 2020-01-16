
library(dplyr)

########################################################################
## April 2019 DMAP data for biomass comp  ########################################################################
setwd("/Users/geraldn/Dropbox/Global_databases/DMAP/DMAP_for_biomass")
l<-list.files(path = ".",pattern=".*species.*txt")

for (i in 1:length(l)){
  # i<-1
  x<-read.table(l[i],sep="\t",header=T)
  x1<-data.frame(t(x[,c(3:length(names(x)))]))
  names(x1)<-x[,2]
  sample_ID<-row.names(x1)
  search<-rep(l[i], times=length(sample_ID))
  x1<-cbind(sample_ID,search,x1)
  x1<-x1[!x1$sample_ID=="Total",]
  # get taxaID and name
  nn<-x[1:2]
  
  if (i==1) {
    data<-x1
    nn2<-nn
  }
  else{
    data<-plyr::rbind.fill(data,x1)
    nn2<-rbind(nn2,nn)
  }
}

#### tidy
dat1<-data %>%  # names(data)
  separate(search,into=c("cruise","area","tax_catagory","gene_domain","pid","coverage","mes"),by="_") %>% 
  dplyr::select(-mes)

nn2<-nn2 %>% 
  filter(!duplicated(ID))
#   names(data)    data$filter_crit         data$sample_ID
write.csv(dat1, file = "/Users/geraldn/Dropbox/Global_databases/DMAP/DMAP_for_biomass/all_seq.csv", row.names = F)
write.csv(nn2, file = "/Users/geraldn/Dropbox/Global_databases/DMAP/DMAP_for_biomass/sp_ids.csv", row.names = F)




############################################################################################
############################################################################################
############################################################################################
setwd("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/work")
###   RS
setwd("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/data/seq")
l<-list.files(path = ".")
l3<-l[which(grepl(1, l))]
l6<-l[which(grepl(60, l))]
l9<-l[which(grepl(90, l))]
for (i in 1:length(l)){
  # i<-1
  x<-read.table(l[i],sep="\t",header=T)
  x1<-data.frame(t(x[,c(3:length(names(x)))]))
  names(x1)<-x[,2]
  sample_ID<-row.names(x1)
  filter_crit<-0
  x1<-cbind(sample_ID,filter_crit,x1)
  x1<-x1[!x1$sample_ID=="Total",]
  if (l[i] %in% l3) {
    x1$filter_crit<-30
  } else if (l[i] %in% l6) {
    x1$filter_crit<-60
  } else { 
    x1$filter_crit<-90 } 
  
  if (i==1) {
    data<-x1
  }
  else{
    data<-rbind.fill(data,x1)
  }
}
#   names(data)    data$filter_crit         data$sample_ID
write.csv(data, file = "/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/data/seq/RS_all_seq.csv", row.names = F)

########################################################################
##  tara  ########################################################################
setwd("/Users/geraldn/Dropbox/Global_databases/DMAP/TARA/seq")
l<-list.files(path = ".")
l3<-l[-which(grepl(60, l))]
l3<-l3[-which(grepl(90, l3))]
l6<-l[which(grepl(60, l))]
l9<-l[which(grepl(90, l))]
for (i in 1:length(l)){
  # i<-1
  x<-read.table(l[i],sep="\t",header=T)
  x1<-data.frame(t(x[,c(3:length(names(x)))]))
  names(x1)<-x[,2]
  sample_ID<-row.names(x1)
  filter_crit<-0
  x1<-cbind(sample_ID,filter_crit,x1)
  x1<-x1[!x1$sample_ID=="Total",]
  if (l[i] %in% l3) {
    x1$filter_crit<-30
  } else if (l[i] %in% l6) {
    x1$filter_crit<-60
  } else { 
    x1$filter_crit<-90 } 
  
  if (i==1) {
    data<-x1
  }
  else{
    data<-rbind.fill(data,x1)
  }
}
#   names(data)    data$filter_crit         data$sample_ID
write.csv(data, file = "/Users/geraldn/Dropbox/Global_databases/DMAP/TARA/tara_all_seq.csv", row.names = F)


#############################################################################################################
#############################################################################################################
library(ggmap)
library(raster)
library(xlsx)
library(dplyr)
################## get layer am map
setwd("/Users/geraldn/Dropbox/Global_databases/GMED/land_distance")
ld <- raster(read.asciigrid("gb_land_distance.asc"))
#############################################################################################################
#####      get distance to chore for all locations
RS<- read.csv("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/RSsampledat.csv")  
RS$Land_Dist <- extract(ld,cbind(RS$Longitude, RS$Latitude))
#   TA<-read.delim(file="/Users/geraldn/Dropbox/Global_databases/DMAP/TARA/TARA_2009-2013_ContextEnviro.txt", encoding='UTF-8', skip=10,
 #               header=F)   coudl use but need to split after import use excel pasted from pangea

TA<- xlsx::read.xlsx("/Users/geraldn/Dropbox/Global_databases/DMAP/TARA/tara.stations.good.xlsx",sheetIndex = 1)
#   names(TA)
TA1<- TA %>%
  select(Station, cat, Date.Time,Latitude,Longitude,Bathy.depth..m.,Distance..km.,Depth..nominal.1,Locality,Locality.2,Locality.3,Locality.4)

TA1$Land_Dist <- raster::extract(ld,cbind(TA1$Longitude, TA1$Latitude))


#You can plot it to get an idea of the data:
plot(ld)
points(RS$Longitude, RS$Latitude, pch=4, cex = .6, col="red")
points(TA$Longitude, TA$Latitude, pch=4, cex = .6)

## export
write.table(RS,"/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/RSsampledat1.csv",row.names=F, sep=",")
write.table(TA1,"/Users/geraldn/Dropbox/Global_databases/DMAP/TARA/tara.stations.good1.csv",row.names=F, sep=",")




