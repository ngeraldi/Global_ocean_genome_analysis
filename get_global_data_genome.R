

library(dplyr)
library(tidyr)

library(ncdf4)  
library(raster)
library(stringr)
library(mregions) 
library(rgeos)
library(rgdal)


dat_genome<-data.table::fread(file = '/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/DMAP_biomass_apr19.csv', sep = ',', header = TRUE)
data<-dat_genome %>% 
  dplyr::select(Latitude,Longitude) %>% 
  dplyr::distinct(Latitude,Longitude, .keep_all = T)

## nead-- human impact; land-dist; primary prod; SST SST variance: human in 100, 500, 1000; biodiversity
##############################################################################################
#####   temp and productivity from !!!!!!!bio_oracle!!!!!!!
setwd("/Users/geraldn/Dropbox/Global_databases/Bio-oracle/Surface_present")   
files<-list.files(patter="*.asc")
#####   get names from file name
namess<-as.data.frame(files)
namess<- namess %>%
  dplyr::rename(file_name=files) %>% 
  mutate(file_name=as.character(file_name))
namess$file_name<-substr(namess$file_name,1,nchar(namess$file_name)-4)
### specify files if needed
ff<-c(23,25,12,14)  #  c(6:16)
f<-files[ff]
cn<-namess$file_name[ff]
###  start loop    i<-1
for (i in 1:length(f)){
  nam<-paste(cn[i])
  x<- raster(read.asciigrid(f[i]))  # names(data)
  data$loopname<- extract(x ,cbind(data$Longitude, data$Latitude))
  q<-length(names(data))-1  # use for names
  names(data) <- c(names(data[,1:q]), paste(nam)) 
}
####### human impact
##   !!!! need to reproject   !!!!!!!!!!!
#  overall OHI  #########
wgs<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#define the mollweide projection coordinate reference system (crs)
moll=crs('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

setwd("/Users/geraldn/Dropbox/Global_databases/OHI/cumulative_impact_one_2013_global_cumul_impact_2013_mol_20150714053146")
OHI_2013<- raster("global_cumul_impact_2013_all_layers_wgs.tif")     #   hist(OHI_2013)  plot(OHI_2013)
#reproject the shapefile to mollweide using spTransform  did in July 2018, dont need to do again
# OHI_2013_wgs <- projectRaster(OHI_2013, crs=wgs)
##   writeRaster(OHI_2013_wgs, filename="global_cumul_impact_2013_all_layers_wgs.tif", format="GTiff", overwrite=TRUE)
#  extract
data$OHI_2013<- extract(OHI_2013 ,cbind(data$Longitude, data$Latitude))   # hist(data$OHI_2013)
### mes with mal deep 81
mes<- extract(OHI_2013 ,cbind(178.5, -28.4)) ######## -28.406167, 179.1413333  14.582000
70.01100
##  gmed depth
setwd("/Users/geraldn/Dropbox/Global_databases/GMED/depth")
global_depth <- raster(read.asciigrid("gb_depth.asc"))
data$global_depth <- extract(global_depth,cbind(data$Longitude, data$Latitude))
#######  land dist
setwd("/Users/geraldn/Dropbox/Global_databases/GMED/land_distance")
land_dist <- raster(read.asciigrid("gb_land_distance.asc"))
data$land_dist <- extract(land_dist,cbind(data$Longitude, data$Latitude))
data1<-data[complete.cases(data$Latitude),]
######  biodiverstiy
shape <- readOGR(dsn = "/Users/geraldn/Dropbox/Global_databases/global_diversity/OBIS_summaries", layer = "summaries") 
#   plot(shape)    summary(shape)     spplot(shape, z="shannon")
crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")  # geographical, datum WGS84 copy from summary
pt <-dplyr::select(data, Longitude, Latitude)
pt<-pt[complete.cases(pt$Longitude),]
coordinates(pt) <- c("Longitude","Latitude")  # se
proj4string(pt) <- crs.geo  
###   pt<-SpatialPoints(pt)  if didn't use previous code
e<-over(pt,shape)  #  names(e)
data1<-bind_cols(data1,e[,c(4,5)])
########################################
#####   human population 
inv_rotate <- function(r) {
  xmin(r) <- 0
  xmax(r) <- 360
  r <- rotate(r)
  xmin(r) <- 0
  xmax(r) <- 360
  r
}
wgs<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
moll=crs('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
#################
setwd("/Users/geraldn/Dropbox/Global_databases/human_pop/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2010")
pop <- raster("gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2010.tif")
pop_rot <- inv_rotate(pop)
# using extract
p100k<-raster::extract(pop, cbind(data1$Longitude, data1$Latitude), buffer = 100000 , fun = sum ,df=T, na.rm=TRUE)# buffer is in m
p500k<-raster::extract(pop, cbind(data1$Longitude, data1$Latitude), buffer = 500000 , fun = sum ,df=T, na.rm=TRUE)# buffer is in m
p1000k<-raster::extract(pop, cbind(data1$Longitude, data1$Latitude), buffer = 1000000 , fun = sum ,df=T, na.rm=TRUE)# buffer is in m
data2 <-data1 %>%
  bind_cols(p100k) %>%
  dplyr::select(-ID) %>%
  dplyr::rename(Pop_in_100km=gpw.v4.population.count.adjusted.to.2015.unwpp.country.totals_2010) %>%
  bind_cols(p500k) %>%
  dplyr::select(-ID) %>%
  dplyr::rename(Pop_in_500km=gpw.v4.population.count.adjusted.to.2015.unwpp.country.totals_2010) %>%
  bind_cols(p1000k) %>%
  dplyr::select(-ID) %>%
  dplyr::rename(Pop_in_1000km=gpw.v4.population.count.adjusted.to.2015.unwpp.country.totals_2010) %>% 
  mutate(Pop_in_100km=replace_na(Pop_in_100km,0), Pop_in_500km=replace_na(Pop_in_500km,0) ,Pop_in_1000km=replace_na(Pop_in_1000km,0)) 
######     get  meow region data      ##################################
##  only coastal regions do not use

#  export table set name depending on which oracle database
setwd("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV")
write.table(data2,"Global_layers_April18.csv", sep=",",row.names=F)


############################################################################################
#########################################
######  get lohghurst biomes    ##################
setwd("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV")
data<-read.table("Global_layers_April18.csv", sep=",",header=T)
shape <- readOGR(dsn = "/Users/geraldn/Dropbox/Global_databases/regions/longhurst_v4_2010", layer = "Longhurst_world_v4_2010") 
#   plot(shape)    summary(shape)     spplot(shape, z="shannon")
crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")  # geographical, datum WGS84 copy from summary
pt <-dplyr::select(data, Longitude, Latitude)
pt<-pt[complete.cases(pt$Longitude),]
coordinates(pt) <- c("Longitude","Latitude")  # se
proj4string(pt) <- crs.geo  
###   pt<-SpatialPoints(pt)  if didn't use previous code
e<-over(pt,shape)  #  names(e)
names(e)[1]<-"lohghurst_biome"
data1<-bind_cols(data,e)
# match ocean
ob_cat<-read.csv(file="/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/lohghurst_biome.csv")
names(ob_cat)[1]<-"lohghurst_biome"
ob_cat<-ob_cat[,c(1,3)]
data2<-data1 %>%    # names(data2)     unique(data2$ocean)
  mutate(lohghurst_biome=as.character(lohghurst_biome)) %>% 
    mutate(lohghurst_biome=replace(lohghurst_biome,lohghurst_biome=="CHIL","HUMB")) %>% 
    left_join(ob_cat) 
## fix malaspina below australia that say antarctic ocean
data2$ocean[data2$ocean=="Antarctic" & data2$Cruise.x=="Mal" & data2$Longitude>145]<-"Pacific"
data2$ocean[data2$ocean=="Antarctic" & data2$Cruise.x=="Mal" & data2$Longitude<145 & data2$Longitude>100]<-"Indian"

#  save
setwd("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV")
write.table(data2,"Global_layers_genome.csv", sep=",",row.names=F)

#data<-read.table("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/DMAP/R/CSV/Global_layers_oct18.csv", sep=",",header=T)


