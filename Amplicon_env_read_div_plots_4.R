
library(dplyr)
library(tidyr) 

##         amplicons  !!!!!!!!!!!!!!!!!!!  
############################################################################################################
##     tara
source("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/R/projects/Global_ocean_genome_analysis/amplicon_stats_tara_source.R")
########################################################################################################
t_rare<- stat_rare[["read_low_sp"]]  #
###   ! seq column tells if rarified or raw reads
t_rare$seq<-"rare"
t_raw<- stat_raw[["read_low_sp"]]  ## names(t_rare)   names(t_raw)  
t_raw$seq<-"raw"
t_all<- t_rare %>% # names(t_all)   names(t_rare)
  bind_rows(t_raw) %>% 
  dplyr::select(Depth,land_dist,Present.Surface.Temperature.Mean,Present.Surface.Temperature.Range,
         Present.Surface.Primary.productivity.Mean, Present.Surface.Primary.productivity.Range,OHI_2013,
         seq,reads,rich,shan,Sample.ID)  %>% 
  mutate(cruise="Tara")
############################################################################################################
##    mal        summary(stat)
source("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/R/projects/Global_ocean_genome_analysis/amplicon_stats_mal_source.R")

#  1- "raw.otu.DNA_read_low"  2- "raw.species.DNA_read_low" 3- "rare.otu.DNA_read_low"   4- "rare.species.DNA_read_low"
m_raw_dna<-stat[["raw.species.DNA_read_low"]]  ## names(m_raw_dna)   names(m_rare_dna) 
m_rare_dna<-stat[["rare.species.DNA_read_low"]] 
m_raw_rna<-stat[["raw.species.RNA_read_low"]]  ## names(m_raw_rna)   names(stat) 
m_rare_rna<-stat[["rare.species.RNA_read_low"]] 
wanted_df<-c(5,13,7,15)   # set what df in stat you want
tab_names<-names(stat[wanted_df])   # get names
mal<-data.table::rbindlist(stat[wanted_df], use.names=T, fill=T, idcol=T)  # names(mal)
mal_all <-mal %>% 
  dplyr::rename(seq=.id)  %>% 
  dplyr::select(Depth,land_dist,Present.Surface.Temperature.Mean,Present.Surface.Temperature.Range,
         Present.Surface.Primary.productivity.Mean, Present.Surface.Primary.productivity.Range,OHI_2013,
         seq,reads,rich,shan,Sample.ID)  %>% 
  mutate(cruise="Mal") %>% 
  mutate(cruise=replace (cruise,!(grepl("deep",Sample.ID) | grepl("M",Sample.ID)), "Profile") ) %>% 
  mutate(cruise=replace (cruise,grepl("deep",Sample.ID), "Deep") ) 
# combine cruises
  dat<- mal_all %>% 
    bind_rows(t_all)

########################################################################################################
####   plot by 7 variables      ##############################################################################
############################################################################################################
   ## specify data      unique(dat_mar$filter_crit)    unique(dat_mar$filter_size)
   #  names(dat)
   # mess<-dat1[dat1$Latitude < -45,]
   r<-c(1,2)  # dependents - in columns reads, rare,  richness  names(dat)
   f<- c(1:7) # variables - 
   varnames<-c("Depth (m)","Distance to land (km x 100)","Surface temperature mean",
               "Surface temperature range","Surface primary productivity mean", 
               "Surface primary productivity range", "Ocean health index")
   m_raw<-tab_names[c(1,3)]
   m_rare<-tab_names[c(2,4)]
    m_dna<-tab_names[c(1,2)]
    m_rna<-tab_names[c(3,4)]
####  no pattern in reads for mar, check terr and then remove     #    Tara_Mar_sp_30f   unique(dat1$Cruise)
    # j<- 11; i<-1   # tab_names
    # tables of significant
    # tara, mal, profile r
    q<-1
    std<-matrix(2, nrow = length(f), ncol = length(q))
    std[c(4)]<-1 #  make 1 if singificant
    smd<-matrix(2, nrow = length(f), ncol = length(q))
    smd[c(2,4,6)]<-1
    smp<-matrix(2, nrow = length(f), ncol = length(q))
    smp[c(2)]<-1
    smr<-matrix(2, nrow = length(f), ncol = length(q))
    smr[c(3,4,6)]<-1
    smdd<-matrix(2, nrow = length(f), ncol = length(q))
    #smdd[c(2,4,6)]<-1
############################################################################################
   par(mfcol=c(length(f),2), mar=c(3.1,3,.5,1), oma=c(2,1,1,0)) 
    sss<-1
   for (j in r){   # dependent
     
     
     for (i in f){  # factor
       # set up for 3 columns
       # removed so everything is rarified
      # if (j==9){      
       #  dat1<-dat %>% 
       #     dplyr::filter(seq=="raw"| seq %in% m_raw)
      # }  else{
         dat1<-dat %>% 
           dplyr::filter(seq=="rare"| seq %in% m_rare) %>% 
            filter(complete.cases(Depth))
         #}
       # set up for plot
       y<-data.frame(dat1[,10])
       x<-data.frame(dat1[,i])
       y<-y[complete.cases(x[,1]),]
       x<-x[complete.cases(x[,1]),]
 # set muliplying for mal
         yyl<-c(0,20); mt<-10; mm<-1

      
       # plot depending on r -  diffent colors, tara dna, mal dna, mal rna
       if (j==1){
         yr<-dat1[dat1$cruise=="Tara",10] / mt # unique(dat$cruise)  psych::describe(xm)
         xr<-dat1[dat1$cruise=="Tara",i]
         # all mal
         ym<-dat1[dat1$seq %in% m_dna,10] / mm  # for mal dna
         xm<-dat1[dat1$seq %in% m_dna,i]
         
         plot(x,y,pch=16,xlab="", ylab="",col="white", ylim=yyl, cex=0.1)
         points(xr,yr,pch=20,col="#fd8d3c", cex=0.2)  # "#fd8d3c","#a6bddb","#3690c0" from map orange and blue
         points(xm,ym,pch=20,col="#3690c0", cex=0.2)     
         
         ##tara
         wd1<-2
         abline(lm(yr~xr), col="#fd8d3c",lty=std[i,sss],lwd=1)# Add linear regression line
         # mal dna
         abline(lm(ym~xm),col="#3690c0",lty=smd[i,sss],lwd=1)# Add linear regression line
       } 
       
       if (j==2){   
       # deep mal  unique(dat1$cruise)
         yr<-dat1[dat1$seq %in% m_dna &dat1$cruise=="Deep",10] / mm  # for mal dna
         xr<-dat1[dat1$seq %in% m_dna & dat1$cruise=="Deep",i]
       # prof mal
       ymp<-dat1[dat1$seq %in% m_dna &dat1$cruise=="Profile",10] / mm  # for mal dna
       xmp<-dat1[dat1$seq %in% m_dna & dat1$cruise=="Profile",i]
       # mall rna
       ymr<-dat1[dat1$seq %in% m_rna,10]   # for mal rna
       xmr<-dat1[dat1$seq %in% m_rna,i]
       range(yr);range(ym);range(ymr) #, 153,23,14
       
       plot(x,y,pch=16,xlab="", ylab="",col="white", ylim=yyl, cex=0.1)
       points(xr,yr,pch=20,col="grey20", cex=0.2)  # "#fd8d3c","#a6bddb","#3690c0" from map orange and blue
       points(xmp,ymp,pch=20,col="blue", cex=0.2)     
       points(xmr,ymr,pch=20,col="grey50", cex=0.2) #  c("#628DDE") between blue
       # mal profile dna
       abline(lm(yr~xr),col="grey20",lty=smdd[i,sss],lwd=1)# Add linear regression lin
       # mal profile dna
       abline(lm(ymp~xmp),col="blue",lty=smp[i,sss],lwd=1)# Add linear regression lin
       
       # mal rna
       abline(lm(ymr~xmr),col="grey50",lty=smr[i,sss],lwd=1.1)# Add linear regression line #747F82c
       
       }
       
     
         mtext(varnames[i], side=1,line=2.0, cex=.8, adj=.5,outer=F)

       
     }
     #sss<-sss+1
   }
   ex<-.8   # temp, depth, latitude, chloro 
   mtext("Richness (Tara/10)", WEST<-2,las=0,line=-.7, cex=ex, at=.5 ,outer=TRUE)
   mtext("Richness", WEST<-2,las=0,line=-29.8, cex=ex, at=.5 ,outer=TRUE)
   #mtext("Temperature (degrees C)", SOUTH<-1,line=-35., cex=ex, at=.5,outer=TRUE)
   #### save as 743 width  991 height    
   
