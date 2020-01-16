
library(dplyr)
library(tidyr)
library(vegan)
library(lme4)
#library(MASS)
library(car)

source("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/R/projects/Global_ocean_genome_analysis/amplicon_stats_mal_source.R")
##   use stat
summary(stat)
########################################################################################################
########################################################################################################
########################################################################################################
#### stats  ######
#
#  1- "raw.otu.DNA_read_low"  2- "raw.species.DNA_read_low" 3- "rare.otu.DNA_read_low"   4- "rare.species.DNA_read_low"
sdat<-stat[["rare.species.DNA_read_low"]]  ## names(stat[[1]])   names(stat)    names(sdat)
########################################################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# if want profile only  ########
sdat <- sdat  %>% 
  filter(grepl("deep",Sample.ID))
########################################################################################################

###scale scale=TRUE divides by the standanrd deviation default does both center and scale
sdat$sTemp.<-scale(sdat$Temp.)
sdat$sDepth<-scale(sdat$Depth)
sdat$sLand_Dist<-scale(sdat$land_dist)
sdat$sLatitude<-scale(sdat$Latitude)
#sdat$sChlorophyll<-scale(sdat$Chlorophyll)
sdat$sLatitude_abs<-scale(sdat$Latitude_abs)
sdat$sPresent.Surface.Temperature.Mean<-scale(sdat$Present.Surface.Temperature.Mean)
sdat$sPresent.Surface.Temperature.Range<-scale(sdat$Present.Surface.Temperature.Range)
sdat$sPresent.Surface.Primary.productivity.Mean<-scale(sdat$Present.Surface.Primary.productivity.Mean)
sdat$sPresent.Surface.Primary.productivity.Range<-scale(sdat$Present.Surface.Primary.productivity.Range)
sdat$sPop_in_100km.Mean<-scale(sdat$Pop_in_100km)
sdat$sPop_in_1000km<-scale(sdat$Pop_in_1000km)
sdat$sOHI_2013<-scale(sdat$OHI_2013)
# fix for deap
sdat$lohghurst_biome <-sdat$Full.name
### plot    untransformed similar  slope=.0000136, signigicant
plot(log(sdat$rich)~sqrt(sdat$reads))
abline(lm(log(sdat$rich)~sqrt(sdat$reads)))  ## summary(lm(sdat$rich~sdat$reads))
#   total reads
plot(log(sdat$rich)~sqrt(sdat$total.read.sample))
abline(lm(log(sdat$rich)~sqrt(sdat$total.read.sample)))  #
#######################################################################################
###   mes around
x<-sdat$reads   # log best
x<-sdat$rich   # no trans best
x<-sdat$shan    # sqrt best    profile only -raw
hist(x)
hist(sqrt(x))
hist(log(x))
shapiro.test(sqrt(x))  #  shapiro.test(log(x+1))
leveneTest(x ~ Crab*Pods*Ulva*Kelp, data=benth)

############################################################################################################
########################################################################################################
####   spatial stats  names(sdat)
library(spdep)
#Create a k=4 nearest neighbor set
sdat<-sdat[complete.cases(sdat$Latitude),]
ssdat<-sdat
sp::coordinates(ssdat) <- ~Longitude+Latitude
#
us.nb4<-knearneigh(coordinates(ssdat), k=8)
us.nb4<-knn2nb(us.nb4)
us.wt4<-nb2listw(us.nb4, style="W")
#
########################################################################################################
########################################################################################################
########################################################################################################
###  richness    names(sdat)     names(sdat2) head

lmer.mod<-lmer(rich~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
    sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
    sPresent.Surface.Primary.productivity.Range+sOHI_2013+(1|lohghurst_biome/Station),data=sdat)##
#
#   mes<- c(AIC(lmer.mod),AIC(lmer.mod2))
mod<-lm(rich~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
          sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
          sPresent.Surface.Primary.productivity.Range+sOHI_2013,ssdat)
#SAR - Lag model
fit.lag<-lagsarlm(mod, ssdat, listw=us.wt4, type="lag", method="MC")
#  summary(fit.lag, Nagelkerke=T)
#SAR - Error model
fit.err<-errorsarlm(mod, ssdat, listw=us.wt4, etype="error", method="MC")
#    summary(fit.err, Nagelkerke=T)
#Spatial Durbin Model
fit.durb<-lagsarlm(mod, ssdat, listw=us.wt4, type="mixed", method="MC")
#   summary(fit.durb, Nagelkerke=T)
#Spatial Durbin Error Model
fit.errdurb<-errorsarlm(mod, ssdat, listw=us.wt4, etype="emixed", method="MC")
#   summary(fit.errdurb, Nagelkerke=T)
#SAC Model  
fit.sac<-sacsarlm(mod, ssdat, listw=us.wt4, type="sac", method="MC")
#    summary(fit.sac, Nagelkerke=T)
#SMA model  
fit.sma<-spautolm(mod, ssdat, listw=us.wt4, family="SMA")
#   summary(fit.sma, Nagelkerke=T)
AICs<-c(AIC(lmer.mod),AIC(mod),AIC(fit.lag), AIC(fit.err), AIC(fit.durb), AIC(fit.errdurb), AIC(fit.sac), AIC(fit.sma))
AICs
#  fit.sma best #  best
###################  test interactions
#           #   summary(fit.sma, Nagelkerke=T)
### check if interactions improve model
mod2<-lm(rich+1~(sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                 sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                 sPresent.Surface.Primary.productivity.Range+sOHI_2013)^2,data=ssdat)##
mod2<-lm(rich~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                         sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                         sPresent.Surface.Primary.productivity.Range+sOHI_2013
           
           ,data=ssdat)

fit.sma<-spautolm(mod2, ssdat, listw=us.wt4, family="SMA")
summary(fit.durb, Nagelkerke=T)
AIC(mod2)   ## 210
Anova(mod2)

piecewiseSEM::sem.model.fits(mod2)     
#piecewiseSEM::sem.fit(mod2,sdat) 
########################################################################################################
########################################################################################################
###  reads

lmer.mod<-lmer(log(reads+1)~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                 sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                 sPresent.Surface.Primary.productivity.Range+sOHI_2013+(1|lohghurst_biome/Station), data=sdat)##
#   summary(lmer.mod)

mod<-lm(log(reads+1)~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
          sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
          sPresent.Surface.Primary.productivity.Range+sOHI_2013,ssdat)
#SAR - Lag model
fit.lag<-lagsarlm(mod, ssdat, listw=us.wt4, type="lag", method="MC")
#  summary(fit.lag, Nagelkerke=T)
#SAR - Error model
fit.err<-errorsarlm(mod, ssdat, listw=us.wt4, etype="error", method="MC")
#    summary(fit.err, Nagelkerke=T)
#Spatial Durbin Model
fit.durb<-lagsarlm(mod, ssdat, listw=us.wt4, type="mixed", method="MC")
#   summary(fit.durb, Nagelkerke=T)
#Spatial Durbin Error Model
fit.errdurb<-errorsarlm(mod, ssdat, listw=us.wt4, etype="emixed", method="MC")
#   summary(fit.errdurb, Nagelkerke=T)
#SAC Model  
fit.sac<-sacsarlm(mod, ssdat, listw=us.wt4, type="sac", method="MC")
#    summary(fit.sac, Nagelkerke=T)
#SMA model  
fit.sma<-spautolm(mod, ssdat, listw=us.wt4, family="SMA")
#   summary(fit.sma, Nagelkerke=T)
AICs<-c(AIC(lmer.mod),AIC(mod),AIC(fit.lag), AIC(fit.err), AIC(fit.durb), AIC(fit.errdurb), AIC(fit.sac), AIC(fit.sma))
AICs   # fit.sma<-spautolm(mod, ssdat, listw=us.wt4, family="SMA")
###################  test interactions
# for reads lmer model best 2077  included interactions higer AIC
### check if interactions improve model


mod_in<-lm(log(reads+1)~(sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                          sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                          sPresent.Surface.Primary.productivity.Range+sOHI_2013)^2,ssdat)

mod_in<-lm(log(reads+1)~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
             sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
             sPresent.Surface.Primary.productivity.Range+sOHI_2013
           +sDepth:sPresent.Surface.Primary.productivity.Range
           +sDepth:sPresent.Surface.Temperature.Range  
             , ssdat)
# no interaction AIC-181 with AIC-250

fit.sma2<-spautolm(mod_in, ssdat, listw=us.wt4, family="SMA")

#AIC(fit.sma2)
summary(fit.sma, Nagelkerke=T)
#  chekc interaction
jtools::interact_plot(lmer.mod, pred = "sPresent.Surface.Temperature.Range", modx = "sDepth", 
                      plot.points = TRUE,interval = TRUE, data=sdat,
                      int.width = 0.95,
                      y.label= "Effect size")  #
jtools::interact_plot(fit.sma2, pred = "sPresent.Surface.Temperature.Range", modx = "sDepth", 
                      plot.points = TRUE,interval = TRUE,
                      int.width = 0.95,
                      y.label= "Effect size")  #
########################################################################################################
###  diversity          x<-sdat$shan    # raw
lmer.mod<-lmer(shan~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                 sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                 sPresent.Surface.Primary.productivity.Range+sOHI_2013+(1|lohghurst_biome/Station), data=sdat)###
#
#   mes<- c(AIC(lmer.mod),AIC(lmer.mod2))
mod<-lm(shan~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
          sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
          sPresent.Surface.Primary.productivity.Range+sOHI_2013,ssdat)
#SAR - Lag model
fit.lag<-lagsarlm(mod, ssdat, listw=us.wt4, type="lag", method="MC")
#  summary(fit.lag, Nagelkerke=T)
#SAR - Error model
fit.err<-errorsarlm(mod, ssdat, listw=us.wt4, etype="error", method="MC")
#    summary(fit.err, Nagelkerke=T)
#Spatial Durbin Model
fit.durb<-lagsarlm(mod, ssdat, listw=us.wt4, type="mixed", method="MC")
#   summary(fit.durb, Nagelkerke=T)
#Spatial Durbin Error Model
fit.errdurb<-errorsarlm(mod, ssdat, listw=us.wt4, etype="emixed", method="MC")
#   summary(fit.errdurb, Nagelkerke=T)
#SAC Model  
fit.sac<-sacsarlm(mod, ssdat, listw=us.wt4, type="sac", method="MC")
#    summary(fit.sac, Nagelkerke=T)
#SMA model  
fit.sma<-spautolm(mod, ssdat, listw=us.wt4, family="SMA")
#   summary(fit.sma, Nagelkerke=T)    Anova(fit.sma)
AICs<-c(AIC(lmer.mod),AIC(mod),AIC(fit.lag), AIC(fit.err), AIC(fit.durb), AIC(fit.errdurb), AIC(fit.sac), AIC(fit.sma))
AICs
###################  test interactions
# for div sma model best 54 AIC and r2 .06
### check if interactions improve model
mod_in<-lm(shan~(sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                  sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                  sPresent.Surface.Primary.productivity.Range+sOHI_2013)^2,ssdat)
mod_in<-lm(shan~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                          sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                          sPresent.Surface.Primary.productivity.Range+sOHI_2013
          
             , ssdat)

fit.sma2<-errorsarlm(mod_in, ssdat, listw=us.wt4, etype="emixed", method="MC")
#
AIC(fit.sma2)
summary(fit.sma2, Nagelkerke=T)
## with interactions  aic 334 and r2 is .27
piecewiseSEM::sem.model.fits(mod_in)
############################################################################################################
########################################################################################################
mod_in<-lm(shan~filter_size_cat+sTemp.+sDepth+sChlorophyll+sLatitude+sDepth:sLatitude+sChlorophyll:sLatitude,ssdat)
fit.sma2<-spautolm(mod_in, ssdat, listw=us.wt4, family="SMA")
##check assumptions
mod<-mod_in
hist(resid(mod))
plot(fitted(mod),resid(mod))
abline(0,0)
qqnorm(resid(mod))
qqline(resid(mod))
plot(mod)###residuals should be random (x=fitted value, y=standardixed res)
plot(mod@frame$sqrt(reads)~fitted(mod))
plot(md$reads,fitted(mod))###resp linear fudncion of teh fitted (x=fitted value, y=response)
qqnorm(mod~resid(mod)|md$reads)##errors are close to normally dist in blocks-linear in plots(x-residuals, y-quantiles of standard norm)
############################################################################################################
########################################################################################################
############################################################################################################
########################################################################################################
## alter filter cat to check significance
sdat<- tdat1 %>%   ## names(sdat)   names(tdat1)
  filter(filter_size_cat != "0.8-20", filter_size_cat != "<0.8")  %>% 
  select(Sample.ID:total.read.sample,filter_size_cat,shan:Latitude_abs)%>% 
  mutate(filter_size_cat= factor(filter_size_cat, levels=c(".5-5","5-20","20-200","180-2000"))) %>% 
  mutate(filter_size_cat2= factor(filter_size_cat, levels=c("5-20","20-200","180-2000",".5-5"))) %>% 
  mutate(filter_size_cat3= factor(filter_size_cat, levels=c("20-200","180-2000",".5-5","5-20"))) %>% 
  mutate(filter_size_cat4= factor(filter_size_cat, levels=c("180-2000",".5-5","5-20","20-200")))
sdat$sTemp.<-scale(sdat$Temp.)
sdat$sDepth<-scale(sdat$Depth)
sdat$sLand_Dist<-scale(sdat$Land_Dist)
sdat$sLatitude<-scale(sdat$Latitude)
sdat$sChlorophyll<-scale(sdat$Chlorophyll)
sdat$sLatitude_abs<-scale(sdat$Latitude_abs)
ssdat<-sdat
sp::coordinates(ssdat) <- ~Longitude+Latitude
us.nb4<-knearneigh(coordinates(ssdat), k=8)
us.nb4<-knn2nb(us.nb4)
us.wt4<-nb2listw(us.nb4, style="W")
#####     reads   #######
mod<-lmer(sqrt(reads)~filter_size_cat+sTemp.+sDepth+sChlorophyll+sLatitude+filter_size_cat:sTemp.+filter_size_cat:sLatitude+sDepth:sLatitude+(1|Province/Station/Depth),data=sdat)
summary(mod)
a<-summary(multcomp::glht(mod, linfct=multcomp::mcp(filter_size_cat="Tukey")))

#####     rich    #######
mod_in<-lm(log(rich)~filter_size_cat+sTemp.+sDepth+sChlorophyll+sLatitude+filter_size_cat:sDepth+filter_size_cat:sTemp.,ssdat)
mod_in<-lm(log(rich)~filter_size_cat2+sTemp.+sDepth+sChlorophyll+sLatitude+filter_size_cat:sDepth+filter_size_cat:sTemp.,ssdat)
mod_in<-lm(log(rich)~filter_size_cat3+sTemp.+sDepth+sChlorophyll+sLatitude+filter_size_cat:sDepth+filter_size_cat:sTemp.,ssdat)
fit.sma2<-spautolm(mod_in, ssdat, listw=us.wt4, family="SMA")
AIC(fit.sma2)
summary(fit.sma2, Nagelkerke=T)
#####    div    #######
mod_in<-lm(shan~filter_size_cat4+sTemp.+sDepth+sChlorophyll+sLatitude+filter_size_cat:sDepth+filter_size_cat:sLatitude+sTemp.:sChlorophyll+sDepth:sLatitude,ssdat)
fit.sma2<-spautolm(mod_in, ssdat, listw=us.wt4, family="SMA")
AIC(fit.sma2)
summary(fit.sma2, Nagelkerke=T)
############################################################################################################
########################################################################################################
