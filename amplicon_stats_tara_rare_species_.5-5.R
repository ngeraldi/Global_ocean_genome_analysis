
library(dplyr)
library(tidyr)
# library(vegan)
# library(edgeR)

source("/Users/geraldn/Dropbox/Documents/KAUST/eDNA/R/projects/Global_ocean_genome_analysis/amplicon_stats_tara_source.R")
## if memory error -  close r, nthen run the following in terminal then open R -- export R_MAX_VSIZE=32000000000
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
#### stats  ######
library(lme4)
library(lmerTest)
#library(MASS)
library(car)
###  shan is my calc, number_of and shannon is from obis
sdat<- stat_rare[["read_low_sp"]]  ## names(sdat)   names(stat_raw)   names(stat_rare)
###scale scale=TRUE divides by the standanrd deviation default does both center and scale
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sdat<-sdat[sdat$filter_size_cat==".5-5",]

sdat$sTemp.<-scale(sdat$Temp.)
sdat$sDepth<-scale(sdat$Depth)
sdat$sLand_Dist<-scale(sdat$land_dist)
sdat$sLatitude<-scale(sdat$Latitude)
sdat$sChlorophyll<-scale(sdat$Chlorophyll)
sdat$sLatitude_abs<-scale(sdat$Latitude_abs)
sdat$sPresent.Surface.Temperature.Mean<-scale(sdat$Present.Surface.Temperature.Mean)
sdat$sPresent.Surface.Temperature.Range<-scale(sdat$Present.Surface.Temperature.Range)
sdat$sPresent.Surface.Primary.productivity.Mean<-scale(sdat$Present.Surface.Primary.productivity.Mean)
sdat$sPresent.Surface.Primary.productivity.Range<-scale(sdat$Present.Surface.Primary.productivity.Range)
sdat$sPop_in_100km.Mean<-scale(sdat$Pop_in_100km)
sdat$sPop_in_1000km<-scale(sdat$Pop_in_1000km)
sdat$sOHI_2013<-scale(sdat$OHI_2013)
### plot    untransformed similar  slope=.0000136, signigicant
plot(log(sdat$rich)~sqrt(sdat$reads))
abline(lm(log(sdat$rich)~sqrt(sdat$reads)))  ## summary(lm(sdat$rich~sdat$reads))


#######################################################################################
###   mes around
x<-sdat$reads   # sqrt best
x<-sdat$rich   # log best

hist(x)
hist(sqrt(x))
hist(log(x))

############################################################################################################
########################################################################################################
####   spatial stats  names(sdat)
library(spdep)
#Create a k=4 nearest neighbor set
ssdat<-sdat
#  ssdat<-ssdat[sdat$filter_size_cat==".5-5",]


sp::coordinates(ssdat) <- ~Longitude+Latitude
#
us.nb4<-knearneigh(coordinates(ssdat), k=8)
us.nb4<-knn2nb(us.nb4)
us.wt4<-nb2listw(us.nb4, style="W")
#
########################################################################################################
########################################################################################################
########################################################################################################
###  richness    names(sdat)     names(sdat2) summary(sdat$filter_size_cat)

lmer.mod<-lmer(log(rich)~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
    sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
    sPresent.Surface.Primary.productivity.Range+sOHI_2013+(1|Locality.2/TStation/Depth),data=sdat)##
#  summary(lmer.mod)
#   mes<- c(AIC(lmer.mod),AIC(lmer.mod2))
mod<-lm(log(rich)~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
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
AICs  # 101.3
###################  test interactions
#         for richnes AIC(fit.errdurb) 80     #   summary(fit.errdurb, Nagelkerke=T)
##check assumptions
mod<-fit.errdurb
hist(resid(mod))
plot(fitted(mod),resid(mod))
abline(0,0)
qqnorm(resid(mod))
qqline(resid(mod))

### check if interactions improve model
mod_in<-lm(log(rich)~(sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
                        sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
                        sPresent.Surface.Primary.productivity.Range+sOHI_2013)^2,sdat) # summary(mod_in)

mod_in<-lm(log(rich)~sDepth+ sLand_Dist + sPresent.Surface.Temperature.Mean+
             sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
             sPresent.Surface.Primary.productivity.Range+ sOHI_2013

         +  sPresent.Surface.Primary.productivity.Mean:sOHI_2013 
         +  sLand_Dist:sPresent.Surface.Primary.productivity.Range 
           
      
           ,sdat)  #sPresent.Surface.Temperature.Range:sDepth+

fit.errdurb2<-errorsarlm(mod_in, ssdat, listw=us.wt4, etype="emixed", method="MC")
summary(fit.errdurb2, Nagelkerke=T)



## plot sst
mod_in<-lm(log(rich)~sPresent.Surface.Temperature.Range
           ,sdat)  #sPresent.Surface.Temperature.Range:sDepth+
fit.errdurb2<-errorsarlm(mod_in, ssdat, listw=us.wt4, etype="emixed", method="MC")
summary(fit.errdurb2, Nagelkerke=T)
predict <- cbind(sdat, predict(fit.errdurb2, pred.type="TS"))

p <- ggplot(predict, aes(sPresent.Surface.Temperature.Range,signal)) +
  geom_point()
p <- p + geom_line(aes(sPresent.Surface.Temperature.Range, fit))
p <- p + geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3)
p

require(ggiraph)
require(ggiraphExtra)
library(spatialreg)
ggiraphExtra::ggPredict(fit.errdurb2,se=TRUE,interactive=TRUE,digits=3, list=us.wt4)

## plot?  length(ssdat$sPresent.Surface.Temperature.Range)
preds <- as.numeric(predict(fit.errdurb2)) # , type="response"
resid.space.and.ndvi <- residuals(fit.errdurb2)
plot(preds~ssdat$sPresent.Surface.Temperature.Range[2:72])
abline(lm(preds~ssdat$sPresent.Surface.Temperature.Range[2:72]))

##### plot interactions
mod_in<-lmer(log(rich)~filter_size_cat*
               sDepth+(1|Locality.2/TStation/Depth),data=sdat)
e <- effects::allEffects(mod_in)
plot(e)
## or test if real works
mod_in<-lm(log(rich)~filter_size_cat:sDepth ,ssdat) 
mod_in2<-errorsarlm(mod_in, ssdat, listw=us.wt4, etype="emixed", method="MC")
e <- effects::allEffects(mod_in2)
plot(e)

AIC(fit.sma2)   ## 80
########################################################################################################
########################################################################################################
############################################################################################################
##### plot interactions

#  sPresent.Surface.Primary.productivity.Mean:sOHI_2013 
#  sLand_Dist:sPresent.Surface.Primary.productivity.Range 

mod_in<-lmer(log(rich)~sPresent.Surface.Primary.productivity.Mean*sOHI_2013+(1|Locality.2/TStation/Depth),data=sdat)
e <- effects::allEffects(mod_in)
plot(e)

mod_in<-lmer(log(rich)~sPresent.Surface.Primary.productivity.Range*sLand_Dist+(1|Locality.2/TStation/Depth),data=sdat)
e <- effects::allEffects(mod_in)
plot(e)


############################################################################################################
########################################################################################################
############################################################################################################
########################################################################################################
## path analaysis
#      source('https://openmx.ssri.psu.edu/software/getOpenMx.R')

library(semPlot)
library(lavaan)

library(OpenMx)
library(tidyverse)
library(knitr)
library(kableExtra)
library(GGally)
padat<-sdat %>% 
  mutate(logrich=log(rich))
  
# specify model     log(rich) or shan  # model with 4 most importnat variables based on z-values
model <-'logrich~sPresent.Surface.Temperature.Mean+
          sPresent.Surface.Primary.productivity.Mean+
            sPresent.Surface.Primary.productivity.Range+sLatitude_abs'
## only sig to make simple plot
model <-'logrich~sDepth+sLand_Dist+sPresent.Surface.Temperature.Mean+
             sPresent.Surface.Temperature.Range+sPresent.Surface.Primary.productivity.Mean+
            sPresent.Surface.Primary.productivity.Range+sOHI_2013+sLatitude_abs'# run model
fit <- cfa(model, data = padat)
# view results
summary(fit, fit.measures = TRUE, standardized=T,rsquare=T)
#  plot
semPlot::semPaths(fit,what="std",layout="circle")
#semPlot::semPaths(fit,"std",layout = 'tree', edge.label.cex=.9, curvePivot = TRUE)
#  check location and change lavesl
semPaths(fit,what="std",nodeLabels=letters[1:5],edgeLabels=1:10,edge.label.cex=1.5,fade=FALSE)
semPlot::semPaths(fit,what="std",layout="circle",nodeLabels=letters[1:5])

# best!!!
labb<-c("richness ","SST mean","PP mean ","PP range"," Latitude ")
labb<-c("richness ", "   Depth  ","Dist. land","SST mean","SST range","PP mean","PP range", "    OHI    "," Latitude ")

semPlot::semPaths(fit,what="std",layout="circle2",nodeLabels=labb
                  ,edge.label.cex=1.1,sizeMan=10,sizeMan2=6)


###  plot correlations of independent var.
GGally::ggcorr(padat[,c(40,41,43,45:48,51,42)], nbreaks = 6, label = T, low = "red3", high = "green3", 
       label_round = 2, name = "Correlation Scale", label_alpha = T, hjust = 0.75) +
  ggtitle(label = "Correlation Plot") +
  theme(plot.title = element_text(hjust = 0.6))





