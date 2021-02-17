##############################################################################################################
# Script and functions to run:
#' * get range overlap (continuous and binary) from PHYLACINE -- HOST GHOST PROJECT 
#'
#'  Ignacio Morales-Castilla, et al.
#'  started October 2018
##############################################################################################################



## to start
rm(list=ls())
options(stringsAsFactors=FALSE)


source("./load_data.R")

require(magrittr)

d2 <- island %<>% mutate(lbin=binomial)
f1 <- d2 %>% filter(Extinct_or_extant=="Extant") %>% group_by(Family) %>% summarize(nExtant=n()) ## 151 families with extant taxa
f2 <- d2 %>% filter(Extinct_or_extant=="Extinct") %>% group_by(Family) %>% summarize(nExtinct=n()) ## 78 familes with extinct taxa
f <- full_join(f1,f2) %>% replace_na(list(nExtant=0,nExtinct=0))
f$prop_ext <- f$nExtinct/(f$nExtinct+f$nExtant)
f %<>% arrange(prop_ext)

View(f)


## load packages
packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','dismo',
                       'doParallel','abind','stringr','rgdal','foreign')
lapply(packs.to.extract,require, character.only=T)



# #list.families=tolower(f2$Family)
# #library("Hmisc")
# #list.families=capitalize(list.families)
# #write.csv(list.families,file='~/GitHub/host-ghost/data/list.families.csv')
# list.families=read.csv('~/GitHub/host-ghost/data/list.families.csv')



## load stack raster with only elephants sps distributions
sps.list <- dir("../data/Ranges/Current/",pattern="*.tif",recursive=T, full.names=T)
sps.list.past<-dir("../data/Ranges/Present_natural/",pattern="*.tif",recursive=T, full.names=T)

sps.csv <- read.csv("../data/Synonymy_table_valid_species_only.csv")

# Carnivora
carn.sps <- subset(sps.csv,Order.1.2=="Carnivora")

# current
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Current/"),"\\."),
                 function(x){return(x[1])}))
pos.carn <- which(namesdir%in%carn.sps$Binomial.1.2)
carn.raster.curr=stack(sps.list[pos.carn])

# past
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Present_natural/"),"\\."),
                          function(x){return(x[1])}))
pos.carn <- which(namesdir%in%carn.sps$Binomial.1.2)
carn.raster.past=stack(sps.list.past[pos.carn])

rich.carn.curr<-calc(carn.raster.curr,sum)
rich.carn.past<-calc(carn.raster.past,sum)

dev.off()
par(mfrow=c(2,1),mar=c(1,1,1,1))
plot(rich.carn.curr)
plot(rich.carn.past)

## save maps for present and past
writeRaster(rich.carn.curr,filename = "../plots_tables/maps/carn.curr.tif",format="GTiff")
writeRaster(rich.carn.past,filename = "../plots_tables/maps/carn.wextinct.tif",format="GTiff")



# Felidae
felid.sps <- subset(sps.csv,Family.1.2=="Felidae")

# current
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Current/"),"\\."),
                 function(x){return(x[1])}))
pos.felid <- which(namesdir%in%felid.sps$Binomial.1.2)
felid.raster.curr=stack(sps.list[pos.felid])

# past
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Present_natural/"),"\\."),
                          function(x){return(x[1])}))
pos.felid <- which(namesdir%in%felid.sps$Binomial.1.2)
felid.raster.past=stack(sps.list.past[pos.felid])

rich.felid.curr<-calc(felid.raster.curr,sum)
rich.felid.past<-calc(felid.raster.past,sum)

dev.off()
par(mfrow=c(2,1),mar=c(1,1,1,1))
plot(rich.felid.curr)
plot(rich.felid.past)

## save maps for present and past
writeRaster(rich.felid.curr,filename = "../plots_tables/maps/felid.curr.tif",format="GTiff")
writeRaster(rich.felid.past,filename = "../plots_tables/maps/felid.wextinct.tif",format="GTiff")


# Rhinocerotidae
rhinos.sps <- subset(sps.csv,Family.1.2=="Rhinocerotidae")

# current
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Current/"),"\\."),
                 function(x){return(x[1])}))
pos.rhinos <- which(namesdir%in%rhinos.sps$Binomial.1.2)
rhinos.raster.curr=stack(sps.list[pos.rhinos])

# past
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Present_natural/"),"\\."),
                          function(x){return(x[1])}))
pos.rhinos <- which(namesdir%in%rhinos.sps$Binomial.1.2)
rhinos.raster.past=stack(sps.list.past[pos.rhinos])

rich.rhinos.curr<-calc(rhinos.raster.curr,sum)
rich.rhinos.past<-calc(rhinos.raster.past,sum)

dev.off()
par(mfrow=c(2,1),mar=c(1,1,1,1))
plot(rich.rhinos.curr)
plot(rich.rhinos.past)

## save maps for present and past
writeRaster(rich.rhinos.curr,filename = "../plots_tables/maps/rhinos.curr.tif",format="GTiff")
writeRaster(rich.rhinos.past,filename = "../plots_tables/maps/rhinos.wextinct.tif",format="GTiff")


# Bovidae
bovids.sps <- subset(sps.csv,Family.1.2=="Bovidae")

# current
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Current/"),"\\."),
                 function(x){return(x[1])}))
pos.bovids <- which(namesdir%in%bovids.sps$Binomial.1.2)
bovids.raster.curr=stack(sps.list[pos.bovids])

# past
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Present_natural/"),"\\."),
                          function(x){return(x[1])}))
pos.bovids <- which(namesdir%in%bovids.sps$Binomial.1.2)
bovids.raster.past=stack(sps.list.past[pos.bovids])

rich.bovids.curr<-calc(bovids.raster.curr,sum)
rich.bovids.past<-calc(bovids.raster.past,sum)

dev.off()
par(mfrow=c(2,1),mar=c(1,1,1,1))
plot(rich.bovids.curr)
plot(rich.bovids.past)

## save maps for present and past
writeRaster(rich.bovids.curr,filename = "../plots_tables/maps/bovids.curr.tif",format="GTiff")
writeRaster(rich.bovids.past,filename = "../plots_tables/maps/bovids.wextinct.tif",format="GTiff")


# Cervidae
cervids.sps <- subset(sps.csv,Family.1.2=="Cervidae")

# current
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Current/"),"\\."),
                 function(x){return(x[1])}))
pos.cervids <- which(namesdir%in%cervids.sps$Binomial.1.2)
cervids.raster.curr=stack(sps.list[pos.cervids])

# past
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Present_natural/"),"\\."),
                          function(x){return(x[1])}))
pos.cervids <- which(namesdir%in%cervids.sps$Binomial.1.2)
cervids.raster.past=stack(sps.list.past[pos.cervids])

rich.cervids.curr<-calc(cervids.raster.curr,sum)
rich.cervids.past<-calc(cervids.raster.past,sum)

dev.off()
par(mfrow=c(2,1),mar=c(1,1,1,1))
plot(rich.cervids.curr)
plot(rich.cervids.past)

## save maps for present and past
writeRaster(rich.cervids.curr,filename = "../plots_tables/maps/cervids.curr.tif",format="GTiff")
writeRaster(rich.cervids.past,filename = "../plots_tables/maps/cervids.wextinct.tif",format="GTiff")



# Ursidae
ursids.sps <- subset(sps.csv,Family.1.2=="Ursidae")

# current
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Current/"),"\\."),
                 function(x){return(x[1])}))
pos.ursids <- which(namesdir%in%ursids.sps$Binomial.1.2)
ursids.raster.curr=stack(sps.list[pos.ursids])

# past
namesdir <- unlist(lapply(strsplit(dir("../data/Ranges/Present_natural/"),"\\."),
                          function(x){return(x[1])}))
pos.ursids <- which(namesdir%in%ursids.sps$Binomial.1.2)
ursids.raster.past=stack(sps.list.past[pos.ursids])

rich.ursids.curr<-calc(ursids.raster.curr,sum)
rich.ursids.past<-calc(ursids.raster.past,sum)

dev.off()
par(mfrow=c(2,1),mar=c(1,1,1,1))
plot(rich.ursids.curr)
plot(rich.ursids.past)

## save maps for present and past
writeRaster(rich.ursids.curr,filename = "../plots_tables/maps/ursids.curr.tif",format="GTiff")
writeRaster(rich.ursids.past,filename = "../plots_tables/maps/ursids.wextinct.tif",format="GTiff")



dev.off()

plot(rich.ursids.past)

plot(rich.ursids.past,
     box = FALSE,
     axes = FALSE,
     col = grey(100:1/100),
     main = "grayscale ursid map")



# ## load stack raster with all current sps distributions
# #setwd("~/Data_Harvard/PHYLACINE/Current/")
# sps.list<-dir("~/Data_Harvard/PHYLACINE/Current/",pattern="*.tif",recursive=T, full.names=T)
# mammal.raster.curr=stack(sps.list)
# #4 minutess

# ## load stack raster with all current sps distributions
# sps.list<-dir("~/Data_Harvard/PHYLACINE/Present_natural/",pattern="*.tif",recursive=T, full.names=T)
# mammal.raster.past=stack(sps.list)


# rich.mamm.curr<-calc(mammal.raster.curr,sum)
# rich.mamm.past<-calc(mammal.raster.past,sum)


# ## save maps for present and past

# writeRaster(rich.carn.curr,filename = "~/yourpath/here/elephants.curr.tif",format="GTiff")
# writeRaster(rich.carn.past,filename = "~/yourpath/here/elephants.wextinct.tif",format="GTiff")

