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


## load packages
packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','dismo',
                       'doParallel','abind','stringr','rgdal','foreign')
lapply(packs.to.extract,require, character.only=T)



# list.families=tolower(f2$Family)
#library("Hmisc")
#list.families=capitalize(list.families)
#write.csv(list.families,file='~/GitHub/host-ghost/data/list.families.csv')
# list.families=read.csv('~/GitHub/host-ghost/data/list.families.csv')


## load stack raster with all sps distributions
#setwd("~/Data_Harvard/PHYLACINE/Current/")
sps.list <- dir("../data/Ranges/Current/",pattern="*.tif",recursive=T, full.names=T)
mammal.raster=stack(sps.list)

names.raster <- names(mammal.raster)
sps.csv <- read.csv("../data/Synonymy_table_valid_species_only.csv")
valid.sps <- subset(sps.csv,Order.1.2%in%c("Carnivora","Cetartiodactyla","Perissodactyla","Primates"))

to.load<-which(names.raster%in%valid.sps$Binomial.1.2)
mammal.raster<-subset(mammal.raster,to.load)

ff<-as.matrix(array(1,dim=c(nlayers(mammal.raster),nlayers(mammal.raster))))
colnames(ff)=names(mammal.raster)
rownames(ff)=names(mammal.raster)
ff[upper.tri(ff)]=0
gg<-which(ff==1,arr.ind = T)
edgelist<-t(apply(gg,1,function(x){return(c(rownames(ff)[x[1]],colnames(ff)[x[2]]))}))


resultss<-apply(edgelist[1:10,],1,function(x){
  
  sps.i<-subset(mammal.raster,x[1])
  sps.j<-subset(mammal.raster,x[2])
  sumsps<-sps.i+sps.j
  overlap.mamm.range=ifelse(length(which(values(sumsps)==2))>0,
                                     length(which(values(sumsps)==2))/sum(values(sps.i)),0)
    
    return(overlap.mamm.range)
  
})


a<-c("a1","a2")
b<-c("b1","b2")
ff=list(a,b)
rast.sub<-subset(mammal.raster,1:4)
names(rast.sub)=c(a,b)
lapply(ff,FUN = function(x){return(subset(rast.sub,x[1])+subset(rast.sub,x[2]))})


## function to compute overlap
overlap.mamm.range=array(NA, dim=c(nlayers(mammal.raster),nlayers(mammal.raster)+1))
row.names(overlap.mamm.range)=names(mammal.raster)
for(i in 1:nlayers(mammal.raster)){
  
  layer.i=mammal.raster[[i]]
  overlap.mamm.range[i,1]=sum(values(layer.i))

  sps.i<-values(layer.i)
  for(j in i:nrow(overlap.mamm.range)){
    print(paste(i,j)) #j=1
    sps.j<-values(mammal.raster[[j]])
    compar.i.j<-rowSums(cbind(sps.i,sps.j))
    overlap.mamm.range[j,i+1]=ifelse(length(which(compar.i.j==2))>0,
           length(which(compar.i.j==2))/sum(values(layer.i)),0)
        
  }
  
}

# Takes a long time - maybe re-do this just for key clades?

## save results
## 
write.csv(overlap.mamm.range,file="../data/overlap.mamm.range.csv")

