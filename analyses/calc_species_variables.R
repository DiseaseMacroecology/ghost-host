# ghost-host: species level metrics
# Max Farrell

# load droptipTree
source("./sim_functions.R")


# Calculate per-species extinction metrics

# Calculating metrics for all species in tree (not subsetting to hosts in hp_list)

if (!file.exists("sp_extinction_metrics.RDS")) {

    dat <- data.frame(binomial=island$binomial[island$Extinct_or_extant=="Extant"])
    
    # Evolutionary Distinctivneess
    
    ED <- evol.distinct(tree,type="equal.splits")
    names(ED)[2] <- "ED"
    
    ED_ext <- evol.distinct(extant_tree,type="equal.splits")
    names(ED_ext)[2] <- "ED_ext"
    
    dat <- left_join(ED_ext, ED)
    dat$ED_gain <- dat$ED_ext - dat$ED
    
    
    # MPD gain
    
    coph <- cophenetic(tree)
    
    # index living species in coph
    ind <- rownames(coph)%in%extant_tree$tip.label
    
    mpd_tot <- apply(coph, 1, mean)
    mpd_tot <- data.frame(Species=names(mpd_tot), mpd_tot=mpd_tot)
    
    mpd_ext <- apply(coph[ind,ind], 1, mean)
    mpd <- data.frame(Species=names(mpd_ext), mpd_ext=mpd_ext)
    
    mpd <- left_join(mpd, mpd_tot)
    mpd$mpd_change <- mpd$mpd_ext - mpd$mpd_tot
    hist(mpd$mpd_change)
    
    dat <- left_join(dat, mpd)
    
    # MNTD gain
    
    mntd_tot <- apply(coph, 1, function(x) min(x[x>0]))
    mntd_tot <- data.frame(Species=names(mntd_tot), mntd_tot=mntd_tot)
    
    mntd_ext <- apply(coph[ind,ind], 1, function(x) min(x[x>0]))
    mntd <- data.frame(Species=names(mntd_ext), mntd_ext=mntd_ext)
    
    mntd <- left_join(mntd, mntd_tot)
    mntd$mntd_change <- mntd$mntd_ext - mntd$mntd_tot
    range(mntd$mntd_change)
    hist(mntd$mntd_change)
    
    dat <- left_join(dat, mntd)
    
    saveRDS(dat,"sp_extinction_metrics.RDS")

} else { dat <- readRDS("sp_extinction_metrics.RDS")}



# Calculate host specificities per parasite

if (!file.exists("host_specificities.rds")) {

    # Calculating clade matrix once for whole tree
    
    clade_mat <- clade.matrix(tree)
    
    host_spec <- data.frame(parasite=sort(unique(hp_list$parasite)))
    
    # pd total
    pd_total <- pd.calc(clade_mat)
    
    host_pd <- sapply(host_spec$parasite, function(para) 
    
          as.numeric(pd.calc(clade_mat, tip.subset=hp_list$host[hp_list$parasite==    para], root.edge=FALSE))
    
          ) 
        
    coph <- cophenetic(tree)
    
    host_mpd <- sapply(host_spec$parasite, function(para) 
    
              mean(coph[rownames(coph)%in%hp_list$host[hp_list$parasite==para],
                        colnames(coph)%in%hp_list$host[hp_list$parasite==para]])
          ) 
        
    min_above_zero <- function(x){
    
      if (length(x)>1) {
        min(x[x>0])
      } else {0}
    
    }
        
    host_mntd <- sapply(host_spec$parasite, function(para) 
    
              min_above_zero(coph[rownames(coph)%in%hp_list$host[hp_list$parasite==   para], 
                colnames(coph)%in%hp_list$host[hp_list$parasite==para]])
          ) 
    
    
    host_maxD <- sapply(host_spec$parasite, function(para) 
    
              max(coph[rownames(coph)%in%hp_list$host[hp_list$parasite==para], 
                colnames(coph)%in%hp_list$host[hp_list$parasite==para]])
          ) 
    
    
    host_spec <- data.frame(parasite=sort(unique(hp_list$parasite)),
                            mpd=host_mpd, mntd=host_mntd, maxD=host_maxD)
    
    host_spec$creepjump <- with(host_spec, mntd/maxD)
    host_spec$creepjump[is.na(host_spec$creepjump)] <- 0
    
    saveRDS(host_spec, "host_specificities.rds")

} else { host_spec <- readRDS("host_specificities.rds")}




# Calculate mean parasite specificities per host species


if (!file.exists("mean_specificities_per_host.rds")) {

    dat_small <- dat[dat$Species%in%hp_list$host,]
    
    mean_mpd <- sapply(dat_small$Species, function(host)
    
        mean(host_spec$mpd [host_spec$parasite%in%hp_list$parasite[hp_list$host==host]])   
    
        )
        
    mean_mntd <- sapply(dat_small$Species, function(host)
    
        mean(host_spec$mntd [host_spec$parasite%in%hp_list$parasite[hp_list$host==host]])   
    
        )
    
    
    mean_maxD <- sapply(dat_small$Species, function(host)
    
        mean(host_spec$maxD [host_spec$parasite%in%hp_list$parasite[hp_list$host==host]])   
    
        )
    
    
    mean_creepjump <- sapply(dat_small$Species, function(host)
    
        mean(host_spec$creepjump [host_spec$parasite%in%hp_list$parasite[hp_list$host==host]])   
    
        )
    
    mean_spec <- data.frame(Species=dat_small$Species,
                                para_mpd=mean_mpd, para_mntd=mean_mntd, 
                                para_maxD=mean_maxD, para_creepjump=mean_creepjump)
            
    saveRDS(mean_spec, "mean_specificities_per_host.rds")

} else { mean_spec <- readRDS("mean_specificities_per_host.rds")}


# Richnesses


if (!file.exists("richness.rds")) {

    comm <- table(hp_list[,c(1:2)])
    
    para_richness <- rowSums(comm)
    host_richness <- colSums(comm)
    
    mean_HSR <- sapply(dat$Species, function(host)
        
          mean(host_richness[colnames(comm)%in%hp_list$parasite[hp_list$host==host]]    )   
    )

    n_SHP <- sapply(dat$Species, function(host)
        
          sum(host_richness[colnames(comm)%in%hp_list$parasite[hp_list$host==host]]==1)   
    )


    # Zyzomys_woodwardi has mean_HSR of 4 with 1 parasite... double check function works:
    # para_richness
    # mean(host_richness[colnames(comm)%in%hp_list$parasite[hp_list$host=="Zyzomys_woodwardi"]])
    # hp_list[hp_list$host=="Zyzomys_woodwardi",]
    # hp_list$parasite[hp_list$host=="Zyzomys_woodwardi"]
    # hp_list[hp_list$parasite==hp_list$parasite[hp_list$host=="Zyzomys_woodwardi"],]

    para_richness[para_richness==max(para_richness)]
    # Homo sapiens...

    HSR <- data.frame(Species=names(mean_HSR), para_meanHSR=mean_HSR)
    PSR <- data.frame(Species=names(para_richness), npara=para_richness)
    SHP <- data.frame(Species=names(n_SHP), nSHP=n_SHP)

    richness <- left_join(HSR,PSR)            
    richness <- left_join(richness,SHP)            

    saveRDS(richness, "richness.rds")

} else { richness <- readRDS("richness.rds")}


