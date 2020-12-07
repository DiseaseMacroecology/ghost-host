# ghost-host: loading data & calculating metrics
# Maxwell J. Farrell

require(ape)
require(picante)
require(caper)
require(dplyr)
require(tidyr)
require(ggplot2)


# Loading Source Data


###################################################
### Faruby et al. Island Extinction Data (2016) ###
###################################################

island <- read.table("../data/Faurby_Island_Endemicity.txt", as.is=T, header=TRUE)
island$binomial <- paste(island$Genus, island$Species, sep="_")

###############################################
### PHYLACINE 1.2 Mammal Phylogeny + Traits ###
###############################################

if(!file.exists("../data/faurby_tree_1.tre")){
  trees <- read.nexus("../data/Complete_tree_5831_Species.nex")
  tree <- trees[[1]]
  length(tree$tip.label)# 5831
  write.tree(tree, "../data/faurby_tree_1.tre")
} else tree <- read.tree("../data/faurby_tree_1.tre")

# Tree of extant species
extinct_sp <- sort(unique(island$binomial[island$Extinct_or_extant=="Extinct"]))
extant_tree <- drop.tip(tree,intersect(tree$tip.label,extinct_sp))

traits <- read.csv("../data/Trait_data.csv", as.is=TRUE)
tax <- traits[,c("Binomial.1.2", "Genus.1.2","Family.1.2","Order.1.2")] 
names(tax) <- c("binomial","genus","family","order")

not.matched <- sort(setdiff(tree$tip.label, tax$binomial))
not.matched.genus <- unique(sub("_.*","", not.matched))
# length(unique(not.matched.genus))# 15
# length(setdiff(not.matched.genus, tax$MSW05_Genus))# 15


######################################
### Host-Parasite Interaction Data ###
######################################

# From Farrell et al. 2020 (in review - figshare repo with data to be published)
hp <- read.csv("../data/hp_list_full.csv")
hp$host <- gsub(" ", "_", hp$host)

# length(sort(setdiff(hp$host, tree$tip.label)))
# missing 115 species in tree, including domestics...

# Subsetting domestic host parasite interactions
# Domestics
domestics <- c(
  "Bos_taurus",
  "Camelus_bactrianus",
  "Camelus_dromedarius",
  "Capra_hircus",
  "Equus_caballus",
  "Ovis_aries",
  "Sus_scrofa",
  "Canis_lupus",
  "Felis_silvestris",
  "Felis_catus",
  "Bubalus_bubalis",
  "Lama_glama",
  "Vicugna_vicugna",
  "Oryctolagus_cuniculus",
  # Others
  "Bison_bison",
  "Equus_asinus", #double check if this is african wild ass or normal donkey
  "Mus_musculus",
  "Rattus_rattus",
  "Rattus_norvegicus",
  "Rangifer_tarandus",
  "Cavia_porcellus")

dom <- hp[hp$host%in%domestics,]

# subsetting for now
hp_list <- hp[hp$host%in%tree$tip.label,]

# interactions lost because of these mismatches
# dim(hp)[1]-dim(hp_list)[1]
# round((dim(hp)[1]-dim(hp_list)[1])/dim(hp)[1],3)

### Checking synonymy table

syn <- read.csv("../data/Synonymy_table_with_unaccepted_species.csv", as.is=T)
missing <- sort(setdiff(hp$host, tree$tip.label))
# setdiff(syn$Binomial.1.2, tree$tip.label)
# length(setdiff(missing,syn$Binomial.1.2))
# intersect(missing,syn$Binomial.1.2)

syn$Elton_binomial <- paste0(syn$EltonTraits.1.0.Genus,"_",syn$EltonTraits.1.0.Species)
# length(intersect(missing,syn$Elton_binomial))/length(missing)
# Success! All missing species except are in Elton_binomial

syn_domestics <- syn[syn$Elton_binomial%in%domestics,]
# syn_domestics$Binomial.1.2[syn_domestics$Binomial.1.2%in%extinct_sp]
# Only Bos_primigenius is "Extinct"

## Setting "Bos_primigenius" to "Bos_taurus" for Phylacine data
## and changing from "Extinct" to "Extant" in island

syn$Binomial.1.2[syn$Binomial.1.2=="Bos_primigenius"] <- "Bos_taurus"
tree$tip.label[tree$tip.label=="Bos_primigenius"] <- "Bos_taurus"
tax$binomial[tax$binomial=="Bos_primigenius"] <- "Bos_taurus"
traits$Binomial.1.2[traits$Binomial.1.2=="Bos_primigenius"] <- "Bos_taurus"

# Tree of extant species
island$Extinct_or_extant[island$binomial=="Bos_primigenius"] <- "Extant"
island$binomial[island$binomial=="Bos_primigenius"] <- "Bos_taurus"

extinct_sp <- sort(unique(island$binomial[island$Extinct_or_extant=="Extinct"]))
extant_tree <- drop.tip(tree,intersect(tree$tip.label,extinct_sp))

# Swapping names
lookup <- setNames(syn$Binomial.1.2,syn$Elton_binomial)
missing <- sort(setdiff(hp$host, tree$tip.label))

# swap names in hp so that Phylacine data can remain consistent
hp$host[hp$host%in%missing] <- lookup[hp$host[hp$host%in%missing]]

# subsetting
hp_list <- hp[hp$host%in%tree$tip.label,]

# interactions lost because of these mismatches
# dim(hp)[1]-dim(hp_list)[1]# 79
# round((dim(hp)[1]-dim(hp_list)[1])/dim(hp)[1],3)# 0.002

# sort(setdiff(hp$host, tree$tip.label))
# only not accepted species names... should be OK to remove

# Removing database indicator from hp_list
hp_list <- unique(hp_list[,c("host","parasite","ParType")])



