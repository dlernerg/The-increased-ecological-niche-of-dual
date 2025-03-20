library(readxl)
library(dplyr)
library(sf)
library(ggplot2)
library(reshape2)
library(phytools)
library(ape)
library(stringr)


GBIFpolygon_groups <- st_read("GBIFpolygon_groups.geojson")
confirmed.dual <- read_csv("confirmed.cvs")

colnames(confirmed.dual) <- confirmed.dual[2,]
confirmed.dual <- confirmed.dual[3:nrow(confirmed.dual),]

GBIF.data <- GBIFpolygon_groups$specie %>% as.data.frame()
GBIF.data$genus <- stringr::word(GBIF.data[,1],1)

GBIFpolygon_groups$genus <- GBIFpolygon_groups$specie %>%
  stringr::word(.,1) 


# species and genera from lists -------------------------------------------

genus.dual <- confirmed.dual$Genus

confirmed.EM <- read_csv("EM_confirmed.csv")

genus.EM <- confirmed.EM$Genus

confirmed.AM <- read_csv("AM_confirmed.csv")

genus.AM <- confirmed.AM$Genus


# list in the phylogeny -------------------------------------------------

tree_Sanchez_etal<-read.tree(file = "SanchezM (2020).tre")

tree <- tree_Sanchez_etal

genus.phylo <- tree$tip.label

genus.phylo.absent <- genus.phylo[which(!(genus.phylo %in% c(genus.EM,genus.AM,genus.dual)))]
genus.phylo.present <- genus.phylo[which((genus.phylo %in% c(genus.EM,genus.AM,genus.dual)))]

tree.new <- drop.tip(tree,genus.phylo.absent)

length(tree.new$tip.label) / length(tree$tip.label)

# which list in GBIF and phylo------------------------------------------------------

GBIF.dual.genus <- unique(GBIF.data[which(GBIF.data$genus %in% genus.dual),1]) %>%
  stringr::word(.,1) %>%
  unique()
  
GBIF.EM.genus <- unique(GBIF.data[which(GBIF.data$genus %in% genus.EM),1]) %>%
  stringr::word(.,1) %>%
  unique()

GBIF.AM.genus <- unique(GBIF.data[which(GBIF.data$genus %in% genus.AM),1]) %>%
  stringr::word(.,1) %>%
  unique()
#list and in phylogeny present in GBIF polygons

GBIF.dual.genus.f <- unique(GBIF.data[which(GBIF.data$genus %in% genus.dual.f),1]) %>%
  stringr::word(.,1) %>%
  unique()

GBIF.EM.genus.f <- unique(GBIF.data[which(GBIF.data$genus %in% genus.EM),1]) %>%
  stringr::word(.,1) %>%
  unique()

GBIF.AM.genus.f <- unique(GBIF.data[which(GBIF.data$genus %in% genus.AM),1]) %>%
  stringr::word(.,1) %>%
  unique()

length(unique(c(GBIF.EM.genus.f,GBIF.AM.genus.f,GBIF.dual.genus.f)))/length(unique(c(genus.EM.f,genus.AM.f,genus.dual.f)))


#update the genus in the phylogeny with the genera present in the GBIF. 

genus.phylo.absent.f <- genus.phylo[which(!(genus.phylo %in% c(GBIF.AM.genus.f,GBIF.EM.genus.f,GBIF.dual.genus.f)))]
genus.phylo.present.f <- genus.phylo[which((genus.phylo %in% c(GBIF.AM.genus.f,GBIF.EM.genus.f,GBIF.dual.genus.f)))]

tree.new <- drop.tip(tree,genus.phylo.absent.f)

length(tree.new$tip.label) / length(tree$tip.label)



