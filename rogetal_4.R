#this script is dependent on rogetal_1.R

# load climate -----------------------------------------------------------------
load("meansBioClim2.RData")
load("P_olsen_df.RData")


# environmental niche ------------------------------------------------------------

library(raster)
sf_use_s2(F)

load("global.inland.RData")

cent <- global.inland$cent %>%
  st_as_sf()

global_inland_points <- as(cent, "Spatial")

proj4string(global_inland_points) <- st_crs(global.inland)$proj4string

GBIFpolygon_groups$genus <- GBIFpolygon_groups$specie %>%
  stringr::word(.,1)

GBIF.polygon.dual <- GBIFpolygon_groups[which(GBIFpolygon_groups$genus %in% GBIF.dual.genus.f),]
GBIF.polygon.AM <- GBIFpolygon_groups[which(GBIFpolygon_groups$genus %in% GBIF.AM.genus.f),]  
GBIF.polygon.EM <- GBIFpolygon_groups[which(GBIFpolygon_groups$genus %in% GBIF.EM.genus.f),]  

values_dual <- st_intersects(global.inland,GBIF.polygon.dual) %>%
  melt(values_dual)
values_EM <- st_intersects(global.inland,GBIF.polygon.EM) %>% 
  melt(values_EM)
values_AM <- st_intersects(global.inland,GBIF.polygon.AM) %>% 
  melt(values_AM)

clim_spec_dual <- matrix(nrow = max(values_dual$value), ncol = ncol(meansBioClim2))
clim_spec_AM <- matrix(nrow = max(values_AM$value), ncol = ncol(meansBioClim2))
clim_spec_EM <- matrix(nrow = max(values_EM$value), ncol = ncol(meansBioClim2))

for (i in 1:max(values_dual$value)){
  
  spec_range <- values_dual$L1[which(values_dual$value == i)]
  clim_spec_dual[i,] <- apply(meansBioClim2[spec_range,], 2, mean, na.rm = T) 
}

clim_spec_dual <- as.data.frame(clim_spec_dual)
clim_spec_dual$spec <- GBIF.polygon.dual$specie
clim_spec_dual$genus <- stringr::word(clim_spec_dual$spec,1)
clim_spec_dual$group <- "dual"

for (i in 1:max(values_EM$value)){
  
  spec_range <- values_EM$L1[which(values_EM$value == i)]
  clim_spec_EM[i,] <- apply(meansBioClim2[spec_range,], 2, mean, na.rm = T) 
}
clim_spec_EM <- as.data.frame(clim_spec_EM)
clim_spec_EM$spec <- GBIF.polygon.EM$specie
clim_spec_EM$genus <- stringr::word(clim_spec_EM$spec,1)
clim_spec_EM$group <- "EM"

for (i in 1:max(values_AM$value)){
  
  spec_range <- values_AM$L1[which(values_AM$value == i)]
  clim_spec_AM[i,] <- apply(meansBioClim2[spec_range,], 2, mean, na.rm = T) 
}

clim_spec_AM <- as.data.frame(clim_spec_AM)
clim_spec_AM$spec <- GBIF.polygon.AM$specie
clim_spec_AM$genus <- stringr::word(clim_spec_AM$spec,1)
clim_spec_AM$group <- "AM"

clim_spec_all <- rbind(clim_spec_dual,clim_spec_AM,clim_spec_EM)
clim_spec_f <- clim_spec_all[which(clim_spec_all$genus %in% rownames(genus_matrix)),]

clim_spec_f <- clim_spec_f[complete.cases(clim_spec_f),-22]

vars_choose <- paste0("V",c(1,2,3,7,12,21,23,24,25))

clim_spec_mean <- clim_spec_f %>%
  group_by(spec) %>%
  # Summarise and check for the number of distinct groups
  summarise(across(vars_choose, mean, na.rm = TRUE),
            distinct_groups = n_distinct(group),
            group = first(group)) %>%
  # Check if there are multiple distinct groups and set to 'dual' if so
  mutate(group = if_else(distinct_groups > 1, "dual", group)) 


clim_spec_upper_quartile <- clim_spec_f %>%
  group_by(spec) %>%
  summarise(across(vars_choose, ~ quantile(.x, probs = 0.25, na.rm = TRUE)),
            distinct_groups = n_distinct(group),
            group = first(group)) %>%
  # Check if there are multiple distinct groups and set to 'dual' if so
  mutate(group = if_else(distinct_groups > 1, "dual", group)) 

clim_spec_range <- clim_spec_f %>%
  group_by(spec) %>%
  summarise(across(vars_choose, ~ max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)),
            distinct_groups = n_distinct(group),
            group = first(group)) %>%
  mutate(group = if_else(distinct_groups > 1, "dual", group)) 


clim_spec_mean$genus <- stringr::word(clim_spec_mean$spec,1)
clim_spec_upper_quartile$genus <- stringr::word(clim_spec_upper_quartile$spec,1)
clim_spec_range$genus <- stringr::word(clim_spec_range$spec,1)


clim_genus_mean <- clim_spec_mean %>% 
  group_by(genus) %>%
  summarise(across(vars_choose,mean,na.rm = T),
            group = first(group))

clim_genus_upper_quartile <- clim_spec_upper_quartile %>%
  group_by(genus) %>%
  summarise(across(vars_choose, ~ quantile(.x, probs = 0.25, na.rm = TRUE)),
            group = first(group))

clim_genus_range <- clim_spec_range %>%
  group_by(genus) %>%
  summarise(across(vars_choose, ~ max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)),
            group = first(group))


clim_genus_mean_f <- clim_genus_mean[which(clim_genus_mean$genus %in% tree.new$tip.label),] %>%
  as.data.frame()

clim_genus_range_f <- clim_genus_range[which(clim_genus_range$genus %in% tree.new$tip.label),] %>%
  as.data.frame()

clim_AM <- clim_genus_mean$V1[which(clim_genus_mean$genus %in% GBIF.AM.genus.f)]
clim_dual <- clim_genus_mean$V1[which(clim_genus_mean$genus %in% GBIF.dual.genus.f)]
clim_EM <- clim_genus_mean$V1[which(clim_genus_mean$genus %in% GBIF.EM.genus.f)]

comp_data <- caper::comparative.data(phy = tree.new, data = clim_genus_range_f, names.col = "genus", vcv = TRUE)

a <- aov(V1 ~ group, data = comp_data$data)
summary(a)
TukeyHSD(a, "group")

clim2plot <- rbind(clim_genus_range,clim_genus_upper_quartile,clim_genus_mean)
clim2plot$test <- rep(c("range","lower_quartile","mean"), each = nrow(clim_genus_range))

clim_plot <- melt(clim_spec_mean)
library(ggpubr)
ggplot(clim_plot[which(clim_plot$variable!="distinct_groups"),], aes(x = group, y = value, group = group, fill = group)) + 
  geom_boxplot() + 
  facet_wrap(~ variable, scales = "free_y") + 
  stat_compare_means(
    aes(label = ..p.signif..), 
    method = "t.test", 
    comparisons = list(c("AM", "EM"), c("AM", "dual"), c("dual", "EM")),
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05),  # Define p-value thresholds
      symbols = c("***", "**", "*")     # Corresponding symbols
    )
  ) + 
  theme_minimal()

# Blomberg K  --------------------------------------------

traits_named <- setNames(comp_data$data$V12, rownames(comp_data$data))
k <- phytools::phylosig(comp_data$phy, traits_named,method="K")
print(k)
