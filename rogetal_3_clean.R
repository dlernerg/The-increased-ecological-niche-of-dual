#this script is dependent on rogetal_1.R

sf_use_s2(FALSE)

st_read("wwf_simple.b.geojson")
clus.globe <- wwf_simpl.b %>% 
  rename(cells = 1) %>%
  mutate(cluster = cells)

trop_geom <- clus.globe %>%
  slice(1, 2, 3, 7, 9, 10, 14) %>%
  st_combine() %>%
  as.data.frame() %>%
  mutate(biome = "trop")

temp_geom <- clus.globe %>%
  slice(4, 5, 8, 12) %>%
  st_combine() %>%
  as.data.frame() %>%
  mutate(biome = "temp")

desert_geom <- clus.globe %>%
  slice(13) %>%
  st_combine() %>%
  as.data.frame() %>%
  mutate(biome = "desert")


bor_geom <- clus.globe %>%
  slice(6,11) %>%
  st_combine() %>%
  as.data.frame() %>%
  mutate(biome = "bor")


clus.globe_sum <- rbind(trop_geom,temp_geom,desert_geom,bor_geom) %>% st_as_sf()

GBIF.polygon.dual <- GBIFpolygon_groups[which(GBIFpolygon_groups$genus %in% GBIF.dual.genus.f),]
GBIF.polygon.AM <- GBIFpolygon_groups[which(GBIFpolygon_groups$genus %in% GBIF.AM.genus.f),]  
GBIF.polygon.EM <- GBIFpolygon_groups[which(GBIFpolygon_groups$genus %in% GBIF.EM.genus.f),]  

GBIF.polygon.dual_comb <- GBIF.polygon.dual %>%
  group_by(genus) %>%
  summarise(geometry = st_combine(polygon))

GBIF.polygon.AM_comb <- GBIF.polygon.AM %>%
  group_by(genus) %>%
  summarise(geometry = st_combine(polygon))

GBIF.polygon.EM_comb <- GBIF.polygon.EM %>%
  group_by(genus) %>%
  summarise(geometry = st_combine(polygon))

biomes_dual <- st_intersects(clus.globe,GBIF.polygon.dual_comb) %>%
  melt()

biomes_AM <- st_intersects(clus.globe,GBIF.polygon.AM_comb) %>% 
  melt()

biomes_EM <- st_intersects(clus.globe,GBIF.polygon.EM_comb) %>% 
  melt()   

colnames(biomes_AM) <- colnames(biomes_dual) <- colnames(biomes_EM) <- c("genus_id","biomes")
biomes_AM$genus <- GBIF.polygon.AM_comb$genus[biomes_AM$genus_id]
biomes_AM$group <- "AM"
biomes_EM$genus <- GBIF.polygon.EM_comb$genus[biomes_EM$genus_id]
biomes_EM$group <- "EM"
biomes_dual$genus <- GBIF.polygon.dual_comb$genus[biomes_dual$genus_id]
biomes_dual$group <- "dual"

library(caper)
library(tidyr)

dual_genus <- unique(biomes_dual$genus)
biomes_AM_filtered <- biomes_AM %>% filter(!genus %in% dual_genus)
biomes_EM_filtered <- biomes_EM %>% filter(!genus %in% dual_genus)

biomes_all <- rbind(
  biomes_AM_filtered,
  biomes_EM_filtered,
  biomes_dual
)

biomes_all$biomes <- paste0("biome",biomes_all$biomes)

biomes_all_wide <- biomes_all %>%
  mutate(presence = 1) %>%  # Create a column to indicate presence
  pivot_wider(
    names_from = biomes,      # Create new columns from the biome values
    values_from = presence,   # Fill new columns with values from 'presence'
    values_fill = list(presence = 0)  # Fill absent biomes with 0
  )


#contingency_table <- table(comp_data$data$group, comp_data$data$biome1)

#####
library(phylolm)
read.csv("names_biome.csv")

biomes_all_wide <- as.data.frame(biomes_all_wide)
biomes_all_wide$genus_id <- as.numeric(rownames(biomes_all_wide))

comp_data <- comparative.data(phy = tree.new, data = biomes_all_wide, 
                              names.col = genus, vcv = TRUE, na.omit = TRUE)       


comp_data$data$group <- as.factor(comp_data$data$group)
comp_data$data$group <- relevel(comp_data$data$group, ref = "AM")  # Change reference group to "EM"


# with dummy vars ---------------------------------------------------------

comp_data$data2 <- comp_data$data %>%
  mutate(group_AM = ifelse(group == "AM", 1, 0),
         group_EM = ifelse(group == "EM", 1, 0),
         group_dual = ifelse(group == "dual", 1, 0))

group_myco <- unique(comp_data$data2$group)
biomes <- names(comp_data$data)[grepl("biome", names(comp_data$data))]  # adjust this as needed

beta_df <- data.frame()
p_val_df <- data.frame()

for (biome in biomes){
  for (group in group_myco){
    
    formula <- paste(biome, "~", paste0("group_",group))
    
    model <- phyloglm(formula, data = comp_data$data2, phy = comp_data$phy)
    summary(model)
    model_normal <- glm(formula, family = binomial(), data = comp_data$data2)
    
    beta <- c(coef(summary(model))[2,1],coef(summary(model_normal))[2,1],biome,group)
    p_val <- c(coef(summary(model))[2,4],coef(summary(model_normal))[2,4],biome,group)
    
    beta_df <- rbind(beta_df, beta)
    p_val_df <- rbind(p_val_df,p_val)
  }
  
}

colnames(beta_df) <- colnames(p_val_df) <- c("pgls","normal","biome","group")

beta_df_long <- beta_df %>% 
  pivot_longer(cols = c(pgls, normal), names_to = "type", values_to = "values")
p_val_df_long <- p_val_df %>% 
  pivot_longer(cols = c(pgls, normal), names_to = "type", values_to = "values_pval")

beta_df_long <- beta_df_long %>%
  mutate(biome_num = as.numeric(str_extract(biome, "\\d+"))) %>%
  mutate(biome_name = names_biome$names[match(biome_num, names_biome$V2)])

beta_df_long$biome_name <- factor(beta_df_long$biome_name, levels = names_biome$names[c(1,2,7,9,14,3,10,13,8,12,4,5,6,11)])


beta_df_long$values <- as.numeric(beta_df_long$values)
p_val_df_long$values_pval <- as.numeric(p_val_df_long$values_pval)

beta_df_long <- left_join(beta_df_long,p_val_df_long, by = c("biome", "group", "type"))
beta_df_long$vjust <- ifelse(beta_df_long$values > 0, -0.3, 0.9)
beta_df_long$type <- as.factor(beta_df_long$type)

ggplot(beta_df_long, aes(x = group, y = values, fill = group, group = type, alpha = type)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ biome_name) + 
  geom_text(aes(label = ifelse(values_pval < 0.06, "*", "")),
            position = position_dodge(width = 0.8)) +  
  labs(y = "Beta Coefficient", title = "Beta Coefficients and Significance across Biomes and Models") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_fill_manual(values = c("normal" = "grey", "PGLS" = "red")) +
  scale_alpha_manual(values = c(0.5, 1)) + # Adjust these values as needed for different alpha levels
  expand_limits(y = c(NA, 1))  # Adjust this value to ensure all asterisks are visible

ggsave(filename = "glm_pgls_individual_allbiomes_normvspgls.svg", dpi = 600)


# poisson, ----------------------------------------------------------------


comp_data_sum <- comp_data$data %>%
  rowwise() %>%
  mutate(total_presence = sum(c_across(starts_with("biome"))),
         presence_trop = sum(c_across(c(biome1, biome2, biome3, biome7,biome9,biome14))),
         presence_temp = sum(c_across(c(biome4,biome5,biome8,biome10,biome12,))),
         presence_bor = sum(c_across(c(biome6,biome11))),
         presence_des = biome13,
  )%>%
  ungroup()


comp_data_sum <- comp_data_sum %>%
  dplyr::select(genus_id, group, total_presence,presence_trop,presence_temp,presence_bor,presence_des) %>%
  as.data.frame()

comp_data_sum$genus <- rownames(comp_data$data)[match(comp_data_sum$genus_id,comp_data$data$genus_id)] 
comp_data_sum$genus <- as.factor(comp_data_sum$genus)

comp_data_new <- comparative.data(phy = tree.new, data = comp_data_sum, 
                                  names.col = genus, vcv = TRUE, na.omit = TRUE)       

comp_data_new$data$group <-  relevel(comp_data_new$data$group, ref = "dual")  

cols_reg <- colnames(comp_data_new$data)[3:7]

plots <- list()
model_save <- list()
for (cols in cols_reg){
  
  formula <- paste(cols, "~ group")
  if (cols == "presence_des"){
    phylo_poisson_model <- phyloglm(formula, data = comp_data_new$data, phy = comp_data_new$phy)
    
  } else{
    phylo_poisson_model <- phyloglm(formula, data = comp_data_new$data, phy = comp_data_new$phy, method = "poisson_GEE")
    
  }
  model_summary <- summary(phylo_poisson_model)$coefficients %>% as.data.frame()
  model_summary$lowerCI <- model_summary$Estimate - 1.96 * model_summary$StdErr
  model_summary$upperCI <- model_summary$Estimate + 1.96 * model_summary$StdErr
  model_summary$term <- rownames(model_summary)
  model_summary$term <- c("Intercept (dual)","AM","EM")
  model_summary$term <- factor(model_summary$term, levels = model_summary$term[3:1])
  model_save[[cols]] <- model_summary
  plot_n <- 
    ggplot(model_summary, aes(x = Estimate, xmin = lowerCI, xmax = upperCI, y = term)) +
    geom_pointrange() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(x = "Effect Size (log odds)", y = "Comparison Group",
         title = cols,
         subtitle = "Reference Group: Dual") 
  
  plots[[cols]] <- plot_n
  
}

library(patchwork)

grid_plots <- (plots[[2]] | plots[[5]]) /
  (plots[[3]] | plots[[4]])

# Combine the first plot with the grid of remaining plots
plot_layout <- plots[[1]] + 
  plot_layout(guides = "collect") +
  grid_plots


ggsave(filename = "poisson_pgls_all.svg", width = 10, height = 6)

