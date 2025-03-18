load("names_biome.RData")
load("wwf_simple.b.RData")
load("GBIFpolygon_groups.RData")

#refer to tree.new from Rogetal_1.R
load(treenew) 

sf_use_s2(F)



boot_fun <- function(clus.globe,GBIF.all){
  
  require(caper)
  require(tidyr)
  require(phylolm)
  require(sf)
  
  GBIF.AM.genus.f <- GBIF.all[1] %>% unlist()
  GBIF.EM.genus.f <- GBIF.all[2] %>% unlist()
  GBIF.dual.genus.f <- GBIF.all[3] %>% unlist()
  
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
  
  
  
  
  biomes_all_wide <- as.data.frame(biomes_all_wide)
  biomes_all_wide$genus_id <- as.numeric(rownames(biomes_all_wide))
  
  comp_data <- comparative.data(phy = tree.new, data = biomes_all_wide, 
                                names.col = genus, vcv = TRUE, na.omit = TRUE)       
  
  
  comp_data$data$group <- as.factor(comp_data$data$group)
  
  
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
    
  }
  
  return(model_save)
  
}


# unconfirmed ratio  ------------------------------------------------------

boot_run <- list()

total_species <- length(GBIF.AM.genus.f) + length(GBIF.EM.genus.f)

sample_AM_size <- floor(length(GBIF.AM.genus.f2)*0.11)
sample_EM_size <- floor(length(GBIF.EM.genus.f2)*0.75)

for (i in 1:100) {
    
  sampled_AM <- sample(GBIF.AM.genus.f2, sample_AM_size, replace = FALSE)
  sampled_EM <- sample(GBIF.EM.genus.f2, sample_EM_size, replace = FALSE)
  
  GBIF.AM.genus.f_new <- setdiff(GBIF.AM.genus.f2, sampled_AM)
  GBIF.EM.genus.f_new <- setdiff(GBIF.EM.genus.f2, sampled_EM)
  GBIF.dual.genus.f_new <- c(GBIF.dual.genus.f,sampled_AM,sampled_EM)
  
  GBIF.all <- list(GBIF.AM.genus.f_new,GBIF.EM.genus.f_new,GBIF.dual.genus.f_new)
  
  boot_run[[i]] <- boot_fun(clus.globe,GBIF.all)  
}


# random ratio  ------------------------------------------------------------

boot_run <- list()

percent_FP <- c(0.05,0.1,0.2)
percent_FN <- c(0.05,0.1,0.2)

NDual_species <- c(GBIF.AM.genus.f,GBIF.EM.genus.f)
Dual_species <- (GBIF.dual.genus.f)

percent_all <- expand.grid(percent_FN,percent_FP)

for(j in 1:9){
  
  percent_FP <- percent_all[j,][1]
  percent_FN <- percent_all[j,][2]
  
  FP_number <- ceiling(percent_FP * length(NDual_species))
  FN_number <- ceiling(percent_FN * length(Dual_species))
  
  for (i in 1:100) {
    
    sampled_FP <- sample(NDual_species, FP_number, replace = FALSE)
    sampled_FN <- sample(Dual_species, FN_number, replace = FALSE)
    
    GBIF.AM.genus.f_new <- c(setdiff(GBIF.AM.genus.f, sampled_FP))
    GBIF.EM.genus.f_new <- setdiff(GBIF.EM.genus.f, sampled_spec)
    GBIF.dual.genus.f_new <- c(GBIF.dual.genus.f,sampled_spec)
    
    GBIF.all <- list(GBIF.AM.genus.f_new,GBIF.EM.genus.f_new,GBIF.dual.genus.f_new)
    
    boot_run[[i]] <- boot_fun(clus.globe,GBIF.all)  
  }
  
}


boot_run <- list(boot_run_positive_0.5,boot_run_positive_0.1,boot_run_positive_0.2)


# plot --------------------------------------------------------------------

boot_run <- list(boot_run_ratio_unconfirmed2,boot_run_positive_0.05)
boot_perc <- c("ratio unconfirmed", "original ratio")

#load the unconfirmed and confirmed summary results from Teste et al. using the 
#script in rogetal_3_clean.R

for (cols in cols_reg){
  
  model_summary_unconfirmed <- model_unconfirmed[[cols]]
  model_summary_confirmed <- model_confirmed[[cols]]
  
  combined_summary <- data.frame()
  for (boot in 1:2){
    boot_choose <- boot_run[[boot]]
    
    estimates_df <- do.call(rbind, lapply(1:100, function(n) {
      # Extract the data frame for each iteration and add an iteration column
      data <- boot_choose[[n]][[cols]]
      data$Iteration <- n
      return(data)
    }))
    
    summary_stats <- estimates_df %>%
      group_by(Row = rep(1:3, times = 100)) %>%  # Identify rows (1, 2, 3) across all iterations
      summarize(
        Estimate = mean(Estimate),
        lowerCI = mean(lowerCI),
        upperCI = mean(upperCI),
        CI_Lower = quantile(Estimate, 0.025),
        CI_Upper = quantile(Estimate, 0.975),
      )
    
    summary_stats$term <- model_summary$term
    summary_stats$source <- paste("Simulation",boot_perc[boot])
    
    combined_summary <- rbind(combined_summary,summary_stats[,c(2,3,4,7,8)])
    
  }
  
  model_summary_confirmed$source <- "Original Data"
  model_summary_unconfirmed$source <- "With Unconfirmed"
  
  combined_summary <- rbind(model_summary_confirmed[,c(1,5,6,7,8)],
                            model_summary_unconfirmed[,c(1,5,6,7,8)],
                            combined_summary)
  
  combined_summary$source <- factor(combined_summary$source, levels = rev(unique(combined_summary$source)))
  
  plot_n <- 
    ggplot(combined_summary, aes(x = Estimate, xmin = lowerCI, xmax = upperCI, y = term, color = source)) +
    geom_pointrange(position = position_dodge(width = 0.5)) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(x = "Effect Size (log odds)", y = "Comparison Group",
         title = cols,
         subtitle = "Reference Group: Dual",
         color = "Source")  # Legend for source
  
  plots[[cols]] <- plot_n
  
  
}


grid_plots <- (plots[[2]] | plots[[5]]) /
  (plots[[3]] | plots[[4]])

# Combine the first plot with the grid of remaining plots
plot_layout <- plots[[1]] + 
  plot_layout(guides = "collect") +
  grid_plots
