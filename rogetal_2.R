library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(ggstar)

# create comparative data for niche conservatism -------------------------------------------------------------------------

genus_matrix <- matrix(0, nrow = length(genus.phylo.present.f), ncol = 3)
colnames(genus_matrix) <- c("EM", "AM", "dual")
rownames(genus_matrix) <- genus.phylo.present.f

# Step 2: Fill the matrix
for (i in 1:length(genus.phylo.present.f)) {
  genus_name <- genus.phylo.present.f[i]
  if (genus_name %in% GBIF.EM.genus.f) genus_matrix[i, "EM"] <- 1
  if (genus_name %in% GBIF.AM.genus.f) genus_matrix[i, "AM"] <- 1
  if (genus_name %in% GBIF.dual.genus.f) genus_matrix[i, "dual"] <- 1
}

genus.df <- as.data.frame(genus_matrix) %>%
  mutate(genus = rownames(.))

group_df <- data.frame(genus = rownames(genus_matrix), 
                       group = apply(genus_matrix, 1, function(x) {
                         ifelse(x[3] == 1, "dual", ifelse(x[2] == 1, "AM", "EM"))
                       }))


comp_data <- caper::comparative.data(phy = tree.new, data = group_df, names.col = "genus", vcv = TRUE)
genus.factor <- as.character(comp_data$data$group)

# plot tree --------------------------------------------------------------------


species_df <- comp_data$data[,-c(1,2)]
colnames(species_df) <- names_biome$names
species_df$genus <- rownames(species_df)

species_clusters <- melt(species_df, id = "genus", variable.name = "Cluster", value.name = "Presence")
species_clusters$Presence[grepl(1, species_clusters$Presence)] = "Present"
species_clusters$Presence <- paste(species_clusters$Cluster,species_clusters$Presence)
species_clusters$Presence[grepl("0", species_clusters$Presence)] = "Absent"
species_clusters$Cluster <- factor(species_clusters$Cluster, levels=unique(species_clusters$Cluster)[c(1,2,7,9,14,3,10,13,8,12,4,5,6,11)])
species_clusters$Presence <- factor(species_clusters$Presence,
                                    levels= c(paste(unique(species_clusters$Cluster)[c(1,2,7,9,14,3,10,13,8,12,4,5,6,11)],"Present"),"Absent"))

species_clusters$Presence <- factor(species_clusters$Presence,
                                    levels= c(paste(unique(species_clusters$Cluster),"Present"),"Absent"))


# Creating a single group label column for plotting
group_df <- data.frame(genus = rownames(genus_matrix), 
                       group = apply(genus_matrix, 1, function(x) {
                         ifelse(x[3] == 1, "dual", ifelse(x[2] == 1, "AM", "EM"))
                       }))



p <- ggtree(tr=comp_data$phy, layout="fan", open.angle=5, size=0.2)
p <- p %<+% group_df + 
  geom_tippoint(aes(color = group), size = 3)


#color_palette <- c(color_palette, "#3f3f3f")
color_palette <- c("#1B4414","cornflowerblue","darkgoldenrod1","blueviolet","white")

p2 <- p +
  geom_fruit(
    data=species_clusters,
    geom=geom_tile,
    mapping=aes(x=Cluster, y=genus, fill=Presence),
    #width=0.1,
    #color="white",
    pwidth=0.3,
    offset=0.1
  ) + 
  scale_fill_manual(
    values=color_palette,
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  ) +
  new_scale_fill()

p2




# mpd ---------------------------------------------------------------------

phydist <- cophenetic(comp_data$phy)
mpd_alpha <- picante::ses.mpd(t(comp_data$data),phydist, null.model = "taxa.labels",runs = 1000, abundance.weighted = F)

mpd_alpha$cluster <- rownames(mpd_alpha)
mpd_alpha$significance <- ifelse(mpd_alpha$mpd.obs.p < 0.075, "clustered", ifelse(mpd_alpha$mpd.obs.p>0.95,"dispersed", "NS"))

mpd.plot <- ggplot(mpd_alpha, aes(cluster,mpd.obs.z))+
  geom_point(aes(col = significance), size = 4) + 
  geom_hline(yintercept = -1.64)+
  geom_hline(yintercept = 1.64)+
  scale_color_manual(values = c("red","black","grey"))+
  #scale_x_continuous(labels = c(1:k), breaks = c(1:k)) +
  theme_classic()+
  theme(legend.position = "none")



# Delta ------------------------------------------------------------
library(ape)
source("code_delta.R")

delta <- function(trait, tree,lambda0,se,sim,thin,burn) {
  
  ar <- ace(trait,tree,type="discret",method="ML",model="ARD")$lik.anc
  
  # deletes the complex part whenever it pops up
  if (class(ar[1,1]) == "complex"){
    ar <- Re(ar)
  }
  
  x  <- nentropy(ar)
  mc1    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mc2    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mchain <- rbind(mc1,mc2)
  deltaA <- mean(mchain[,2]/mchain[,1])
  
  return(deltaA)
}

deltaA <- delta(genus.factor,comp_data$phy,0.1,0.0589,10000,10,100)

random_delta <- rep(NA,100)

for (i in 1:100){
  rtrait <- sample(genus.factor)
  random_delta[i] <- delta(rtrait,comp_data$phy,0.1,0.0589,1000,10,100)
}

random_delta <- random_delta[1:i-1]
p_value <- sum(random_delta>deltaA)/length(random_delta)
boxplot(random_delta)
abline(h=deltaA,col="red")

