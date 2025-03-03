# Scripts and Data for "The-increased-ecological-niche-of-dual mycorrhizal tree species"

## Description

This repository contains the R scripts and necessary data used for the analysis in the paper **The increased ecological niche of dual mycorrhizal tree species** by Rog et al. These scripts are designed to reproduce the analyses described in the paper.

### Scripts

1. **rogetal_1.R**:  
   This script performs the initial data cleaning and preprocessing. It reads in the `confirmed.xlsx`, `AM_confirmed.xlsx`, and `EM_confirmed.xlsx` files, combines the confirmed tree species with a phylogenetic tree from Sanchez Martinez et al., (2020) and outputs cleaned data that will be used in subsequent analysis.

2. **rogetal_2.R**:  
   This script conducts the phylogenetic analyses observed in figure 1 (phylogenetic map plotting, MPD and delta statistics), using the preprocessed data from `rogetal_1.R` file. 

3. **rogetal_3.R**:  
   This script runs the biome overlap data observed in Figure 2. It uses the cleaned data and `names_biome.RData` and `wwf_simple.b.RData` to generate binomial and poisson models of the occupided space for each of the mycorrhizal associating tree species groups.

4. **rogetal_4.R**:  
   This script performs the final analysis of the paper (Figure 3) - an analysis on the niche space. It combines the data generated in rogetal_1.R with environmental data from WorldClim data (Hijmans et al. 2015) and Olsen P data (McDowell et al. 2023).

5. **rogetal_5_sims2**
   This script performs the simulations of the paper. One needs to run rogetal_3.R before running the simulations, and load the confirmed and unconfirmed model summary dataframes as "model_unconfirmed" and "model_confirmed".
     
### Data Files

1. **confirmed.xlsx**:  
   This Excel file contains confirmed dual mycorrhizal genera from Teste et al., 2020.

2. **AM_confirmed.xlsx**:  
   This Excel file contains confirmed AM mycorrhizal genera from Soudzilovskaia et al., 2020.

3. **EM_confirmed.xlsx**:  
   This Excel file contains confirmed EM mycorrhizal genera from Soudzilovskaia et al., 2020.

4. **global.inland.RData**:  
   RData file containing grid cells of the world map using the 'dggridR' package which employs the ISEA Discrete Global Grid System for equal-sized grid projection on a 2-dimensional plane.

5. **wwf_simple.b.RData**:  
   Contains polygons for the wwf biomes, whereby the polygons were smoothened, as described in Lerner et al., 2023.

6. **names_biomes.RData**:  
   WWF Biome names for the polygons from wwf_simple.b.RData


