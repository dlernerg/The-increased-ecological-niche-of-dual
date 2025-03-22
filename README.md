# Scripts and Data for "The increased environmental niche of dual mycorrhizal tree species"

- **Author**: David Lerner
- **Email**: dlernerg@gmail.com

## Version

- **Current Version**: v2.2
- **Release Date**: March 22, 2025

## Description

This repository contains the R scripts used for the analysis in the paper **The increased environmental niche of dual mycorrhizal tree species** by Rog et al. These scripts are designed to reproduce the analyses described in the paper.

The necessary data can be found in a supplementary DOI: https://zenodo.org/records/15065494

## R Version

This project was developed using R version 4.4.2.

## R Package Versions

This project uses the following R packages and their respective versions:

- `ggplot2`: 3.3.5
- `phytools`: 2.4.4
- `ape`: 5.8.1
- `sf`: 1.0.19
- `raster`: 3.6.31
- `ggtree`: 3.14.0
- `ggnewscale`: 0.5.1
- `ggstar`: 1.0.4
- `caper`: 1.0.3
- `phylolm`: 2.5.6

### Scripts

1. **rogetal_1.R**:  
   This script performs the initial data cleaning and preprocessing. It reads in the `confirmed.csv`, `AM_confirmed.csv`, and `EM_confirmed.csv` files, combines the confirmed tree species with a phylogenetic tree from Sanchez Martinez et al., (2020) and outputs cleaned data that will be used in subsequent analysis.

All the following scripts are dependent on rogetal_1.R

2. **rogetal_2.R**:  
   This script conducts the phylogenetic analyses observed in figure 1 (phylogenetic map plotting, MPD and delta statistics), using the preprocessed data from `rogetal_1.R` file. 

3. **rogetal_3.R**:  
   This script runs the biome overlap data observed in Figure 2. It uses the cleaned data and `names_biome.csv` and `wwf_simple.b.geojson` to generate binomial and poisson models of the occupided space for each of the mycorrhizal associating tree species groups.

4. **rogetal_4.R**:  
   This script performs the final analysis of the paper (Figure 3) - an analysis on the niche space. It combines the data generated in rogetal_1.R with environmental data from WorldClim data (Hijmans et al. 2015) and Olsen P data (McDowell et al. 2023).

5. **rogetal_5_sims2**
   This script performs the simulations of the paper. It uses the ouptut from rogetal_3.R with either the generated confirmed and unconfirmed model summary dataframes "model_confirmed" and "model_unconfirmed", respectively.
     
### Data Files

1. **confirmed.csv**:  
   List of confirmed dual mycorrhizal genera from Teste et al., 2020.

2. **AM_confirmed.csv**:  
   List of confirmed AM mycorrhizal genera from Soudzilovskaia et al., 2020.

3. **EM_confirmed.csv**:  
   List of confirmed EM mycorrhizal genera from Soudzilovskaia et al., 2020.

4. **global.inland.geojson**:  
   A shapefile containing grid cell polygons of the world map, created using the 'dggridR' package, which employs the ISEA Discrete Global Grid System for equal-sized grid projection on a 2-dimensional plane (described in detail in Lerner et al., (2023)).

5. **wwf_simple.b.csv**:  
   Contains smoothened polygons for the WWF biomes, as described in Lerner et al., 2023.

6. **names_biomes.csv**:  
   Table containing the names of WWF Biome, aligned with the biome polygons in the "wwf_simple.b.geojson" shapefile. 

7. **GBIFpolygon_groups.geojson**
   This shapefile containing the species distributions data from Lerner et al., 2023. Each record represents a species' population, with the following columns: 1. unique identifier index, 2. species ID, 3. distribution polygons.  

8. **SanchezM (2020).tre**
   The phylogenetic tree obtained from Sanchez Martinez et al., 2020

9. **code_delta.R**
    R code to run delta statistic obtained from Borges et al., 2019

10. **meansBioClim2.csv**
    This file contains the average BioClim variables (Hijmans et al., 2005) for each world grid cell. Each row corresponds to a specific grid cell from the global.inland.geojson shapefile, with the 19 columns representing the mean values of the BioClim variables for that cell.
    
11. **P_olsen_df.csv**
    A list of average olsen phosporous measurements (McDowell et al., 2023). The order of the values correspond to the order of the world grid cells (obtained from global.inland.geojson) 

