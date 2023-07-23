Script to reproduce Skeels and Yaxley 2023, "Functional endemism captures hotspots of unique phenotypes and restricted ranges"
Contact: alexander.skeels@gmail.com

# Set up
### Libraries
```{r}
# load in all the packages we need
library(ape)
library(raster)
library(phyloregion)
library(Matrix)
library(sp)
library(sf)
library(phylobase)
library(lmtest)
library(car)
library(maptools)
library(rgeos)
library(ggplot2)
library(scico)
library(gridExtra)
library(foreach)
library(doParallel)
```

### Directory
```{r}
# make the directory with all the data you working directory
proj_dir <- "Functional_endemism"
setwd(proj_dir)
```

# Data
All data is publically available from the following resources:

Spatial data is from Jetz W., Thomas G.H., Joy J.B., Hartmann K., Mooers A.O. 2012. The global diversity of birds in space and time. Nature. 491:444–448, based on the data from BirdLife international (https://www.birdlife.org/)

Phylogenetic data is the majority rules consensus tree produced in B. C. Weeks, S. Naeem, J. R. Lasky, and J. A. Tobias, “Diversity and extinction risk are inversely related at a global scale,” Ecology Letters, vol. 25, no. 3, pp. 697–707, 2022, from the posterior distribution of Jetz W., Thomas G.H., Joy J.B., Hartmann K., Mooers A.O. 2012. The global diversity of birds in space and time. Nature. 491:444–448.

Trait data is the the Avonet database from J. A. Tobias, C. Sheard, A. L. Pigot, et al., “Avonet: morphological, ecological and geographical data for all birds,” Ecology Letters, vol. 25, no. 3, pp. 581–597, 2022.

The AVONET1_BirdLife.csv and singe_bird_phylo.tre are the exact files from those resources and is not hosted on this repository, while the bird_spatial_sparse_matrix.rds is a transofmration of the original spatial data from BirdLife international into sparse matrix of the presence and absence data which we include in this GitHub repository for ease of access.

### Load Data
```{r}
#Load in the spatial data for all terrestrial birds, using the Jetz taxonomy. Ranges are projected onto a 110 km2 equal area projection. 
bird_sp <- readRDS("data/bird_spatial_sparse_matrix.rds")

# load trait data
bird_traits <- read.csv("data/AVONET1_BirdLife.csv")

# load phylogenetic data
bird_phy <- read.tree("data/singe_bird_phylo.tre")

# load map data for plotting
data("wrld_simpl")
```

### Match Data
```{r}
# Remove any species from the trait dataset that aren't in a community
bird_traits$Species1 <- gsub(" ", "_", bird_traits$Species1)
bird_traits <- bird_traits[bird_traits$Species1 %in% colnames(bird_sp),]
rownames(bird_traits) <- bird_traits$Species1

# subset 11 traits of interest
bird_traits <- bird_traits[, c(11:21)]

# subset species with trait and phylo information
bird_sp_coords <- bird_sp[, 1:2]
bird_sp <- bird_sp[, which(colnames(bird_sp) %in% rownames(bird_traits) & colnames(bird_sp) %in% bird_phy$tip.label)]

# check they match
dim(bird_traits)
dim(bird_sp)

# subset phylo to have species with trait and spatial
bird_phy <- drop.tip(bird_phy, bird_phy$tip.label[which(!bird_phy$tip.label %in% colnames(bird_sp))])

# write the subset tree to file
write.tree(bird_phy, file="output/clean_bird_phy.tree")
```

### Transform and Scale Data
```{r}
# log traits of interest
bird_traits_l <- log(bird_traits)

# scale traits of interest
bird_traits_s <- scale(bird_traits_l, center = T)

# create a sparse matrix representation of the presence absence data for calculations of FE and PE - saves a lot of time!!!!
bird_sparse <- as(as.matrix(bird_sp), "sparseMatrix")  
rownames(bird_sparse) <- paste0("site", 1:40680)

# save the presence/absence data
saveRDS(bird_sparse , file="output/bird_sparse_matrix.rds")

```

# Analyses
### Compare algorithms for producing a dendrogram-based representation of functional trait space 
```{r}
# run PCA to reduce colinearity among functional trait axes
bird_traits_pca <- prcomp(bird_traits_s)
summary(bird_traits_pca)

# Characterize trait space by getting Euclidean distances between species in trait space
bird_traits_dist <- dist(bird_traits_pca$x)

# Following Maire et al 2015 (doi/10.1111/geb.12299) and Mouchet et al 2008 (doi/10.1111/j.0030-1299.2008.16594.x)
# compare alternative methods to characterise functional trait space using dendrogram-based methods
dendro_1 <- as.phylo(hclust(bird_traits_dist, method="complete"))
dendro_2 <- as.phylo(hclust(bird_traits_dist, method="average"))
dendro_3 <- as.phylo(hclust(bird_traits_dist, method="mcquitty"))
dendro_4 <- as.phylo(hclust(bird_traits_dist, method="median"))
dendro_5 <- as.phylo(hclust(bird_traits_dist, method="centroid"))
dendro_6 <- as.phylo(hclust(bird_traits_dist, method="ward.D"))

# Get cophenetic distances on the dendrogram to compare to original trait distances
dendro_1_co <- cophenetic(dendro_1)
dendro_2_co <- cophenetic(dendro_2)
dendro_3_co <- cophenetic(dendro_3)
dendro_4_co <- cophenetic(dendro_4)
dendro_5_co <- cophenetic(dendro_5)
dendro_6_co <- cophenetic(dendro_6)

# standardsise the distances to have the same maximum value as the observed distances
dendro_1_co_std <- (dendro_1_co/max(dendro_1_co, na.rm=T))*max(bird_traits_dist)
dendro_2_co_std <- (dendro_2_co/max(dendro_2_co, na.rm=T))*max(bird_traits_dist)
dendro_3_co_std <- (dendro_3_co/max(dendro_3_co, na.rm=T))*max(bird_traits_dist)
dendro_4_co_std <- (dendro_4_co/max(dendro_4_co, na.rm=T))*max(bird_traits_dist)
dendro_5_co_std <- (dendro_5_co/max(dendro_5_co, na.rm=T))*max(bird_traits_dist)
dendro_6_co_std <- (dendro_6_co/max(dendro_6_co, na.rm=T))*max(bird_traits_dist)

# make sure everything in order for correlation!
bird_traits_dist <- as.matrix(bird_traits_dist)
bird_traits_dist_order <- bird_traits_dist[order(rownames(bird_traits_dist)), order(colnames(bird_traits_dist))] 
dendro_1_co_std_order <- dendro_1_co_std[order(rownames(dendro_1_co_std)), order(colnames(dendro_1_co_std))]
dendro_2_co_std_order <- dendro_2_co_std[order(rownames(dendro_2_co_std)), order(colnames(dendro_2_co_std))]
dendro_3_co_std_order <- dendro_3_co_std[order(rownames(dendro_3_co_std)), order(colnames(dendro_3_co_std))]
dendro_4_co_std_order <- dendro_4_co_std[order(rownames(dendro_4_co_std)), order(colnames(dendro_4_co_std))]
dendro_5_co_std_order <- dendro_5_co_std[order(rownames(dendro_5_co_std)), order(colnames(dendro_5_co_std))]
dendro_6_co_std_order <- dendro_6_co_std[order(rownames(dendro_6_co_std)), order(colnames(dendro_6_co_std))]

# calculate mean squared deviation (mSD) - smaller values are better
mean((bird_traits_dist_order - dendro_1_co_std_order )^2)
mean((bird_traits_dist_order - dendro_2_co_std_order )^2) # UPGMA is best
mean((bird_traits_dist_order - dendro_3_co_std_order )^2)
mean((bird_traits_dist_order - dendro_4_co_std_order )^2)
mean((bird_traits_dist_order - dendro_5_co_std_order )^2)
mean((bird_traits_dist_order - dendro_6_co_std_order )^2)


# also see if there if the alternative algorithms have a downstream effect on estimates of FE
FE_dendro_1 <- phylo_endemism(bird_sparse, dendro_1,  weighted=T)
FE_dendro_2 <- phylo_endemism(bird_sparse, dendro_2,  weighted=T)
FE_dendro_3 <- phylo_endemism(bird_sparse, dendro_3,  weighted=T)
FE_dendro_4 <- phylo_endemism(bird_sparse, dendro_4,  weighted=T)
FE_dendro_5 <- phylo_endemism(bird_sparse, dendro_5,  weighted=T)
FE_dendro_6 <- phylo_endemism(bird_sparse, dendro_6,  weighted=T)

FEs <- cbind(FE_dendro_1,
             FE_dendro_2,
             FE_dendro_3,
             FE_dendro_4,
             FE_dendro_5,
             FE_dendro_6)

cor(FEs) # they all are very highly correlated (r > 0.99) except UPGMC which is less so (r < 0.85)

# Save the best dendrogram
write.tree(dendro_2, file="output/UPGMA_functional_tree_birds.tree")
```

### Calculate Endemism Metrics
```{r}
# calculate FE and PE using the phylo_endemism function in the phyloregion package
FE_observed <- phylo_endemism(bird_sparse, dendro_2,  weighted=T)
PE_observed <- phylo_endemism(bird_sparse, bird_phy,  weighted=T)

# save the output
saveRDS(FE_observed, file='output/UPGMA_functional_endemism_global.rds')
saveRDS(PE_observed, file='output/phylogenetic_endemism_global.rds')
```

### Species Richness
```{r}
# calculate species richness
richness <- rowSums(bird_sp, na.rm=T)
```

### Create null FE and PE
```{r}
#clean up a bit
rm(list=ls())
gc()

# just load required items to save memory
bird_sparse <- readRDS("output/bird_sparse_matrix.rds")
bird_phylo <- read.tree("output/clean_bird_phy.tree")
bird_dendro <- read.tree("output/UPGMA_functional_tree_birds.tree")
FE_observed <- readRDS('output/UPGMA_functional_endemism_global.rds')
PE_observed <- readRDS('output/phylogenetic_endemism_global.rds')

# We are using a tip-shuffle null model to maintain species distributions and range sizes but randomize the functional trait value
# number of cores to utilise on your machine
num_cores <- 50 # set this to a suitable number for your machine

# Set up a parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# create an empty matrix to fill with null model data
# each column is a simulation, each row is a site

PE_null_mat <- matrix(ncol=1000, nrow=40680)
FE_null_mat <- matrix(ncol=1000, nrow=40680)
# run 1000 times - in parallel
foreach(i = 1:1000, .packages = c("phyloregion", "ape")) %dopar% {
  set.seed(i) # set seed for reproducibility - should give the exact same results as we present in the manuscript
  
  # get names
  bird_dendro_names <- bird_dendro$tip.label
  bird_phy_names <- bird_phylo$tip.label
  
  # shuffle em
  bird_dendro_names_shuffle <- bird_dendro_names[sample(1:length(bird_dendro_names))]
  bird_phy_names_shuffle <- bird_phy_names[sample(1:length(bird_phy_names))]
  
  # create new dendro/phy
  bird_dendro_shuffle <- bird_dendro
  bird_dendro_shuffle$tip.label <- bird_dendro_names_shuffle
  
  # give em the shuffled names
  bird_phy_shuffle <- bird_phylo
  bird_phy_shuffle$tip.label <- bird_phy_names_shuffle
  
  # calculate Functional Endemism on shuffled
  FE_sim <- phylo_endemism(bird_sparse, bird_dendro_shuffle,  weighted=T)
  FE_null_mat[,i] <- FE_sim
  
  # save the FE null
  saveRDS(FE_sim, file=file.path('output', 'FE_ses', paste(i, 'functional_endemism_global.rds', sep="_")))
          
 # calculate Phylogenetic Endemism
  PE_sim <- phylo_endemism(bird_sparse, bird_phy_shuffle,  weighted=T)
  PE_null_mat[,i] <- PE_sim
  
  # save the PE null
  saveRDS(PE_sim, file=file.path('output', 'PE_ses', paste(i, 'phylogenetic_endemism_global.rds', sep="_")))
}

# Stop the parallel backend
stopCluster(cl)

# again clean up
gc()
```
### Calculate the standardised effect sizes
```{r}
PE_null_mat <- matrix(ncol=1000, nrow=40680)
FE_null_mat <- matrix(ncol=1000, nrow=40680)

for(i in 1:1000){
  
  FE_sim <- readRDS(file=file.path('output', 'FE_ses', paste(i, 'functional_endemism_global.rds', sep="_")))
  PE_sim <- readRDS(file=file.path('output', 'PE_ses', paste(i, 'phylogenetic_endemism_global.rds', sep="_")))
  FE_null_mat[,i] <-  FE_sim 
  PE_null_mat[,i] <-  PE_sim 
  
}

#standardised effect sizes
PE_ses <-  (PE_observed - apply(PE_null_mat, 1, mean)) / apply(PE_null_mat, 1, sd)
FE_ses <-  (FE_observed - apply(FE_null_mat, 1, mean)) / apply(FE_null_mat, 1, sd)

# NA values are when the standard deviation is 0 - set these values to 0
PE_ses[which(is.na(PE_ses))] <- 0
FE_ses[which(is.na(FE_ses))] <- 0

saveRDS(PE_ses, file="output/PE_ses.rds")
saveRDS(FE_ses, file="output/FE_ses.rds")

gc()
```

### Bind metrics in a data frame
```{r}

# just load required items
bird_sparse <- readRDS("output/bird_sparse_matrix.rds")
bird_phylo <- read.tree("output/clean_bird_phy.tree")
bird_dendro <- read.tree("output/UPGMA_functional_tree_birds.tree")
FE_observed <- readRDS('output/UPGMA_functional_endemism_global.rds')
PE_observed <- readRDS('output/phylogenetic_endemism_global.rds')
PE_ses <- readRDS("output/PE_ses.rds")
FE_ses <- readRDS("output/FE_ses.rds")

# Get the coordinates for teh data frame
bird_sp_df <- readRDS("output/bird_sp_coords.rds")

# Get richness
richness <- rowSums(bird_sparse)

# add each metric to the spatial data
bird_sp_df$FE_obs <- FE_observed
bird_sp_df$PE_obs <- PE_observed
bird_sp_df$FE_ses <- FE_ses
bird_sp_df$PE_ses <- PE_ses
bird_sp_df$SR <- richness

# add scaled and centered metrics 
bird_sp_df$FE_obs_sc <- scale(bird_sp_df$FE_obs)
bird_sp_df$PE_obs_sc <- scale(bird_sp_df$PE_obs)
bird_sp_df$FE_ses_sc <- scale(bird_sp_df$FE_ses)
bird_sp_df$PE_ses_sc <- scale(bird_sp_df$PE_ses)
```

### Linear Models
```{r}

# What are the relationships between FE, PE, and richness
m1 <-  lm(FE_obs_sc ~ PE_obs_sc + SR, data=bird_sp_df)
m2 <-  lm(FE_obs_sc ~ PE_obs_sc, data=bird_sp_df)
m3 <-  lm(FE_obs_sc ~ SR, data=bird_sp_df)

# Variance inflation factors
vif(m1) # SR and PE are intercorrelated

# compare model fit with ANOVA
anova(m1, m2, m3) # M2 is the best fit

bird_sp_df$res_obs <- residuals(m2)
bird_sp_df$std_res_obs <- rstandard(m2)

# So FE and PE are very strongly correlated and are also correlated with SR, such that SR doesn't explain additional variance beyond PE

# This means we should look at richness corrected measures of FE - the SES
m1_ses <- lm(FE_ses_sc ~ PE_ses_sc + SR, data=bird_sp_df)
m2_ses <- lm(FE_ses_sc ~ PE_ses_sc,data=bird_sp_df)
m3_ses <- lm(FE_ses_sc ~ SR, data=bird_sp_df)

# Variance inflation factors
vif(m1_ses) # SR and PE are not intercorrelated

# compare model fit
anova(m1_ses, m2_ses, m3_ses) # M2 is the best fit

# model summary
summary(m2_ses)

# get standardised residuals
bird_sp_df$res <- residuals(m2_ses)
bird_sp_df$std_res <- rstandard(m2_ses)

# finally check the relationship between FE and FEses
m4 <-  lm(FE_ses_sc ~ FE_obs_sc, data=bird_sp_df)
summary(m4)
```

# Visualisation
### Make rasters
```{r}
# raster template
ea_template <- raster(ncol=360, nrow=113, crs=CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'),
                      xmn = -17319287,
                      ymn = -7153821,
                      xmx = 17319287,
                      ymx = 7208893)

# initialise raster
FE_ras <- ea_template
PE_ras <- ea_template
RE_ras <- ea_template

# fill rasters
values(PE_ras)  <- PE_ses
values(FE_ras)  <- FE_ses
values(RE_ras)  <- c(bird_sp_df$std_res)

# equal-area version world outline
data("wrld_simpl")
wrld_simpl_behr <- spTransform(wrld_simpl, 
                               CRSobj = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"))
wrld_simpl_behr_buff <- gBuffer(wrld_simpl_behr, width=30000)
wrld_simpl_behr_simp <- gSimplify(wrld_simpl_behr_buff, 0.1)
wrld_simpl_behr_simp <- gUnaryUnion(wrld_simpl_behr_simp)

# make sf class
wrld_simpl_behr_simp_sf <- st_as_sf(wrld_simpl_behr_simp)

# maks the rasters - as we're not interested in the oceans
PE_ras <- mask(PE_ras, wrld_simpl_behr_simp)
FE_ras <- mask(FE_ras, wrld_simpl_behr_simp)
RE_ras <- mask(RE_ras, wrld_simpl_behr_simp)
```

### Make Figure 1
```{r}
# devtools::install_github("https://github.com/coolbutuseless/svgparser")
silhouette_kiwi <- svgparser::read_svg("data/kiwi.svg")
NZ_points <- bird_sp_df[which(bird_sp_df$std_res > 7),] # all std res > 7 are found in New Zealand

p1 <- ggplot(bird_sp_df, aes(x=PE_ses, y=FE_ses, colour=SR))+
  geom_point(alpha=0.3)+
  scale_colour_scico(palette = "lapaz", direction=1)+
  stat_smooth(method="lm", colour="black")+
  annotate("text", x = -2, y = 21, 
           label = expression(paste("R"^2, "= 0.545; ", beta, " = 0.738", "\u00B1", 0.003, "; P<0.0001")), size=6) +
  labs(title = "(a)") +
  theme_bw()+
  stat_ellipse(data =  NZ_points, aes(x=PE_ses, y=FE_ses), colour = "#36454F",  size=1) +
  annotation_custom(silhouette_kiwi, xmin = 0, xmax = 3, ymin = 16, ymax = 19)

p1
# Plot FE map with world borders
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

limit <- max(abs(values(FE_ras)), na.rm=T) * c(-1, 1)

p2 <- ggplot() +
  geom_raster(data = as.data.frame(FE_ras, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico(palette = "vik", direction=1, na.value="white", limit = limit) +
  labs(title = "(b) FEses") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank())+
  labs(fill = "", x="", y="")+
  theme_bw()#+
  #geom_path(data=circleFun(c(16640370,-4774845),2000000,npoints = 100),aes(x,y), colour = "#36454F",  size=1)

p2

# Plot PE
limit <- max(abs(values(PE_ras)), na.rm=T) * c(-1, 1)

p3 <- ggplot() +
  geom_raster(data = as.data.frame(PE_ras, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico(palette = "vik", direction=1, na.value="white", limit = limit) +
  labs(title = "(c) PEses") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank())+
  labs(fill = "", x="", y="")+
  theme_bw()#+
  #geom_path(data=circleFun(c(16640370,-4774845),2000000,npoints = 100),aes(x,y), colour = "#36454F",  size=1)

p3

# Plot FE residuals map with world borders

limit <- max(abs(values(RE_ras)), na.rm=T) * c(-1, 1)

p4 <- ggplot() +
  geom_raster(data = as.data.frame(RE_ras, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico(palette = "cork", direction=-1, na.value="white", limit = limit) +
  labs(title = "(d) FEses ~ PEses std. residuals") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        x.axis.text = element_blank(),
        x.axis.title = element_blank())+
  theme_bw()+
  labs(fill = "", x="", y="")#+
  #geom_path(data=circleFun(c(16640370,-4774845),2000000,npoints = 100),aes(x,y), colour = "#36454F",  size=1)

grid.arrange(p1, p2, p3,p4, layout_matrix = matrix(c(1,1,1,2,3,4), nrow=3))
```

# Supplementary Materials

```{r}

### Make Figure S1
pS1_a <- ggplot(bird_sp_df, aes(x=PE_obs, y=FE_obs))+
  geom_point(alpha=0.3)+
  scale_colour_scico(palette = "lapaz", direction=1)+
  stat_smooth(method="lm", colour="red")+
  labs(title = "(a)") +
  theme_bw()+
  labs(x="Phylogenetic endemism", y="Functional endemism")


pS1_b <- ggplot(bird_sp_df, aes(x=SR, y=FE_obs))+
  geom_point(alpha=0.3)+
  scale_colour_scico(palette = "lapaz", direction=1)+
  stat_smooth(method="lm", colour="red")+
  labs(title = "(b)") +
  theme_bw()+
  labs(x="Species richness", y="Functional endemism")

grid.arrange(pS1_a, pS1_b)

### Make Figure S2

# raster template
ea_template <- raster(ncol=360, nrow=113, crs=CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'),
                      xmn = -17319287,
                      ymn = -7153821,
                      xmx = 17319287,
                      ymx = 7208893)

# initialise raster
FE_obs_ras <- ea_template
PE_obs_ras <- ea_template
SR_obs_ras <- ea_template

# fill rasters
values(PE_obs_ras)  <- PE_observed
values(FE_obs_ras)  <- FE_observed
values(SR_obs_ras)  <- richness

# maks the rasters - as we're not interested in th eoceans
PE_obs_ras <- mask(PE_obs_ras, wrld_simpl_behr_simp)
FE_obs_ras <- mask(FE_obs_ras, wrld_simpl_behr_simp)
SR_obs_ras <- mask(SR_obs_ras, wrld_simpl_behr_simp)

# Plot FE map 
limit <- max(abs(values(FE_obs_ras)), na.rm=T) * c(-1, 1)

pS2_a <- ggplot() +
  geom_raster(data = as.data.frame(FE_obs_ras, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico(palette = "vik", direction=1, na.value="white", limit = limit) +
  labs(title = "(a)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank())+
  labs(fill = "FE", x="", y="")+
  theme_bw()


# Plot FE map 
limit <- max(abs(values(PE_obs_ras)), na.rm=T) * c(-1, 1)

pS2_b <- ggplot() +
  geom_raster(data = as.data.frame(PE_obs_ras, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico(palette = "vik", direction=1, na.value="white", limit = limit) +
  labs(title = "(b)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank())+
  labs(fill = "PE", x="", y="")+
  theme_bw()

# Plot FE map 
limit <- max(abs(values(SR_obs_ras)), na.rm=T) * c(-1, 1)

pS2_c <- ggplot() +
  geom_raster(data = as.data.frame(SR_obs_ras, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico(palette = "vik", direction=1, na.value="white", limit = limit) +
  labs(title = "(c)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank())+
  labs(fill = "SR", x="", y="")+
  theme_bw()

grid.arrange(pS2_a, pS2_b, pS2_c)

### Make Figure S3

pS3 <- ggplot(bird_sp_df, aes(x=SR, y=FE_ses))+
  geom_point(alpha=0.3)+
  scale_colour_scico(palette = "lapaz", direction=1)+
  stat_smooth(method="lm", colour="red")+
  theme_bw()+
  labs(x="Species richness", y="Functional endemism (SES)")

### Make Figure S4

pS4 <- ggplot(bird_sp_df, aes(x=FE_obs, y=FE_ses))+
  geom_point(alpha=0.3)+
  scale_colour_scico(palette = "lapaz", direction=1)+
  stat_smooth(method="lm", colour="red")+
  theme_bw()+
  labs(x="Functional Endemism", y="Functional endemism (SES)")


### Make Figure S5
RE_ras_pos <- Which(RE_ras > 2)
RE_ras_neg <- Which(RE_ras < -2)
values(RE_ras_pos)[which(values(RE_ras_pos)==0)] <- NA
values(RE_ras_neg)[which(values(RE_ras_neg)==0)] <- NA

pS5_a <- ggplot() +
  geom_raster(data = as.data.frame(RE_ras_pos, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico_d(palette = "vik", direction=1, na.value="white") +
  labs(title = "(a)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        x.axis.text = element_blank(),
        x.axis.title = element_blank())+
  theme_bw()+
  labs(fill = "std. residual > 2", x="", y="")#+

pS5_b <- ggplot() +
  geom_raster(data = as.data.frame(RE_ras_neg, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black", linewidth=0.001)+
  scale_fill_scico_d(palette = "vik", direction=1, na.value="white") +
  labs(title = "(b)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        x.axis.text = element_blank(),
        x.axis.title = element_blank())+
  theme_bw()+
  labs(fill = "std. residual < -2", x="", y="")#+

grid.arrange(pS5_a, pS5_b)
```

