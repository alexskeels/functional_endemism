Function Endemism
================
Alexander Skeels & Keaghan Yaxley
2023-05-04

# Set up

### Libraries

``` r
library(ape)
library(raster)
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
```

### Directory

``` r
proj_dir <- "Functional_endemism"
```

### Functions

``` r
calc_PE <- function(tree, sites_x_tips,presence=c("presence","abundance","probability")) {
  
  # add code to check that the values are correct for the presence type:
  # 0 or 1 for presence - this calculates PE (Rosauer et al 2009)
  # from 0 to 1 for probability - this calculates model weighted PE (Rosauer, in prep)
  # any value for abundance - this calculation is equivalent to BED (Cadotte & Davies 2010)
  
  #default value for presence
  if (is.na(presence)) {presence="presence"}
  
  # change to a phylobase phylo4 object
  if (class(tree) == "phylo") {tree <- phylo4(tree)}
  
  sites_x_branches <- data.frame(cbind(rep(0,nrow(sites_x_tips))))
  
  for (i in 1:nTips(tree)) {
    sites_x_branches[,i] <- sites_x_tips[,which(labels(tree)[i]==names(sites_x_tips))]
    names( sites_x_branches)[i] <- labels(tree)[i]
  }
  
  rm(sites_x_tips); #gc()
  branch.labels <- as.character(labels(tree))
  branch.count <- length(labels(tree))
  
  # add names and occupancy columns for internal branches
  for (i in (nTips(tree)+1):branch.count) {
    branch.labels[i] <- paste("b",i,sep="")
    desc <- as.integer(descendants(tree,i, type="tips"))
    if (presence=="abundance") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=sum))
    } else if (presence=="presence") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=max))
    } else if (presence=="probability") {
      branch_col <- as.numeric(apply(sites_x_branches[,desc],MARGIN=1,FUN=parent.prob))
    }
    sites_x_branches[,i] <- branch_col
    names(sites_x_branches[i]) <- branch.labels[i]
    #cat(i,branch.labels[i],length(desc),"\n")
  #  gc(verbose=F)
  }
  
  #scale columns (branches) to sum to 1
  sites_x_branches <- apply(sites_x_branches,MARGIN=2,FUN=scale.to,1)
  
  #now scale branches to sum to their length
  branch.lengths <- as.numeric(edgeLength(tree,1:branch.count))
  branch.lengths[is.na(branch.lengths)] <- 0
  for (i in 1:length(branch.lengths)) {
    sites_x_branches[,i] <- sites_x_branches[,i] * branch.lengths[i]
  }
  
  PE.vec <- apply(sites_x_branches,MARGIN=1,FUN=sum,na.rm=T)
  PE <- data.frame(cbind(1:nrow(sites_x_branches),PE.vec))
  names(PE) <- c("site","PE")
  return(PE)
}

parent.prob <- function(probabilities) {
  # probabilities is a vector of values between 0 and 1
  # add code to check values are of correct type!
  parent.prob <- 1 - prod(1-probabilities)
  return(parent.prob)
}

scale.to <- function(vec,vec.sum) {
  #mat is a vector
  #this function rescales each the vector values to sum to 'vec.sum'
  vec.tot <- sum(vec,na.rm=TRUE)
  if (vec.tot > 0) {
    vec.out <- vec.sum*vec/vec.tot
  } else {
    vec.out <- rep(0,times=length(vec))  #should columns for tips / branches with no occurrences be removed?
  }
  return(vec.out)
}

weighted_endemism <- function(x){
  #species_richness <- rowSums(x)
  proportions_matrix <- apply(x, 2, FUN=function(y){y <- y/sum(y)})
  WE <- rowSums(proportions_matrix, na.rm=T)
  WE_df <- data.frame(site =1:length(WE), WE=WE)
  return(WE_df)
}
```

# Data

### Load Data

``` r
#Load in the spatial data for all terrestrial birds, using the Jetz taxonomy. Ranges are projected onto a 110 km2 equal area projection. 
bird_sp <- read.csv("data/bird_spatial_data.csv")

# load trait data
bird_traits <- read.csv("data/AVONET1_BirdLife.csv")

# load phylogenetic data
bird_phy <- read.tree("data/singe_bird_phylo.tre")

# load map data
data("wrld_simpl")
```

### Match Data

``` r
# Remove any species from the trait dataset that aren't in a community
bird_traits$Species1 <- gsub(" ", "_", bird_traits$Species1)
bird_traits <- bird_traits[bird_traits$Species1 %in% colnames(bird_sp),]
rownames(bird_traits) <- bird_traits$Species1

# subset traits of interest
bird_traits <- bird_traits[, c(11:21)]

# subset species with trait and phylo information
bird_sp_coords <- bird_sp[, 1:2]
bird_sp <- bird_sp[, which(colnames(bird_sp) %in% rownames(bird_traits) & colnames(bird_sp) %in% bird_phy$tip.label)]

# check they match
dim(bird_traits)
dim(bird_sp)

# subset phylo to have species with trait and spatial
bird_phy <- drop.tip(bird_phy, bird_phy$tip.label[which(!bird_phy$tip.label %in% colnames(bird_sp))])
```

### Transform and Scale Data

``` r
# log traits of interest
bird_traits_l <- log(bird_traits)

# scale traits of interest
bird_traits_s <- scale(bird_traits_l, center = T)
```

# Analyses

### Run PCoA

``` r
# run PCoA
bird_traits_pca <- prcomp(bird_traits_s)
summary(bird_traits_pca)

# Get distances between species
bird_traits_dist <- dist(bird_traits_pca$x)

# Cluster distances to get a dendrogram
bird_traits_dendro <- as.phylo(hclust(bird_traits_dist))

# Save the dendrogram
write.tree(bird_traits_dendro, file="output/functional_tree_birds.tree")
```

### Calculate Endemism Metrics

``` r
# calculate Functional Endemism
FE <- calc_PE(bird_traits_dendro, bird_sp, "presence")
saveRDS(FE, file='output/functional_endemism_global.rds')

# calculate Weighted Endemism
WE <- weighted_endemism(bird_sp)
saveRDS(WE, file='output/weighted_endemism_global.rds')

# calculate Phylogenetic Endemism
PE <- calc_PE(bird_phy, bird_sp, "presence")
saveRDS(PE, file='output/phylogenetic_endemism_global.rds')
```

### Species Richness

``` r
richness <- rowSums(bird_sp[, 3:ncol(bird_sp)], na.rm=T)
```

### Bind metrics

``` r
# Get the coordinates
bird_sp_df <- bird_sp[, 1:2]

# add each metric to the spatial data
bird_sp_df$FE <- FE$PE
bird_sp_df$WE <- WE$WE
bird_sp_df$PE <- PE$PE
bird_sp_df$SR <- richness

# add scaled and centered metrics 
bird_sp_df$FE_sc <- scale(bird_sp_df$FE)
bird_sp_df$WE_sc <- scale(bird_sp_df$WE)
bird_sp_df$PE_sc <- scale(bird_sp_df$PE)
bird_sp_df$SR_sc <- scale(bird_sp_df$SR)
```

### Linear Models

``` r
# analysis
m1 <- lm(FE_sc ~ WE_sc + PE_sc, data=bird_sp_df)
m2 <- lm(FE_sc ~ WE_sc,         data=bird_sp_df)
m3 <- lm(FE_sc ~ PE_sc,         data=bird_sp_df)
m4 <- lm(FE_sc ~ 1,             data=bird_sp_df)

# Variance inflation factors
vif(m1)

# see which model is the best fit
lrtest(m1,m2)
lrtest(m1,m3)
lrtest(m1,m4)

# model summary
summary(m3)

# get standardised residuals
bird_sp_df$res <- residuals(m3)
bird_sp_df$std_res <- rstandard(m3)
```

# Visualisation

### Make rasters

``` r
# equal area raster  template
ea_template <- raster(ncol=360, nrow=113, crs=CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs'),
                      xmn = -17319287,
                      ymn = -7153821,
                      xmx = 17319287,
                      ymx = 7208893)

# initialise rasters
FE_ras <- ea_template
WE_ras <- ea_template
PE_ras <- ea_template
SR_ras <- ea_template
RE_ras <- ea_template

# fill rasters
values(SR_ras)  <- bird_sp_df$SR
values(FE_ras)  <- bird_sp_df$FE
values(WE_ras)  <- bird_sp_df$WE
values(PE_ras)  <- bird_sp_df$PE
values(RE_ras)  <- c(bird_sp_df$std_res)

# equal-area version world outline
wrld_simpl_behr <- spTransform(wrld_simpl, 
                               CRSobj = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"))
wrld_simpl_behr_buff <- gBuffer(wrld_simpl_behr, width=0.75)
wrld_simpl_behr_simp <- gSimplify(wrld_simpl_behr_buff, 0.1)
wrld_simpl_behr_simp <- gUnaryUnion(wrld_simpl_behr_simp)

# make sf class
wrld_simpl_behr_simp_sf <- st_as_sf(wrld_simpl_behr_simp)
```

### Make Figure 1

``` r
# devtools::install_github("https://github.com/coolbutuseless/svgparser")
silhouette_kiwi <- svgparser::read_svg("data/kiwi.svg")

p1 <- ggplot(bird_sp_df, aes(x=PE, y=FE, colour=SR))+
  geom_point(alpha=0.5)+
  scale_colour_scico(palette = "acton", direction=-1)+
  stat_smooth(method="lm", colour="black")+
  annotate("text", x = 20, y = 4.5, 
           label = expression(paste("R"^2, "= 0.96; ", beta, " = 0.98", "\u00B1", 0.001, "; P<0.0001")), size=4) +
  labs(title = "(a)") +
  theme_bw()+
  geom_point(data =  data.frame(PE = 15.81307, FE = 2.378979, SR=0), aes(x = PE, y = FE), shape = 21, size = 5, colour = "#FF69B4", fill=NA, stroke=2) +
  annotation_custom(silhouette_kiwi, xmin = 10, xmax = 15, ymin = 2.5, ymax = 3.5)

# Plot FE map with world borders
p2 <- ggplot() +
  geom_raster(data = as.data.frame(FE_ras, xy=T), aes(x = x, y = y, fill = layer)) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black")+
  scale_fill_scico(palette = "batlowW", direction=-1) +
  labs(title = "(b) FE") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank())+
  labs(fill = "", x="", y="")+
  theme_bw()+
  geom_point(data =  data.frame(x = 16161452 , y = -5216668, SR=0), aes(x = x, y = y), shape = 21, size = 5, colour = "#FF69B4", fill=NA, stroke=2) 
       

# Plot FE residuals map with world borders
RE_df <- as.data.frame(RE_ras > 2 , xy=T)
RE_df$layer <- as.numeric(RE_df$layer )
p3 <- ggplot() +
  geom_raster(data = RE_df, aes(x = x, y = y, fill = as.factor(layer))) +
  geom_sf(data = wrld_simpl_behr_simp_sf, fill = NA, color = "black")+
  scale_fill_manual(values = c("white", "#CC5500" ), guide = guide_legend(reverse = TRUE)) +
  labs(title = "(c) FE~ PE std. residuals > 2") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        x.axis.text = element_blank(),
        x.axis.title = element_blank())+
  theme_bw()+
  labs(fill = "", x="", y="")+
  geom_point(data =  data.frame(x = 16161452 , y = -5216668, SR=0), aes(x = x, y = y), shape = 21, size = 5, colour = "#FF69B4", fill=NA, stroke=2) 

# Arrange maps side by side
grid.arrange(p1, p2, p3, layout_matrix = matrix(c(1,1,2,3), nrow=2))
```
