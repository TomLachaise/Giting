#### Packages ####
#
setwd("C:/Users/tomla/Dropbox/RootFun/Data_scripts/Exploratories/Stats/FunRoot") # Change address to set to source file location
#
source("HighstatLibV12.R")
source("ggpub.R")
require(readr)
require(ggplot2)
require(gghighlight)
require(data.table)
require(GrapheR)
require(tidyr)
require(doBy)
require(dplyr)
require(plotrix) # std.error()
require(robustHD)
require(scales)
require(formattable)
require(gridExtra)
require(magicfor)
require(GGally)
require(rgl)
require(car)
require(crunch)
require(ggthemes)
require(scales)  # for percentage scales

#### Data loading ####

# Load the datasets of LM results

# Across Explo
# CWM Raw
CWM_Results_Across_Raw_coef  <- read.table(
  "CWM_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
# PCA SLA
PCA.SLA_Results_Across_Raw_coef  <- read.table(
  "PCA.SLA_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
# PCA SLA other axes
PCA.SLA_norot_Results_Across_Raw_coef  <- read.table(
  "PCA.SLA_norot_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",") 
# PCA Belowground
PCA.Belowground_Results_Across_Raw_coef <- read.table(
  "PCA.Belowground_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",") 
# PCA SLA guild
PCA.SLA.Guild_Results_Across_Raw_coef <- read.table(
  "PCA.SLA.Guild_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
# PCA Belowground guild
PCA.Belowground.Guild_Results_Across_Raw_coef <- read.table(
  "PCA.Belowground.Guild_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
# pPCA Belowground
PCA.p_SLA.Guild_Results_Across_Raw_coef  <- read.table(
  "PCA.p_SLA.Guild_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
# Guild
Guild_Results_Across_Raw_coef <- read.table(
  "Guild_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
# CWM per slice
CWM_per_slice <- read.table(
  "CWM_per_slice.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")

#### 1) Forest plot of coefficients ####

#### 1.1) CWM Raw ####

levels(CWM_Results_Across_Raw_coef$trait)

# Trait raw
trait_raw_names  <-  c("CWM_LEDA_SLA",
                       "CWM_Morpho_Root.weight.ratio_mean", 
                       "CWM_Klim_BBRsum",                              
                       "CWM_Morpho_Root.tissue.density.Rose_mean",
                       "CWM_Morpho_Specific.root.length_mean",
                       "CWM_Morpho_Root.branching.intensity_mean",
                       "CWM_Morpho_Root.1Order.average.diameter_mean",
                       "CWM_HM_col",    
                       "CWM_Iso_Root_N_pc_mean",          
                       "CWM_RDepth_RD50_mean"
                       )

# filter out phylo traits

CWM_Results_Across_Raw_coef_Raw <- dplyr::filter(CWM_Results_Across_Raw_coef, trait %in% trait_raw_names)
CWM_Results_Across_Raw_coef_Raw <- droplevels(CWM_Results_Across_Raw_coef_Raw)

# Remove the intercept and calculate the bounds for the x scale
CWM_Results_Across_Raw_coef_Raw <- dplyr::filter(CWM_Results_Across_Raw_coef_Raw, term != "(Intercept)")
CWM_Results_Across_Raw_coef_Raw <- data.table(CWM_Results_Across_Raw_coef_Raw)
CWM_Results_Across_Raw_coef_Raw[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
CWM_Results_Across_Raw_coef_Raw[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the CWM traits names
trait.name <- CWM_Results_Across_Raw_coef_Raw$trait
trait.name <- unique(trait.name)
levels(trait.name)

# Rearrange the order of CWM traits from left to right
CWM_Results_Across_Raw_coef_Raw$trait <- factor(CWM_Results_Across_Raw_coef_Raw$trait,
          levels= c("CWM_LEDA_SLA",                                 "CWM_Morpho_Root.weight.ratio_mean",
                    "CWM_Klim_BBRsum",                              "CWM_Morpho_Root.tissue.density.Rose_mean",
                    "CWM_Morpho_Specific.root.length_mean",         "CWM_Morpho_Root.branching.intensity_mean",
                    "CWM_Morpho_Root.1Order.average.diameter_mean", "CWM_HM_col",
                    "CWM_Iso_Root_N_pc_mean",                       "CWM_RDepth_RD50_mean"))

# The delta R2 are in FD_models_Env_drivers_traits.R

trait.labels<-c("CWM Specific leaf area \n R² = 0.53", "CWM Root weight ratio \n R² = 0.59",
                "CWM Bud-bank size \n R² = 0.35", "CWM Root tissue density \n R² = 0.21",
                "CWM Specific root length \n R² = 0.26", "CWM Root branching intensity \n R² = 0.44",
                "CWM Fine roots diameter \n R² = 0.44", "CWM Mycorrhizal colonization \n R² = 0.41",
                "CWM Fine roots %N \n R² = 0.25", "CWM Rooting depth 50% \n R² = 0.48")

names(trait.labels) <- levels(CWM_Results_Across_Raw_coef_Raw$trait)

# Give explicit names to environments
levels(CWM_Results_Across_Raw_coef_Raw$term)
levels(CWM_Results_Across_Raw_coef_Raw$term) <- c("Intercept",
                            "Region Hainich", "Region Schorfheide",
                            "C:N", "delta-N-15", "Land-use-intensity",
                            "NH4","NO3", "N:P", "pH", "Phosphorus",
                            "Sand content", "Moisture")

# Remove unused level (Intercept)
term <- droplevels(CWM_Results_Across_Raw_coef_Raw$term)
CWM_Results_Across_Raw_coef_Raw$term <- droplevels(CWM_Results_Across_Raw_coef_Raw$term)

# Reorganize per order of importance (larger spatial grain at top, soil chemicals at bottom)
 CWM_Results_Across_Raw_coef_Raw$term <- factor(term, levels = c(
   "delta-N-15",
   "NH4",
   "NO3",
   "Phosphorus",
   "N:P",
   "C:N",
   "pH",
   "Sand content",
   "Moisture",
   "Land-use-intensity",
   "Region Schorfheide",
   "Region Hainich"))
 
# Define colors based on significance from the p.values
CWM_Results_Across_Raw_coef_Raw$direction[CWM_Results_Across_Raw_coef_Raw$estimate > 0]   <-"pos"
CWM_Results_Across_Raw_coef_Raw$direction[CWM_Results_Across_Raw_coef_Raw$estimate < 0]   <-"neg"
CWM_Results_Across_Raw_coef_Raw$direction[CWM_Results_Across_Raw_coef_Raw$p.value > 0.05] <-"NS"

# Forest Plot #

png("p1_Forest_CWM_Raw.png", width = 3200, height = 4000, units = 'px', res = 300)

p1_forest_frequentist = ggplot(data=CWM_Results_Across_Raw_coef_Raw,
                               aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  # Color based on significane and direction
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol = 5, nrow = 2, scales = "free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=9,face="bold"), # Title theme
        axis.text.y = element_text(size=9), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=7), # Estimates strength labels
        axis.title = element_text(size=9,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 9, hjust=0, vjust = 1, angle=180, face="bold"),
        strip.text.x = element_text(size = 9, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p1_forest_frequentist

dev.off()

#### 1.2) CWM Phylo ####

# Trait phylo
trait_phylo_names <- c("CWM_p_Klim_BBRsum_est",                              
                       "CWM_p_LEDA_Height_est",                             
                       "CWM_p_LEDA_SLA_est",                                
                       "CWM_p_Morpho_Root.1Order.average.diameter_mean_est",
                       "CWM_p_Morpho_Root.tissue.density.Rose_mean_est",    
                       "CWM_p_Morpho_Root.weight.ratio_mean",           
                       "CWM_p_Morpho_Specific.root.length_mean_est",        
                       "CWM_p_RDepth_Maximum.rooting.depth_mean_est",       
                       "CWM_p_Seed.weight" )

# filter out phylo traits

CWM_Results_Across_Raw_coef_Phylo <- dplyr::filter(CWM_Results_Across_Raw_coef, trait %in% trait_phylo_names)

# Remove the intercept and calculate the bounds for the x scale
CWM_Results_Across_Raw_coef_Phylo <- dplyr::filter(CWM_Results_Across_Raw_coef_Phylo, term != "(Intercept)")
CWM_Results_Across_Raw_coef_Phylo <- data.table(CWM_Results_Across_Raw_coef_Phylo)
CWM_Results_Across_Raw_coef_Phylo[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
CWM_Results_Across_Raw_coef_Phylo[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name<-CWM_Results_Across_Raw_coef_Phylo$trait
trait.name<-unique(trait.name)
trait.name

# Rearrange the order of traits from aboveground to belowground
CWM_Results_Across_Raw_coef_Phylo$trait <- factor(CWM_Results_Across_Raw_coef_Phylo$trait,
                levels= c("CWM_p_LEDA_SLA_est","CWM_p_LEDA_Height_est","CWM_p_Seed.weight",
                          "CWM_p_Morpho_Specific.root.length_mean_est","CWM_p_Morpho_Root.1Order.average.diameter_mean_est",
                          "CWM_p_Morpho_Root.tissue.density.Rose_mean_est","CWM_p_RDepth_Maximum.rooting.depth_mean_est",
                          "CWM_p_Morpho_Root.weight.ratio_mean","CWM_p_Klim_BBRsum_est"))

# The delta R2 are 

# Phylo R squared
trait.labels<-c("Specific leaf area \n 0.43","Height \n 0.25", "Seed weight \n 0.07",
                "Specific root length \n 0.32","Fine roots diameter \n 0.15",
                "Root tissue density \n 0.22","Maximum rooting depth \n 0.65",
                "Root weight ratio \n 0.37","Bud bank size \n 0.43")

names(trait.labels)<-levels(CWM_Results_Across_Raw_coef_Phylo$trait)

# Give explicit names to traits
CWM_Results_Across_Raw_coef_Phylo$term<-as.factor(CWM_Results_Across_Raw_coef_Phylo$term)
levels(CWM_Results_Across_Raw_coef_Phylo$term)
levels(CWM_Results_Across_Raw_coef_Phylo$term)<- c("Intercept",
                                         "Region Hainich","Region Schorfheide",
                                         #"Aspect Eastness", "Aspect Northness", "Bacteria",
                                         "Bulk density",
                                         "C/N ratio", "d15N soil",
                                         "Fertilization",
                                         "Grazing","Mowing","NH4","NO3", "N/P ratio",
                                         "pH", "Silt", "Slope",
                                         "Soil Phosphorus", "Soil moisture", "Air temperature" #"Topographic Wetness Index"
                                         )

# OPTIONAL: Remove quadratic effects
# CWM_Results_Across_Raw_coef_Phylo<-CWM_Results_Across_Raw_coef_Phylo[!grepl("?", CWM_Results_Across_Raw_coef_Phylo$term),]

# Remove unused level (Intercept)
term <- droplevels(CWM_Results_Across_Raw_coef_Phylo$term)
CWM_Results_Across_Raw_coef_Phylo$term <- droplevels(CWM_Results_Across_Raw_coef_Phylo$term)

# Reorganize per order of importance (most response at top, least at bottom)
levels(CWM_Results_Across_Raw_coef_Phylo$term)
CWM_Results_Across_Raw_coef_Phylo$term <- factor(term, levels = c(
  "Region Schorfheide",
  "Region Hainich",
  "Slope",
  "Air temperature",
  "Soil moisture",
  "C/N ratio",
  "N/P ratio",
  "Soil Phosphorus",
  "d15N soil",
  "NH4",
  "NO3",
  "pH",
  "Silt",
  "Bulk density",
  "Grazing",
  "Mowing",
  "Fertilization"))

# Define colors based on significance from the p.values
CWM_Results_Across_Raw_coef_Phylo$direction[CWM_Results_Across_Raw_coef_Phylo$estimate > 0]    <-"pos"
CWM_Results_Across_Raw_coef_Phylo$direction[CWM_Results_Across_Raw_coef_Phylo$estimate < 0]    <-"neg"
CWM_Results_Across_Raw_coef_Phylo$direction[CWM_Results_Across_Raw_coef_Phylo$p.value  > 0.05] <-"NS"

# Forest Plot #

png("p2_Forest_CWM_Phylo.png", width = 4000, height = 2500, units = 'px', res = 300)

p2_forest_frequentist = ggplot(data=CWM_Results_Across_Raw_coef_Phylo,
                               aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.2)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.2)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.2)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title=element_text(size=6,face="bold"), # Title theme
        axis.text.y=element_text(size=6), # Traits labels
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=4), # Estimates labels
        axis.title=element_text(size=6,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 6, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 6, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p2_forest_frequentist

dev.off()

#### 1.3) PCA SLA ####

PCA.SLA_Results_Across_Raw_coef
# Trait raw
trait_raw_names  <-  c("PCA.SLA_PC1", "PCA.SLA_PC2")

# filter out phylo traits
PCA.SLA_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
PCA.SLA_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA_Results_Across_Raw_coef, term != "(Intercept)")
PCA.SLA_Results_Across_Raw_coef <- data.table(PCA.SLA_Results_Across_Raw_coef)
PCA.SLA_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
PCA.SLA_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- PCA.SLA_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
PCA.SLA_Results_Across_Raw_coef$trait <- factor(PCA.SLA_Results_Across_Raw_coef$trait,
                                          levels= c("PCA.SLA_PC1", "PCA.SLA_PC2"))
# The delta R2 are 
# Imputed
trait.labels<-c("Dimension 1 \n Collaboration \n R² = 0.46",
                "Dimension 2 \n Conservation and regeneration \n R² = 0.67")
#
names(trait.labels) <- levels(PCA.SLA_Results_Across_Raw_coef$trait)
# Remove unused level (Intercept)
PCA.SLA_Results_Across_Raw_coef$term <- as.factor(PCA.SLA_Results_Across_Raw_coef$term)
PCA.SLA_Results_Across_Raw_coef$term <- droplevels(PCA.SLA_Results_Across_Raw_coef$term)

# Give explicit names to Env
levels(PCA.SLA_Results_Across_Raw_coef$term)
levels(PCA.SLA_Results_Across_Raw_coef$term)<- c(
  "Region Hainich","Region Schorfheide",
       "C:N", "delta-N-15", "Land-use-intensity", "NH4", "NO3", "N:P",
       "pH", "Phosphorus", "Sand content", "Moisture")
# Reorganize per order of importance (global to proximal)
PCA.SLA_Results_Across_Raw_coef$term <- ordered(
  PCA.SLA_Results_Across_Raw_coef$term,
  levels = c("delta-N-15", "NH4", "NO3", "Phosphorus", "N:P", "C:N", "pH", "Sand content",
             "Moisture", "Land-use-intensity", "Region Schorfheide", "Region Hainich"))

# Define colors based on significance from the p.values
PCA.SLA_Results_Across_Raw_coef$direction[PCA.SLA_Results_Across_Raw_coef$estimate > 0]   <- "pos"
PCA.SLA_Results_Across_Raw_coef$direction[PCA.SLA_Results_Across_Raw_coef$estimate < 0]   <- "neg"
PCA.SLA_Results_Across_Raw_coef$direction[PCA.SLA_Results_Across_Raw_coef$p.value > 0.05] <- "NS"

# Forest Plot #

png("p3_Forest.plot_PCA.SLA.png", width = 1500, height = 2000, units = 'px', res = 300)

p3_Forest.plot_PCA.SLA = ggplot(data=PCA.SLA_Results_Across_Raw_coef,
                               aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.3)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=9,face="bold"), # Title theme
        axis.text.y = element_text(size=9), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=9), # Estimates labels
        axis.title = element_text(size=9,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 9, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 9, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p3_Forest.plot_PCA.SLA

dev.off()
#### 1.3.2) PCA SLA_norot other axes ####

PCA.SLA_norot_Results_Across_Raw_coef

# Remove the intercept and calculate the bounds for the x scale
PCA.SLA_norot_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA_norot_Results_Across_Raw_coef, term != "(Intercept)")
PCA.SLA_norot_Results_Across_Raw_coef <- data.table(PCA.SLA_norot_Results_Across_Raw_coef)
PCA.SLA_norot_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
PCA.SLA_norot_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]


# Keep the axes for up to 90% of variance explained (93.5%) so Axes 3 to 6
PCA.SLA_norot_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA_norot_Results_Across_Raw_coef, trait != "CWM_SLA_norot_Axis7")
PCA.SLA_norot_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA_norot_Results_Across_Raw_coef, trait != "CWM_SLA_norot_Axis8")
PCA.SLA_norot_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA_norot_Results_Across_Raw_coef, trait != "CWM_SLA_norot_Axis9")
PCA.SLA_norot_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA_norot_Results_Across_Raw_coef, trait != "CWM_SLA_norot_Axis10")

# Save the spatial scales names
trait.name <- PCA.SLA_norot_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
PCA.SLA_norot_Results_Across_Raw_coef$trait <- factor(PCA.SLA_norot_Results_Across_Raw_coef$trait,
                                                levels= c("CWM_SLA_norot_Axis3", "CWM_SLA_norot_Axis4",
                                                          "CWM_SLA_norot_Axis5", "CWM_SLA_norot_Axis6"))
# The delta R2 are 
# Imputed
trait.labels<-c("Dimension 3 \n Deeper rooting buds \n R² = 0.21",
                "Dimension 4 \n Low root tissue density \n R² = 0.07",
                "Dimension 5 \n \n R² = 0.17",
                "Dimension 6 \n \n R² = 0.19")
#
names(trait.labels) <- levels(PCA.SLA_norot_Results_Across_Raw_coef$trait)
# Remove unused level (Intercept)
PCA.SLA_norot_Results_Across_Raw_coef$term <- as.factor(PCA.SLA_norot_Results_Across_Raw_coef$term)
PCA.SLA_norot_Results_Across_Raw_coef$term <- droplevels(PCA.SLA_norot_Results_Across_Raw_coef$term)

# Give explicit names to Env
levels(PCA.SLA_norot_Results_Across_Raw_coef$term)
levels(PCA.SLA_norot_Results_Across_Raw_coef$term)<- c(
  "Region Hainich", "Region Schorfheide",
  "C:N", "Land-use-intensity", "NH4", "NO3", "N:P",
  "pH", "Phosphorus", "Sand content", "Moisture")
# Reorganize per order of importance (global to proximal)
PCA.SLA_norot_Results_Across_Raw_coef$term <- ordered(
  PCA.SLA_norot_Results_Across_Raw_coef$term,
  levels = c("NH4", "NO3", "Phosphorus", "N:P", "C:N", "pH", "Sand content",
             "Moisture", "Land-use-intensity", "Region Schorfheide", "Region Hainich"))

# Define colors based on significance from the p.values
PCA.SLA_norot_Results_Across_Raw_coef$direction[PCA.SLA_norot_Results_Across_Raw_coef$estimate > 0]   <- "pos"
PCA.SLA_norot_Results_Across_Raw_coef$direction[PCA.SLA_norot_Results_Across_Raw_coef$estimate < 0]   <- "neg"
PCA.SLA_norot_Results_Across_Raw_coef$direction[PCA.SLA_norot_Results_Across_Raw_coef$p.value > 0.05] <- "NS"

# Forest Plot #

png("p3_2_Forest.plot_PCA.SLA_other_axes_norot.png", width = 2500, height = 2000, units = 'px', res = 300)

p3_2_Forest.plot_PCA.SLA_other_axes_norot = ggplot(data=PCA.SLA_norot_Results_Across_Raw_coef,
                                aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.3)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=9,face="bold"), # Title theme
        axis.text.y = element_text(size=9), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=9), # Estimates labels
        axis.title = element_text(size=9,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 9, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 9, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p3_2_Forest.plot_PCA.SLA_other_axes_norot

dev.off()

#### 1.6) PCA Belowground ####

PCA.Belowground_Results_Across_Raw_coef
# Trait raw
trait_raw_names  <-  c("PCA.Belowground_PC1", "PCA.Belowground_PC2")

# filter out phylo traits
PCA.Belowground_Results_Across_Raw_coef <- dplyr::filter(PCA.Belowground_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
PCA.Belowground_Results_Across_Raw_coef <- dplyr::filter(PCA.Belowground_Results_Across_Raw_coef, term != "(Intercept)")
PCA.Belowground_Results_Across_Raw_coef <- data.table(PCA.Belowground_Results_Across_Raw_coef)
PCA.Belowground_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
PCA.Belowground_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- PCA.Belowground_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
PCA.Belowground_Results_Across_Raw_coef$trait <- factor(PCA.Belowground_Results_Across_Raw_coef$trait,
                                                levels= c("PCA.Belowground_PC1", "PCA.Belowground_PC2"))
# The delta R2 are 
# Imputed
trait.labels<-c("Dimension 1 \n Collaboration \n R² = 0.48", "Dimension 2 \n Conservation-regeneration \n R² = 0.58")
#
names(trait.labels) <- levels(PCA.Belowground_Results_Across_Raw_coef$trait)
# Remove unused level (Intercept)
PCA.Belowground_Results_Across_Raw_coef$term <- as.factor(PCA.Belowground_Results_Across_Raw_coef$term)
PCA.Belowground_Results_Across_Raw_coef$term <- droplevels(PCA.Belowground_Results_Across_Raw_coef$term)

# Give explicit names to Env
levels(PCA.Belowground_Results_Across_Raw_coef$term)
levels(PCA.Belowground_Results_Across_Raw_coef$term)<- c(
  "Region Hainich","Region Schorfheide",
  "C:N", "delta-N-15", "Land-use-intensity", "NH4", "NO3", "N:P",
  "pH", "Moisture")
# Reorganize per order of importance (global to proximal)
PCA.Belowground_Results_Across_Raw_coef$term <- ordered(
  PCA.Belowground_Results_Across_Raw_coef$term,
  levels = c("delta-N-15", "NH4", "NO3",  "N:P", "C:N", "pH",
             "Moisture", "Land-use-intensity", "Region Schorfheide", "Region Hainich"))

# Define colors based on significance from the p.values
PCA.Belowground_Results_Across_Raw_coef$direction[PCA.Belowground_Results_Across_Raw_coef$estimate > 0]   <-"pos"
PCA.Belowground_Results_Across_Raw_coef$direction[PCA.Belowground_Results_Across_Raw_coef$estimate < 0]   <-"neg"
PCA.Belowground_Results_Across_Raw_coef$direction[PCA.Belowground_Results_Across_Raw_coef$p.value > 0.05] <-"NS"

# Forest Plot #

png("p4_Forest.plot_PCA.Belowground.png", width = 1500, height = 2000, units = 'px', res = 300)

p4_Forest.plot_PCA.Belowground = ggplot(data=PCA.Belowground_Results_Across_Raw_coef,
                                   aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.3)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=9,face="bold"), # Title theme
        axis.text.y = element_text(size=9), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=9), # Estimates labels
        axis.title = element_text(size=9,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 9, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 9, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p4_Forest.plot_PCA.Belowground

dev.off()

#### 1.7) PCA SLA Guild ####

PCA.SLA.Guild_Results_Across_Raw_coef
# Trait raw
trait_raw_names  <-  c("CWM_SLA_guild_Axis1", "CWM_SLA_guild_Axis2")

# filter out phylo traits
PCA.SLA.Guild_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA.Guild_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
PCA.SLA.Guild_Results_Across_Raw_coef <- dplyr::filter(PCA.SLA.Guild_Results_Across_Raw_coef, term != "(Intercept)")
PCA.SLA.Guild_Results_Across_Raw_coef <- data.table(PCA.SLA.Guild_Results_Across_Raw_coef)
PCA.SLA.Guild_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
PCA.SLA.Guild_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- PCA.SLA.Guild_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
PCA.SLA.Guild_Results_Across_Raw_coef$trait <- factor(PCA.SLA.Guild_Results_Across_Raw_coef$trait,
                                                levels= c("CWM_SLA_guild_Axis1", "CWM_SLA_guild_Axis2"))
# The delta R2 are in FD_models_Env_drivers_traits
trait.labels <- c("Dimension 1 \n Collaboration \n R² = 0.45", "Dimension 2 \n Conservation-regeneration \n R² = 0.52")
#
names(trait.labels) <- levels(PCA.SLA.Guild_Results_Across_Raw_coef$trait)
# Remove unused level (Intercept)
PCA.SLA.Guild_Results_Across_Raw_coef$term <- as.factor(PCA.SLA.Guild_Results_Across_Raw_coef$term)
PCA.SLA.Guild_Results_Across_Raw_coef$term <- droplevels(PCA.SLA.Guild_Results_Across_Raw_coef$term)

# Give explicit names to Env
levels(PCA.SLA.Guild_Results_Across_Raw_coef$term)
levels(PCA.SLA.Guild_Results_Across_Raw_coef$term) <- c(
  "Region Hainich","Region Schorfheide",
  "C:N", "delta-N-15", "Land-use-intensity", "NH4", "NO3", "N:P",
  "pH", "Phosphorus", "Moisture")

# Reorganize per order of importance (global to proximal)
PCA.SLA.Guild_Results_Across_Raw_coef$term <- ordered(
  PCA.SLA.Guild_Results_Across_Raw_coef$term,
  levels = c("delta-N-15", "NH4", "NO3", "Phosphorus", "N:P", "C:N", "pH", "Moisture",
             "Land-use-intensity",
             "Region Schorfheide", "Region Hainich"))

# Define colors based on significance from the p.values
PCA.SLA.Guild_Results_Across_Raw_coef$direction[PCA.SLA.Guild_Results_Across_Raw_coef$estimate > 0]   <-"pos"
PCA.SLA.Guild_Results_Across_Raw_coef$direction[PCA.SLA.Guild_Results_Across_Raw_coef$estimate < 0]   <-"neg"
PCA.SLA.Guild_Results_Across_Raw_coef$direction[PCA.SLA.Guild_Results_Across_Raw_coef$p.value > 0.05] <-"NS"

# Forest Plot #

png("p5_Forest.plot_PCA.SLA.Guild.png", width = 1500, height = 2000, units = 'px', res = 300)

p5_Forest.plot_PCA.SLA.Guild = ggplot(data=PCA.SLA.Guild_Results_Across_Raw_coef,
           aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.3)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=9,face="bold"), # Title theme
        axis.text.y = element_text(size=9), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=9), # Estimates labels
        axis.title = element_text(size=9,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 9, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 9, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p5_Forest.plot_PCA.SLA.Guild

dev.off()

#### 1.8) PCA Belowground Guild ####

PCA.Belowground.Guild_Results_Across_Raw_coef
# Trait raw
trait_raw_names  <-  c("CWM_Belowground_guild_norot_Axis1", "CWM_Belowground_guild_norot_Axis2")

# filter out phylo traits
PCA.Belowground.Guild_Results_Across_Raw_coef <- dplyr::filter(PCA.Belowground.Guild_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
PCA.Belowground.Guild_Results_Across_Raw_coef <- dplyr::filter(PCA.Belowground.Guild_Results_Across_Raw_coef, term != "(Intercept)")
PCA.Belowground.Guild_Results_Across_Raw_coef <- data.table(PCA.Belowground.Guild_Results_Across_Raw_coef)
PCA.Belowground.Guild_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
PCA.Belowground.Guild_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- PCA.Belowground.Guild_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
PCA.Belowground.Guild_Results_Across_Raw_coef$trait <- factor(PCA.Belowground.Guild_Results_Across_Raw_coef$trait,
                                                      levels= c("CWM_Belowground_guild_norot_Axis1",
                                                                "CWM_Belowground_guild_norot_Axis2"))
# The delta R2 are in FD_models_Env_drivers_traits
trait.labels <- c("Dimension 1 \n Do-it-yourself \n R² = 0.40", "Dimension 2 \n Conservation and regeneration \n R² = 0.52")
#
names(trait.labels) <- levels(PCA.Belowground.Guild_Results_Across_Raw_coef$trait)
# Remove unused level (Intercept)
PCA.Belowground.Guild_Results_Across_Raw_coef$term <- as.factor(PCA.Belowground.Guild_Results_Across_Raw_coef$term)
PCA.Belowground.Guild_Results_Across_Raw_coef$term <- droplevels(PCA.Belowground.Guild_Results_Across_Raw_coef$term)

# Give explicit names to Env
levels(PCA.Belowground.Guild_Results_Across_Raw_coef$term)
levels(PCA.Belowground.Guild_Results_Across_Raw_coef$term)<- c(
  "Region Hainich","Region Schorfheide",
  "C:N", "delta-N-15", "Land-use-intensity", "NH4", "NO3", "N:P",
  "pH", "Moisture")

# Reorganize per order of importance (global to proximal)
PCA.Belowground.Guild_Results_Across_Raw_coef$term <- ordered(
  PCA.Belowground.Guild_Results_Across_Raw_coef$term,
  levels = c("delta-N-15", "NH4", "NO3","N:P", "C:N", "pH", "Moisture",
             "Land-use-intensity",
             "Region Schorfheide", "Region Hainich"))

# Define colors based on significance from the p.values
PCA.Belowground.Guild_Results_Across_Raw_coef$direction[PCA.Belowground.Guild_Results_Across_Raw_coef$estimate > 0]   <-"pos"
PCA.Belowground.Guild_Results_Across_Raw_coef$direction[PCA.Belowground.Guild_Results_Across_Raw_coef$estimate < 0]   <-"neg"
PCA.Belowground.Guild_Results_Across_Raw_coef$direction[PCA.Belowground.Guild_Results_Across_Raw_coef$p.value > 0.05] <-"NS"

# Forest Plot #

png("p6_Forest.plot_PCA.Belowground.Guild.png", width = 1500, height = 2000, units = 'px', res = 300)

p6_Forest.plot_PCA.Belowground.Guild = ggplot(data=PCA.Belowground.Guild_Results_Across_Raw_coef,
                                      aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.3)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=9,face="bold"), # Title theme
        axis.text.y = element_text(size=9), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=9), # Estimates labels
        axis.title = element_text(size=9,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 9, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 9, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p6_Forest.plot_PCA.Belowground.Guild

dev.off()

#### 1.9) PCA_p Belowground ####
PCA.p_SLA.Guild_Results_Across_Raw_coef
# Trait raw
trait_raw_names  <-  c("CWM_p_ext_PCA_Axis1", "CWM_p_ext_PCA_Axis2")

# filter out phylo traits
PCA.p_SLA.Guild_Results_Across_Raw_coef <- dplyr::filter(PCA.p_SLA.Guild_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
PCA.p_SLA.Guild_Results_Across_Raw_coef <- dplyr::filter(PCA.p_SLA.Guild_Results_Across_Raw_coef, term != "(Intercept)")
PCA.p_SLA.Guild_Results_Across_Raw_coef <- data.table(PCA.p_SLA.Guild_Results_Across_Raw_coef)
PCA.p_SLA.Guild_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
PCA.p_SLA.Guild_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- PCA.p_SLA.Guild_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
PCA.p_SLA.Guild_Results_Across_Raw_coef$trait <- factor(PCA.p_SLA.Guild_Results_Across_Raw_coef$trait,
                                                levels= c("CWM_p_ext_PCA_Axis1","CWM_p_ext_PCA_Axis2"))
# The delta R2 are 
# Imputed
trait.labels<-c("pCWM PC1 \n R² = 0.36","pCWM PC2 \n R² = 0.45")
#
names(trait.labels) <- levels(PCA.p_SLA.Guild_Results_Across_Raw_coef$trait)

# Remove unused level (Intercept)
PCA.p_SLA.Guild_Results_Across_Raw_coef$term <-as.factor(PCA.p_SLA.Guild_Results_Across_Raw_coef$term)
PCA.p_SLA.Guild_Results_Across_Raw_coef$term <- droplevels(PCA.p_SLA.Guild_Results_Across_Raw_coef$term)

# Give explicit names to Env
levels(PCA.p_SLA.Guild_Results_Across_Raw_coef$term)
levels(PCA.p_SLA.Guild_Results_Across_Raw_coef$term)<- c(
  "Region Hainich", "Region Schorfheide",
  "C:N", "delta-N-15", "Land-use-intensity", "NH4", "NO3", "N:P",
  "pH", "Moisture")
# Reorganize per order of importance (global to proximal)
PCA.p_SLA.Guild_Results_Across_Raw_coef$term <- ordered(
  PCA.p_SLA.Guild_Results_Across_Raw_coef$term,
  levels = c("delta-N-15", "NH4", "NO3", "N:P", "C:N", "pH",
              "Moisture","Land-use-intensity", "Region Schorfheide", "Region Hainich"))

# Define colors based on significance from the p.values
PCA.p_SLA.Guild_Results_Across_Raw_coef$direction[PCA.p_SLA.Guild_Results_Across_Raw_coef$estimate > 0]   <-"pos"
PCA.p_SLA.Guild_Results_Across_Raw_coef$direction[PCA.p_SLA.Guild_Results_Across_Raw_coef$estimate < 0]   <-"neg"
PCA.p_SLA.Guild_Results_Across_Raw_coef$direction[PCA.p_SLA.Guild_Results_Across_Raw_coef$p.value > 0.05] <-"NS"

# Forest Plot #

png("p6_Forest.plot_PCA.p_Belowground.png", width = 1500, height = 2000, units = 'px', res = 300)

p6_Forest.plot_PCA.p_Belowground = ggplot(data=PCA.p_SLA.Guild_Results_Across_Raw_coef,
                                   aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.3)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=9,face="bold"), # Title theme
        axis.text.y = element_text(size=9), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=9), # Estimates labels
        axis.title = element_text(size=9,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 9, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 9, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p6_Forest.plot_PCA.p_Belowground

dev.off()

#### 1.9) Guild ####

Guild_Results_Across_Raw_coef
# Trait raw
trait_raw_names  <-  c("CWM_Guild.forb", "CWM_Guild.grass", "CWM_Guild.herbaceous.legume")

# filter out phylo traits
Guild_Results_Across_Raw_coef <- dplyr::filter(Guild_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
Guild_Results_Across_Raw_coef <- dplyr::filter(Guild_Results_Across_Raw_coef, term != "(Intercept)")
Guild_Results_Across_Raw_coef <- data.table(Guild_Results_Across_Raw_coef)
Guild_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
Guild_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- Guild_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
Guild_Results_Across_Raw_coef$trait <- factor(Guild_Results_Across_Raw_coef$trait,
                                             levels= c("CWM_Guild.forb", "CWM_Guild.grass", "CWM_Guild.herbaceous.legume"))
# The delta R2 are 
# Imputed
trait.labels <- c("Proportion of forbs \n R² = 0.23", "Proportion of Poales \n R² = 0.28", "Proportion of Fabaceae \n R² = 0.27")
#
names(trait.labels) <- levels(Guild_Results_Across_Raw_coef$trait)

# Remove unused level (Intercept)
Guild_Results_Across_Raw_coef$term <- as.factor(Guild_Results_Across_Raw_coef$term)
Guild_Results_Across_Raw_coef$term <- droplevels(Guild_Results_Across_Raw_coef$term)

# Give explicit names to Env
levels(Guild_Results_Across_Raw_coef$term)
levels(Guild_Results_Across_Raw_coef$term)<- c(
  "Region Hainich", "Region Schorfheide",
  "C:N", "delta-N-15", "Land-use-intensity", "NH4", "NO3","N:P",
  "pH", "Sand content", "Moisture")
# Reorganize per order of importance (global to proximal)
Guild_Results_Across_Raw_coef$term <- ordered(
  Guild_Results_Across_Raw_coef$term,
  levels = c("delta-N-15", "NH4", "NO3", "N:P", "C:N", "pH", "Sand content",
             "Moisture", "Land-use-intensity", "Region Schorfheide", "Region Hainich"))

# Define colors based on significance from the p.values
Guild_Results_Across_Raw_coef$direction[Guild_Results_Across_Raw_coef$estimate > 0]   <-"pos"
Guild_Results_Across_Raw_coef$direction[Guild_Results_Across_Raw_coef$estimate < 0]   <-"neg"
Guild_Results_Across_Raw_coef$direction[Guild_Results_Across_Raw_coef$p.value > 0.05] <-"NS"

# Forest Plot #

png("p7_Forest.plot_Guild.png", width = 1500, height = 2000, units = 'px', res = 300)

p7_Forest.plot_Guild = ggplot(data=Guild_Results_Across_Raw_coef,
                                   aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the point for estimates
  geom_pointrange(aes(col=direction), size = 0.3)+
  
  # Replot the points for those covered by the band
  geom_pointrange(aes(col=direction), size = 0.3)+ 
  scale_color_manual(values=c("red","grey42","blue"))+
  
  facet_wrap(~trait, ncol=9, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.3)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title = element_text(size=7,face="bold"), # Title theme
        axis.text.y = element_text(size=7), # Traits labels
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=7), # Estimates labels
        axis.title = element_text(size=7,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 7, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 7, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p7_Forest.plot_Guild

dev.off()

#### MS Fig 2B: CWM_per_slice ####

# Bar colors based on direction
CWM_per_slice$direction[CWM_per_slice$CWM_dif > 0]   <-"pos"
CWM_per_slice$direction[CWM_per_slice$CWM_dif < 0]   <-"neg"

# Subset for each plot we want

# Top left: Specific root length and Bud-bank size
CWM_Top_left <- dplyr::filter(CWM_per_slice, trait %in% c('CWM Specific root length','CWM Bud-bank size') &
                                Slice == "Top_left_conservative_DIY")
# Top right: Root-weight-ratio and Mycorrhizal colonization
CWM_Top_right <- dplyr::filter(CWM_per_slice, trait %in% c("CWM Root weight ratio", "CWM Mycorrhizal colonization") &
                                Slice == "Top_right_conservative_myco")
# Bottom left: Specific leaf area and Branching intensity
CWM_Bottom_left <- dplyr::filter(CWM_per_slice, trait %in% c("CWM Specific leaf area", "CWM Branching intensity") &
                                 Slice == "Bottom_left_acquisitive_DIY")
# Bottom right: Rooting depth 50% and Fine roots diameter
CWM_Bottom_right <- dplyr::filter(CWM_per_slice, trait %in% c("CWM Rooting depth 50%", "CWM Fine roots %N") &
                                   Slice == "Bottom_right_acquisitive_myco")
CWM_list <- list(CWM_Top_left, CWM_Top_right, CWM_Bottom_left, CWM_Bottom_right)
#
# Plot
# Top left

png("p_CWM_Top_left.png", width = 850, height = 850, units = 'px', res = 300)

p_CWM_Top_left <- ggplot(CWM_Top_left, aes(x=trait, y=CWM, fill=direction)) +
     geom_bar(stat="identity")+ scale_fill_manual(values=c("blue","red"))+ labs(y= 'Difference to global CWM')+
     theme_minimal()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
           axis.line = element_line(colour = "black"), axis.title.x = element_blank(),
           axis.text.x = element_text(color="black") , legend.position = "none")+
  scale_x_discrete(labels = c("CWM\nSpecific root length",
                              "CWM\nBud-bank size"))
p_CWM_Top_left
dev.off()

# Top right
png("p_CWM_Top_right.png", width = 850, height = 850, units = 'px', res = 300)

p_CWM_Top_right <- ggplot(CWM_Top_right, aes(x=trait, y=CWM, fill=direction)) +
  geom_bar(stat="identity")+ scale_fill_manual(values=c("blue","red"))+ labs(y= 'Difference to global CWM')+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.title.x = element_blank(),
        axis.text.x = element_text(color="black") , legend.position = "none")+
  scale_x_discrete(labels = c("CWM\nRoot weight ratio",
                              "CWM\nMycorrhizal\ncolonization"))

p_CWM_Top_right
dev.off()

# Bottom left
png("p_CWM_Bottom_left.png", width = 850, height = 850, units = 'px', res = 300)

p_CWM_Bottom_left <- ggplot(CWM_Bottom_left, aes(x=trait, y=CWM, fill=direction)) +
  geom_bar(stat="identity")+ scale_fill_manual(values=c("blue","red"))+ labs(y= 'Difference to global CWM')+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.title.x = element_blank(),
        axis.text.x = element_text(color="black") , legend.position = "none")+
  scale_x_discrete(labels = c("CWM\nSpecific leaf area",
                              "CWM\nBranching intensity"))
p_CWM_Bottom_left
dev.off()
# Bottom right
png("p_CWM_Bottom_right.png", width = 850, height = 850, units = 'px', res = 300)

p_CWM_Bottom_right <- ggplot(CWM_Bottom_right, aes(x=trait, y=CWM, fill=direction)) +
  geom_bar(stat="identity")+ scale_fill_manual(values=c("blue","red"))+ labs(y= 'Difference to global CWM')+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.title.x = element_blank(),
        axis.text.x = element_text(color="black") , legend.position = "none")+
  scale_x_discrete(labels = c("CWM\nRooting depth 50%",
                              "CWM\nFine roots %N"))
p_CWM_Bottom_right
dev.off()


#### 2) Combined Estimates Plots CWM CWV superimposed ####
#### 2.1) Raw CWM CWV ####
#
CWM_Results_Across_Raw_coef$Type <- "CWM"
CWV_Results_Across_Raw_coef$Type <- "CWV"
#
CWM_Results_Across_Raw_coef$trait = gsub(pattern = "CWM_", replacement = "", x = CWM_Results_Across_Raw_coef$trait)
CWV_Results_Across_Raw_coef$trait = gsub(pattern = "CWV_", replacement = "", x = CWV_Results_Across_Raw_coef$trait)
#
CWM_CWV_Results_Across_Raw_coef <- bind_rows(CWM_Results_Across_Raw_coef, CWV_Results_Across_Raw_coef)
#
levels(CWM_CWV_Results_Across_Raw_coef$term)
CWM_CWV_Results_Across_Raw_coef$trait<-as.factor(CWM_CWV_Results_Across_Raw_coef$trait)
levels(CWM_CWV_Results_Across_Raw_coef$trait)
CWM_CWV_Results_Across_Raw_coef$Type<-as.factor(CWM_CWV_Results_Across_Raw_coef$Type)

# Trait raw
levels(CWM_CWV_Results_Across_Raw_coef$trait)
trait_raw_names  <-  c("Klim_BBRsum_est",                              
                       "LEDA_Height_est",                             
                       "LEDA_SLA_est",                                
                       "Morpho_Root.1Order.average.diameter_mean_est",
                       "Morpho_Root.tissue.density.Rose_mean_est",    
                       "Morpho_Root.weight.ratio_mean",           
                       "Morpho_Specific.root.length_mean_est",        
                       "RDepth_Maximum.rooting.depth_mean_est",       
                       "Seed.weight")

# filter out phylo traits
CWM_CWV_Results_Across_Raw_coef <- dplyr::filter(CWM_CWV_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
CWM_CWV_Results_Across_Raw_coef <- dplyr::filter(CWM_CWV_Results_Across_Raw_coef, term != "(Intercept)")
CWM_CWV_Results_Across_Raw_coef <- data.table(CWM_CWV_Results_Across_Raw_coef)
CWM_CWV_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
CWM_CWV_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- CWM_CWV_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
CWM_CWV_Results_Across_Raw_coef$trait <- factor(CWM_CWV_Results_Across_Raw_coef$trait,
     levels= c("LEDA_SLA_est", "LEDA_Height_est", "Seed.weight",
               "Morpho_Specific.root.length_mean_est", "Morpho_Root.1Order.average.diameter_mean_est",
               "Morpho_Root.tissue.density.Rose_mean_est", "RDepth_Maximum.rooting.depth_mean_est",
               "Morpho_Root.weight.ratio_mean", "Klim_BBRsum_est"))

# The delta R2 are
trait.labels<-c("Specific leaf area \n R² CWM 0.47 \n R² CWV 0.42",    
                "Height \n 0.22 \n 0.14", "Seed weight \n 0.44 \n 0.30",
                "Specific root length \n 0.31 \n 0.24",  "Fine roots diameter \n 0.50 \n 0.15",
                "Root tissue density \n 0.19 \n 0.25",   "Maximum rooting depth \n 0.66 \n 0.24",
                "Root weight ratio \n 0.58 \n 0.52",     "Bud bank size \n 0.44 \n 0.38")

names(trait.labels) <- levels(CWM_CWV_Results_Across_Raw_coef$trait)

# Give explicit names to traits
CWM_CWV_Results_Across_Raw_coef$term <- as.factor(CWM_CWV_Results_Across_Raw_coef$term)
levels(CWM_CWV_Results_Across_Raw_coef$term)
levels(CWM_CWV_Results_Across_Raw_coef$term) <- c("Intercept",
                                                 "Region Hainich","Region Schorfheide",
                                                 #"Aspect Eastness", "Aspect Northness", "Bacteria",
                                                 "Bulk density", "C/N ratio",
                                                 "d15N soil","Fertilization",
                                                 "Grazing","Mowing","NH4","NO3","N/P ratio",
                                                 "pH", "Silt", "Slope",
                                                 "Soil moisture", "Temperature",
                                                 "Soil Phosphorus" # "Topographic Wetness Index"
                                                 )
# Remove unused level (Intercept)
term <- droplevels(CWM_CWV_Results_Across_Raw_coef$term)
CWM_CWV_Results_Across_Raw_coef$term <- droplevels(CWM_CWV_Results_Across_Raw_coef$term)

# Reorganize per order of importance (most response at top, least at bottom)
levels(CWM_CWV_Results_Across_Raw_coef$term)
CWM_CWV_Results_Across_Raw_coef$term <- factor(term, levels = c(
  "Region Schorfheide",
  "Region Hainich",
  "Temperature",
  "Soil moisture",
  #"Topographic Wetness Index",
  #"Aspect Eastness",
  #"Aspect Northness",
  "Slope",
  "Bacteria",
  "C/N ratio",
  "N/P ratio",
  "Soil Phosphorus",
  "d15N soil",
  "NH4",
  "NO3",
  "pH",
  "Silt",
  "Bulk density",
  "Grazing",
  "Mowing",
  "Fertilization"))
# Define colors based on significance from the p.values
CWM_CWV_Results_Across_Raw_coef$direction[CWM_CWV_Results_Across_Raw_coef$estimate > 0]   <-"pos"
CWM_CWV_Results_Across_Raw_coef$direction[CWM_CWV_Results_Across_Raw_coef$estimate < 0]   <-"neg"
CWM_CWV_Results_Across_Raw_coef$direction[CWM_CWV_Results_Across_Raw_coef$p.value > 0.05] <-"NS"

CWM_CWV_Results_Across_Raw_coef$agg_direction<-paste(CWM_CWV_Results_Across_Raw_coef$Type,
                                                     CWM_CWV_Results_Across_Raw_coef$direction,
                                                     sep="_")
CWM_CWV_Results_Across_Raw_coef$agg_direction <- as.factor(CWM_CWV_Results_Across_Raw_coef$agg_direction)
levels(CWM_CWV_Results_Across_Raw_coef$agg_direction)

png("p6_Forest_CWM_CWV_Raw.png", width = 4000, height = 2500, units = 'px', res = 300)

p6_forest_frequentist = ggplot(data=CWM_CWV_Results_Across_Raw_coef,
                               aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the points for estimates
  geom_pointrange(aes(col=agg_direction, shape= Type), size = 0.6)+ 
  scale_color_manual(values=c("red","grey42","blue",
                              "orangered4","grey42","royalblue4"))+
  scale_shape_manual(values=c(1, 6))+
  
  
  facet_wrap(~trait, ncol=10, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.2)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=agg_direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title=element_text(size=6,face="bold"), # Title theme
        axis.text.y=element_text(size=6), # Traits labels
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=4), # Estimates labels
        axis.title=element_text(size=6,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 6, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 6, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p6_forest_frequentist

dev.off()


#### 2.2) Phylo CWM CWV ####
CWM_Results_Across_Raw_coef  <- read.table("CWM_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
CWV_Results_Across_Raw_coef  <- read.table("CWV_Results_Across_Raw_coef.txt",  sep = "\t", na.strings = "NA", h = T, dec = ",")
#
CWM_Results_Across_Raw_coef$Type <- "CWM"
CWV_Results_Across_Raw_coef$Type <- "CWV"
#
CWM_Results_Across_Raw_coef$trait = gsub(pattern = "CWM_p_", replacement = "", x = CWM_Results_Across_Raw_coef$trait)
CWV_Results_Across_Raw_coef$trait = gsub(pattern = "CWV_p_", replacement = "", x = CWV_Results_Across_Raw_coef$trait)
#
CWM_CWV_Results_Across_Raw_coef <- bind_rows(CWM_Results_Across_Raw_coef, CWV_Results_Across_Raw_coef)
#
levels(CWM_CWV_Results_Across_Raw_coef$term)
CWM_CWV_Results_Across_Raw_coef$trait<-as.factor(CWM_CWV_Results_Across_Raw_coef$trait)
levels(CWM_CWV_Results_Across_Raw_coef$trait)
CWM_CWV_Results_Across_Raw_coef$Type<-as.factor(CWM_CWV_Results_Across_Raw_coef$Type)

# Trait raw
levels(CWM_CWV_Results_Across_Raw_coef$trait)
trait_raw_names  <-  c("Klim_BBRsum_est",                              
                       "LEDA_Height_est",                             
                       "LEDA_SLA_est",                                
                       "Morpho_Root.1Order.average.diameter_mean_est",
                       "Morpho_Root.tissue.density.Rose_mean_est",    
                       "Morpho_Root.weight.ratio_mean",           
                       "Morpho_Specific.root.length_mean_est",        
                       "RDepth_Maximum.rooting.depth_mean_est",       
                       "Seed.weight")

# filter out phylo traits
CWM_CWV_Results_Across_Raw_coef <- dplyr::filter(CWM_CWV_Results_Across_Raw_coef, trait %in% trait_raw_names)

# Remove the intercept and calculate the bounds for the x scale
CWM_CWV_Results_Across_Raw_coef <- dplyr::filter(CWM_CWV_Results_Across_Raw_coef, term != "(Intercept)")
CWM_CWV_Results_Across_Raw_coef <- data.table(CWM_CWV_Results_Across_Raw_coef)
CWM_CWV_Results_Across_Raw_coef[,y_min := -(max(abs(estimate))+max(std.error)), by = term]
CWM_CWV_Results_Across_Raw_coef[,y_max := max(abs(estimate))+max(std.error), by = term]

# Save the spatial scales names
trait.name <- CWM_CWV_Results_Across_Raw_coef$trait
trait.name <- unique(trait.name)
trait.name

# Rearrange the order of spatial scales from global to local from left to right
CWM_CWV_Results_Across_Raw_coef$trait <- factor(CWM_CWV_Results_Across_Raw_coef$trait,
             levels= c("LEDA_SLA_est","LEDA_Height_est","Seed.weight",
                       "Morpho_Specific.root.length_mean_est","Morpho_Root.1Order.average.diameter_mean_est",
                       "Morpho_Root.tissue.density.Rose_mean_est","RDepth_Maximum.rooting.depth_mean_est",
                       "Morpho_Root.weight.ratio_mean","Klim_BBRsum_est"))

# The delta R2 are
trait.labels<-c("Specific leaf area \n R² pCWM 0.43 \n R² pCWV 0.39",    
                "Height \n 0.25 \n 0.18", "Seed weight \n 0.07 \n 0.33",
                "Specific root length \n 0.32 \n 0.22",  "Fine roots diameter \n 0.15 \n 0.25",
                "Root tissue density \n 0.22 \n 0.25",   "Maximum rooting depth \n 0.65 \n 0.30",
                "Root weight ratio \n 0.37 \n 0.31",     "Bud bank size \n 0.43 \n 0.29")

names(trait.labels) <- levels(CWM_CWV_Results_Across_Raw_coef$trait)

# Give explicit names to traits
CWM_CWV_Results_Across_Raw_coef$term<-as.factor(CWM_CWV_Results_Across_Raw_coef$term)
levels(CWM_CWV_Results_Across_Raw_coef$term)
levels(CWM_CWV_Results_Across_Raw_coef$term) <- c("Intercept",
                                                  "Region Hainich","Region Schorfheide",
                                                  #"Aspect Eastness", "Aspect Northness", "Bacteria",
                                                  "Bulk density", "C/N ratio",
                                                  "d15N soil","Fertilization",
                                                  "Grazing","Mowing","NH4","NO3","N/P ratio",
                                                  "pH", "Silt", "Slope",
                                                  "Soil moisture", "Temperature",
                                                  "Soil Phosphorus" #, "Topographic Wetness Index"
                                                  )

# Remove unused level (Intercept)
term <- droplevels(CWM_CWV_Results_Across_Raw_coef$term)
CWM_CWV_Results_Across_Raw_coef$term <- droplevels(CWM_CWV_Results_Across_Raw_coef$term)

# Reorganize per order of importance (most response at top, least at bottom)
levels(CWM_CWV_Results_Across_Raw_coef$term)
CWM_CWV_Results_Across_Raw_coef$term <- factor(term, levels = c(
  "Region Schorfheide",
  "Region Hainich",
  "Temperature",
  "Soil moisture",
  #"Topographic Wetness Index",
  #"Aspect Eastness",
  #"Aspect Northness",
  "Slope",
  #"Bacteria",
  "C/N ratio",
  "N/P ratio",
  "Soil Phosphorus",
  "d15N soil",
  "NH4",
  "NO3",
  "pH",
  "Silt",
  "Bulk density",
  "Grazing",
  "Mowing",
  "Fertilization"))
# Define colors based on significance from the p.values
CWM_CWV_Results_Across_Raw_coef$direction[CWM_CWV_Results_Across_Raw_coef$estimate > 0]   <- "pos"
CWM_CWV_Results_Across_Raw_coef$direction[CWM_CWV_Results_Across_Raw_coef$estimate < 0]   <- "neg"
CWM_CWV_Results_Across_Raw_coef$direction[CWM_CWV_Results_Across_Raw_coef$p.value > 0.05] <- "NS"

CWM_CWV_Results_Across_Raw_coef$agg_direction<-paste(CWM_CWV_Results_Across_Raw_coef$Type,
                                                     CWM_CWV_Results_Across_Raw_coef$direction,
                                                     sep="_")
CWM_CWV_Results_Across_Raw_coef$agg_direction <- as.factor(CWM_CWV_Results_Across_Raw_coef$agg_direction)
levels(CWM_CWV_Results_Across_Raw_coef$agg_direction)

png("p7_Forest_CWM_CWV_Phylo.png", width = 4000, height = 2500, units = 'px', res = 300)

p7_forest_frequentist = ggplot(data=CWM_CWV_Results_Across_Raw_coef,
                               aes(x = term,y = estimate, ymin= estimate - std.error, ymax= estimate + std.error))+
  
  # Plot the points for estimates
  geom_pointrange(aes(col=agg_direction, shape= Type), size = 0.6)+ 
  scale_color_manual(values=c("red","grey42","blue",
                              "orangered4","grey42","royalblue4"))+
  scale_shape_manual(values=c(1, 6))+
  
  
  facet_wrap(~trait, ncol=10, scales="free_x", labeller = labeller(trait=trait.labels)) +
  geom_hline(yintercept = 0, linetype=2, size = 0.2)+ # Vertical line at 0
  xlab('')+ ylab("")+
  #ggtitle(paste("Regression coefficients of predictor traits explaining species success across spatial scales"))+
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error,col=agg_direction),width=0.5,cex=0.4)+ 
  
  theme(plot.title=element_text(size=6,face="bold"), # Title theme
        axis.text.y=element_text(size=6), # Traits labels
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=4), # Estimates labels
        axis.title=element_text(size=6,face="bold"), # Title size
        legend.position = "none",
        strip.text.y = element_text(size = 6, hjust=0,vjust = 1,angle=180,face="bold"),
        strip.text.x = element_text(size = 6, colour = "black"), # Facets labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.border = element_rect(fill = NA, colour = "black"))+
  
  coord_flip()
p7_forest_frequentist

dev.off()

#### 3.1) Scatterplots with prediction - Raw ####

# !!! Run FD_Env_drivers_traits_models first to get the models objects and data frames !!!
df.list <- split(FD_new_models.melt, FD_new_models.melt$trait)

# Keep the retained predictors
columns_names<-c(    "Slope", "Aspect_Northness", "TWI",
                     "Ta_10","SM_10",
                     "pH", "bulk_density_mean",
                     "G_std", "M_std", "F_std", 
                     "CN_ratio", "Total_P_v", "NP_ratio", "Silt_v",
                     "NH4_mean_v", "NO3_mean_v",
                     "bacteria_total_mean_v",
                     "d15N_soil_mean_v")

# Rename the retained predictors
new_names<-c("Slope","Aspect","Topographic Wetness Index","Temperature","Soil moisture",
             "pH","Bulk density",
             "Grazing","Mowing","Fertilization",
             "CN ratio", "Soil Phosphorus","NP ratio","Silt",
             "NH4","NO3",
             "Bacteria","d15N soil")
#
names(new_names) <- columns_names

# # # # #

Plot_My_Raw <- function(df, model, my_lab=NULL){
  
  for(i in columns_names){
    
    new.data <- df
    new.tomerge <- new.data # Additional rows for prediction : duplicate the dataset to permit poly predictions
    drop<-c("value","Exploratory",i) # Turn the other predictors to 0 (they are scaled centered)
    new.data[,!(names(new.data) %in% drop)]<-0

    # Merge
    new.data2<-new.data
    new.data2$predicted <- predict(CWM_Results_Across_Raw[[paste0(model)]], newdata = new.data2 ,type = "response")
    new.data2 <- new.data2[1:nrow(df),]
    
    p_multi_pred <- ggplot(new.data2, aes_string(x = i, y = "value" , color= "Exploratory") ) +
      geom_point(aes(color= Exploratory)) +
      geom_line(aes(y = predicted, color=Exploratory), position = 'jitter', size = 1)+
      xlab(paste(new_names[i]))+
      ylab(paste0(my_lab))+
      theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
            text = element_text(size=12),legend.position = "none")
    put(p_multi_pred)
  }
  
  magic_result()
  
}

# Run function
#
names(df.list)
#
require(magicfor)
magic_for(silent = TRUE)

# Bud bank size
model = "CWM_Klim_BBsum_est"
my_lab = "Bud bank size"
df = df.list[["CWM_Klim_BBsum_est"]]
plot_CWM_Klim_BBsum_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_Klim_BBsum_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_Klim_BBsum_est[["p_multi_pred"]], ncol=4))
dev.off()

# Height
model = "CWM_LEDA_Height_est"
my_lab = "Height"
df = df.list[["CWM_LEDA_Height_est"]]
plot_CWM_LEDA_Height_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_LEDA_Height_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_LEDA_Height_est[["p_multi_pred"]], ncol=4))
dev.off()

# SLA
model = "CWM_LEDA_SLA_est"
my_lab = "Specific leaf area"
df = df.list[["CWM_LEDA_SLA_est"]]
plot_CWM_LEDA_SLA_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_LEDA_SLA_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_LEDA_SLA_est[["p_multi_pred"]], ncol=4))
dev.off()

# Fine roots diameter
model = "CWM_Morpho_Root.1Order.average.diameter_mean_est"
my_lab = "Fine roots diameter"
df = df.list[["CWM_Morpho_Root.1Order.average.diameter_mean_est"]]
plot_CWM_Morpho_Root.1Order.average.diameter_mean_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_Morpho_Root.1Order.average.diameter_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_Morpho_Root.1Order.average.diameter_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Root tissue density
model = "CWM_Morpho_Root.tissue.density.Rose_mean_est"
my_lab = "Root tissue density"
df = df.list[["CWM_Morpho_Root.tissue.density.Rose_mean_est"]]
plot_CWM_Morpho_Root.tissue.density.Rose_mean_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_Morpho_Root.tissue.density.Rose_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_Morpho_Root.tissue.density.Rose_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Root weight ratio
model = "CWM_Morpho_Root.weight.ratio_mean_est"
my_lab = "Root weight ratio"
df = df.list[["CWM_Morpho_Root.weight.ratio_mean_est"]]
plot_CWM_Morpho_Root.weight.ratio_mean_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_Morpho_Root.weight.ratio_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_Morpho_Root.weight.ratio_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Specific root length
model = "CWM_Morpho_Specific.root.length_mean_est"
my_lab = "Specific root length"
df = df.list[["CWM_Morpho_Specific.root.length_mean_est"]]
plot_CWM_Morpho_Specific.root.length_mean_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_Morpho_Specific.root.length_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_Morpho_Specific.root.length_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Maximum rooting depth
model = "CWM_RDepth_Maximum.rooting.depth_mean_est"
my_lab = "Maximum rooting depth"
df = df.list[["CWM_RDepth_Maximum.rooting.depth_mean_est"]]
plot_CWM_RDepth_Maximum.rooting.depth_mean_est <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_RDepth_Maximum.rooting.depth_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_RDepth_Maximum.rooting.depth_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Seed weight
model = "CWM_Seed.weight"
my_lab = "Seed weight"
df = df.list[["CWM_Seed.weight"]]
plot_CWM_Seed.weight <- Plot_My_Raw(df, model, my_lab)
png("p_scatterRaw_CWM_Seed.weight.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(plot_CWM_Seed.weight[["p_multi_pred"]], ncol=4))
dev.off()
#

#### 3.2) Scatterplots with prediction - Poly ####

# !!! Run FD_Env_drivers_traits_models first to get the models objects and data frames !!!
df.list <- split(FD_new_models.melt, FD_new_models.melt$trait)

# Keep the retained predictors
columns_names<-c("Slope","Aspect","pH","G_std","M_std","F_std","bulk_density_mean",
                 "CN_ratio_mean","Total_P","CS_ratio","PS_ratio","Silt","NH4_mean","NO3_mean",
                 "bacteria_total_mean","fungi_mean","d15N_soil_mean","TWI")

# Rename the retained predictors
new_names<-c("Slope","Aspect","pH","Grazing","Mowing","Fertilization","Bulk density",
             "CN ratio", "Soil Phosphorus","CS ratio","PS ratio","Silt","NH4","NO3",
             "Bacteria","Fungi","d15N soil","TWI")
#
names(new_names) <- columns_names

# # # # #

Plot_My_Poly <- function(df, model, my_lab=NULL){
  
  for(i in columns_names){
    
    new.data <- df
    new.tomerge <- new.data # Additional rows for prediction : duplicate the dataset to permit poly predictions
    drop<-c("value","Exploratory",i) # Turn the other predictors to 0 (they are scaled centered)
    new.data[,!(names(new.data) %in% drop)]<-0
    
    # Merge
    new.data2<-rbind(new.data, new.tomerge)
    new.data2$predicted <- predict(Results_Across_Poly[[paste0(model)]], newdata = new.data2 ,type = "response")
    new.data2 <- new.data2[1:nrow(df),]
    
    p_multi_pred <- ggplot(new.data2, aes_string(x = i, y = "value") ) +
      geom_point(aes(colour= Exploratory)) +
      geom_smooth(aes(y = predicted), method = 'loess',  size = 1, se = FALSE, colour="black")+
      xlab(paste(new_names[i]))+
      ylab(paste0(my_lab))+
      theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
            text = element_text(size=12),legend.position = "none")
    put(p_multi_pred)
  }
  
  magic_result()
  
}

# Run function
#
names(df.list)
#
require(magicfor)
magic_for(silent = TRUE)

# Bud bank size
model = "CWM_Klim_BBsum_est"
my_lab = "Bud bank size"
df = df.list[["CWM_Klim_BBsum_est"]]
p_CWM_Klim_BBsum_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_Klim_BBsum_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_Klim_BBsum_est[["p_multi_pred"]], ncol=4))
dev.off()

# Height
model = "CWM_LEDA_Height_est"
my_lab = "Height"
df = df.list[["CWM_LEDA_Height_est"]]
p_CWM_LEDA_Height_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_LEDA_Height_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_LEDA_Height_est[["p_multi_pred"]], ncol=4))
dev.off()

# SLA
model = "CWM_LEDA_SLA_est"
my_lab = "Specific leaf area"
df = df.list[["CWM_LEDA_SLA_est"]]
p_CWM_LEDA_SLA_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_LEDA_SLA_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_LEDA_SLA_est[["p_multi_pred"]], ncol=4))
dev.off()

# Fine roots diameter
model = "CWM_Morpho_Root.1Order.average.diameter_mean_est"
my_lab = "Fine roots diameter"
df = df.list[["CWM_Morpho_Root.1Order.average.diameter_mean_est"]]
p_CWM_Morpho_Root.1Order.average.diameter_mean_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_Morpho_Root.1Order.average.diameter_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_Morpho_Root.1Order.average.diameter_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Root tissue density
model = "CWM_Morpho_Root.tissue.density.Rose_mean_est"
my_lab = "Root tissue density"
df = df.list[["CWM_Morpho_Root.tissue.density.Rose_mean_est"]]
p_CWM_Morpho_Root.tissue.density.Rose_mean_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_Morpho_Root.tissue.density.Rose_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_Morpho_Root.tissue.density.Rose_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Root weight ratio
model = "CWM_Morpho_Root.weight.ratio_mean_est"
my_lab = "Root weight ratio"
df = df.list[["CWM_Morpho_Root.weight.ratio_mean_est"]]
p_CWM_Morpho_Root.weight.ratio_mean_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_Morpho_Root.weight.ratio_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_Morpho_Root.weight.ratio_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Specific root length
model = "CWM_Morpho_Specific.root.length_mean_est"
my_lab = "Specific root length"
df = df.list[["CWM_Morpho_Specific.root.length_mean_est"]]
p_CWM_Morpho_Specific.root.length_mean_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_Morpho_Specific.root.length_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_Morpho_Specific.root.length_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Maximum rooting depth
model = "CWM_RDepth_Maximum.rooting.depth_mean_est"
my_lab = "Maximum rooting depth"
df = df.list[["CWM_RDepth_Maximum.rooting.depth_mean_est"]]
p_CWM_RDepth_Maximum.rooting.depth_mean_est <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_RDepth_Maximum.rooting.depth_mean_est.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_RDepth_Maximum.rooting.depth_mean_est[["p_multi_pred"]], ncol=4))
dev.off()

# Seed weight
model = "CWM_Seed.weight"
my_lab = "Seed weight"
df = df.list[["CWM_Seed.weight"]]
p_CWM_Seed.weight <- Plot_My_Poly(df, model, my_lab)
png("p_scatterPoly_CWM_Seed.weight.png", width = 5000, height = 2500, units = 'px', res = 300)
do.call("grid.arrange", c(p_CWM_Seed.weight[["p_multi_pred"]], ncol=4))
dev.off()
#

#### 4) Relative Importance of predictors (relaimpo) barplots ####

#### 4.1) CWM ralaimpo barplot ####

#View(CWM_relaimpo_sweep_result)

# Filter out phylo coef
CWM_no_phylo_relaimpo_sweep_result <- dplyr::filter(CWM_relaimpo_sweep_result, !grepl("CWM_p_", Trait))
CWM_no_phylo_relaimpo_sweep_result$Env<-as.factor(CWM_no_phylo_relaimpo_sweep_result$Env)
#
names(CWM_no_phylo_relaimpo_sweep_result)
levels(CWM_no_phylo_relaimpo_sweep_result$Env)
levels(CWM_no_phylo_relaimpo_sweep_result$Env)<- c("Exploratory",
                           "Aspect Eastness","Aspect Northness",
                           "Bacterial abundance", "Bulk density",
                           "CN ratio", "d15N soil",
                           "Fertilization","Grazing","Mowing",
                           "NH4","NO3","NP ratio","pH","Silt",
                           "Slope","Soil moisture","Temperature",
                           "Soil phosphorus", "Topographical Wetness Index")

# Sum across traits
CWM_no_phylo_relaimpo_sweep_result_sum <- CWM_no_phylo_relaimpo_sweep_result %>% group_by(Env) %>% summarize_at("relaimpo.result.lmg", list(~sum(.)))
CWM_no_phylo_relaimpo_sweep_result_sum <- CWM_no_phylo_relaimpo_sweep_result_sum[order(CWM_no_phylo_relaimpo_sweep_result_sum$relaimpo.result.lmg, decreasing = TRUE),]  
CWM_relaimpo<-CWM_no_phylo_relaimpo_sweep_result_sum
#
png("p_CWM_relaimpo.png", width = 2000, height = 1200, units = 'px', res = 300)

ggplot(CWM_relaimpo, aes(x=reorder(Env, relaimpo.result.lmg), y=relaimpo.result.lmg)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Environmental driver")+
  ylab("R² partitioned by averaging over orders (lmg) for CWM")+
  theme(
    strip.text.x = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 6, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 6, hjust = 0, colour = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))+
  coord_flip()


dev.off()

#### 4.1.2) CWM relaimpo predictors grouped by category ####
levels(CWM_no_phylo_relaimpo_sweep_result_sum$Env)
# Create category
CWM_no_phylo_relaimpo_sweep_result_sum$Cat_Env <- CWM_no_phylo_relaimpo_sweep_result_sum$Env
CWM_no_phylo_relaimpo_sweep_result_sum$Cat_Env<-plyr::revalue(CWM_no_phylo_relaimpo_sweep_result_sum$Cat_Env,
               c(# Geography
                 "Aspect Eastness"="Geography",
                 "Aspect Northness"="Geography",
                 "Slope"="Geography",
                 "Temperature"="Geography",
                 "Topographical Wetness Index"="Geography",
                 # Soil
                 "Soil moisture"="Soil",
                 "Bacterial abundance"="Soil",
                 "Bulk density"="Soil",
                 "CN ratio"="Soil",
                 "d15N soil"="Soil",
                 "NH4"="Soil",
                 "NO3"="Soil",
                 "NP ratio"="Soil",
                 "pH"="Soil",
                 "Silt"="Soil",
                 "Soil phosphorus"="Soil",
                 # Land-use
                 "Mowing"="Land use",
                 "Grazing"="Land use",
                 "Fertilization"="Land use"
                 ))


#
png("p_CWM_Cat_relaimpo.png", width = 2000, height = 1200, units = 'px', res = 300)

ggplot(CWM_no_phylo_relaimpo_sweep_result_sum, aes(x=reorder(Cat_Env, relaimpo.result.lmg), y=relaimpo.result.lmg)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Environmental driver")+
  ylab("R² partitioned by averaging over orders (lmg) for CWM")+
  theme(
    strip.text.x = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 6, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 6, hjust = 0, colour = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))+
  coord_flip()


dev.off()

#### 4.2) CWV ralaimpo barplot ####

#View(CWV_relaimpo_sweep_result)

# Filter out phylo coef
CWV_no_phylo_relaimpo_sweep_result <- dplyr::filter(CWV_relaimpo_sweep_result, !grepl("CWV_p_", Trait))
CWV_no_phylo_relaimpo_sweep_result$Env<-as.factor(CWV_no_phylo_relaimpo_sweep_result$Env)
#
names(CWV_no_phylo_relaimpo_sweep_result)
levels(CWV_no_phylo_relaimpo_sweep_result$Env)
levels(CWV_no_phylo_relaimpo_sweep_result$Env) <- c("Exploratory",
                                                    "Aspect Eastness","Aspect Northness",
                                                    "Bacterial abundance", "Bulk density",
                                                    "CN ratio", "d15N soil",
                                                    "Fertilization","Grazing","Mowing",
                                                    "NH4","NO3","NP ratio","pH","Silt",
                                                    "Slope","Soil moisture","Temperature",
                                                    "Soil phosphorus", "Topographical Wetness Index")

# Sum across traits
CWV_no_phylo_relaimpo_sweep_result_sum <- CWV_no_phylo_relaimpo_sweep_result %>% group_by(Env) %>% summarize_at("relaimpo.result.lmg", list(~sum(.)))
CWV_no_phylo_relaimpo_sweep_result_sum <- CWV_no_phylo_relaimpo_sweep_result_sum[order(CWV_no_phylo_relaimpo_sweep_result_sum$relaimpo.result.lmg, decreasing = TRUE),]  
CWV_relaimpo<-CWV_no_phylo_relaimpo_sweep_result_sum
#
png("p_CWV_relaimpo.png", width = 2000, height = 1200, units = 'px', res = 300)

ggplot(CWV_relaimpo, aes(x=reorder(Env, relaimpo.result.lmg), y=relaimpo.result.lmg)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Environmental driver")+
  ylab("R² partitioned by averaging over orders (lmg) for CWV")+
  theme(
    strip.text.x = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 6, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 6, hjust = 0, colour = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))+
  coord_flip()

dev.off()

#### 4.2.2) CWV relaimpo predictors grouped by category ####
levels(CWV_no_phylo_relaimpo_sweep_result_sum$Env)
# Create category
CWV_no_phylo_relaimpo_sweep_result_sum$Cat_Env <- CWV_no_phylo_relaimpo_sweep_result_sum$Env
CWV_no_phylo_relaimpo_sweep_result_sum$Cat_Env<-plyr::revalue(CWV_no_phylo_relaimpo_sweep_result_sum$Cat_Env,
                                                              c(# Geography
                                                                "Aspect Eastness"="Geography",
                                                                "Aspect Northness"="Geography",
                                                                "Slope"="Geography",
                                                                "Temperature"="Geography",
                                                                "Topographical Wetness Index"="Geography",
                                                                # Soil
                                                                "Soil moisture"="Soil",
                                                                "Bacterial abundance"="Soil",
                                                                "Bulk density"="Soil",
                                                                "CN ratio"="Soil",
                                                                "d15N soil"="Soil",
                                                                "NH4"="Soil",
                                                                "NO3"="Soil",
                                                                "NP ratio"="Soil",
                                                                "pH"="Soil",
                                                                "Silt"="Soil",
                                                                "Soil phosphorus"="Soil",
                                                                # Land-use
                                                                "Mowing"="Land use",
                                                                "Grazing"="Land use",
                                                                "Fertilization"="Land use"
                                                              ))


#
png("p_CWV_Cat_relaimpo.png", width = 2000, height = 1200, units = 'px', res = 300)

ggplot(CWV_no_phylo_relaimpo_sweep_result_sum, aes(x=reorder(Cat_Env, relaimpo.result.lmg), y=relaimpo.result.lmg)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Environmental driver")+
  ylab("R² partitioned by averaging over orders (lmg) for CWV")+
  theme(
    strip.text.x = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 6, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 6, hjust = 0, colour = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))+
  coord_flip()


dev.off()

#### 4.3) CWM CWV relaimpo stacked barplot ####
CWM_relaimpo_sweep_result$Type <- "CWM"
CWV_relaimpo_sweep_result$Type <- "CWV"

CWM_CWV_relaimpo <- bind_rows(CWM_relaimpo_sweep_result, CWV_relaimpo_sweep_result)
#

# Filter out phylo coef
CWM_CWV_relaimpo <- dplyr::filter(CWM_CWV_relaimpo, !grepl("_p_", Trait))
CWM_CWV_relaimpo$Env<-as.factor(CWM_CWV_relaimpo$Env)
#
names(CWM_CWV_relaimpo)
levels(CWM_CWV_relaimpo$Env)
levels(CWM_CWV_relaimpo$Env) <- c("Exploratory",
                                  "Aspect Eastness", "Aspect Northness",
                                  "Bacterial abundance", "Bulk density",
                                  "CN ratio", "d15N soil",
                                  "Fertilization","Grazing","Mowing",
                                  "NH4","NO3","NP ratio","pH","Silt",
                                  "Slope","Soil moisture","Temperature",
                                  "Soil phosphorus", "Topographical Wetness Index")

# Sum across traits
CWM_CWV_relaimpo_sum <- CWM_CWV_relaimpo %>% group_by(Env, Type) %>% summarize_at("relaimpo.result.lmg", list(~sum(.)))

#
png("p_CWM_CWV_relaimpo.png", width = 2000, height = 1200, units = 'px', res = 300)

ggplot(CWM_CWV_relaimpo_sum, aes(x=reorder(Env, relaimpo.result.lmg), y=relaimpo.result.lmg, fill=Type)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Environmental driver")+
  ylab("R² partitioned by averaging over orders (lmg) for CWV")+
  theme(
    strip.text.x = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 6, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 6, hjust = 0, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    plot.title = element_blank(),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))+
  coord_flip()

dev.off()

#### 4.3.2) CWM CWV relaimpo stacked barplot grouped by category ####
levels(CWM_CWV_relaimpo_sum$Env)
# Create category
CWM_CWV_relaimpo_sum$Cat_Env <- CWM_CWV_relaimpo_sum$Env
CWM_CWV_relaimpo_sum$Cat_Env<-plyr::revalue(CWM_CWV_relaimpo_sum$Cat_Env,
                                                              c(# Geography
                                                                "Aspect Eastness"="Geography",
                                                                "Aspect Northness"="Geography",
                                                                "Slope"="Geography",
                                                                "Temperature"="Geography",
                                                                "Topographical Wetness Index"="Geography",
                                                                # Soil
                                                                "Soil moisture"="Soil",
                                                                "Bacterial abundance"="Soil",
                                                                "Bulk density"="Soil",
                                                                "CN ratio"="Soil",
                                                                "d15N soil"="Soil",
                                                                "NH4"="Soil",
                                                                "NO3"="Soil",
                                                                "NP ratio"="Soil",
                                                                "pH"="Soil",
                                                                "Silt"="Soil",
                                                                "Soil phosphorus"="Soil",
                                                                # Land-use
                                                                "Mowing"="Land use",
                                                                "Grazing"="Land use",
                                                                "Fertilization"="Land use"
                                                              ))


#
png("p_CWM_CWV_Cat_relaimpo.png", width = 2000, height = 1200, units = 'px', res = 300)

ggplot(CWM_CWV_relaimpo_sum, aes(x=reorder(Cat_Env, relaimpo.result.lmg), y=relaimpo.result.lmg, fill=Type)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Environmental driver")+
  ylab("R² partitioned by averaging over orders (lmg) for CWV")+
  theme(
    strip.text.x = element_text(size = 6, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 6, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 6, hjust = 0, colour = "black"),
    legend.position = "right",
    plot.title = element_blank(),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))+
  coord_flip()


dev.off()


