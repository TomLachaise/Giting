#### Packages ####
setwd("C:/Users/tomla/Dropbox/RootFun/Data_scripts/Exploratories/Stats/FunRoot") # Change address to set to source file location
#
source("HighstatLibV12.R")
source("ggpub.R")
source("R_rainclouds.R")
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
require(Hmisc)
require(corrplot)
require(corrr)
require(naniar)
require(factoextra)
require(cowplot)
require(rlist)
require(regclass)
require(ggthemes)

#### Data loading ####

# PCA Scores
#"Above_PC1", "Above_PC2",                                         
#"Below_PC1", "Below_PC2", "Below_PC3", "Below_PC4",
#
## Faith index
#"PD", "CWFaith",
#
## Quadratic Entropy
#"TD.QE", "PD.QE", "FD_Above.QE", "FD_Below.QE"

# Load the dataset of SEM results

# Eco proc
SEM.results_env_proc <- read.table("SEM.results_env_proc.txt", sep = "\t", na.strings = "NA", h = T, dec = ",")

# Load the dataset of Div EPs
Div_EPs <- read.table("Div_EPs_for_models.txt", sep = "\t", na.strings = "NA", h = T, dec = ",")

#
Exploratory<-unique(Div_EPs$Exploratory)

# Redefine LUI names

# 2008 - 2019 #
names(Div_EPs)[names(Div_EPs) == 'LUI_2008_2018']    <- 'LUI'
names(Div_EPs)[names(Div_EPs) == 'G_std_2008_2018']  <- 'G_std'
names(Div_EPs)[names(Div_EPs) == 'M_std_2008_2018']  <- 'M_std'
names(Div_EPs)[names(Div_EPs) == 'F_std_2008_2018']  <- 'F_std'
names(Div_EPs)[names(Div_EPs) == 'MF_2008_2018']     <- 'MF_std'


#### I) Eco proc Heatmaps ####

names(SEM.results_env_proc)
head(SEM.results_env_proc)

# Give readable names to every Ecosystem process
levels(SEM.results_env_proc$lhs)
# Keep absolute indirect effect
SEM.results_env_proc$lhs = gsub(pattern = "_indirect_abs", replacement = "", x = SEM.results_env_proc$lhs)
# Remove raw indirect effect (with both directions)
SEM.results_env_proc<- dplyr::filter(SEM.results_env_proc, !grepl('_indirect', lhs))
#
SEM.results_env_proc$lhs <- as.factor(SEM.results_env_proc$lhs)
levels(SEM.results_env_proc$lhs)
levels(SEM.results_env_proc$lhs)<- c("Aboveground traits PC1", "Aboveground traits PC2", "Aboveground traits PC3", "Aboveground traits PC4",
                                     "Aboveground Functional Richness","Belowground Functional Richness", "Soils PC1", "Soils PC2",
                                     "Belowground traits PC1", "Belowground traits PC2", "Belowground traits PC3", "Belowground traits PC4",
                                     "Aboveground productivity",
                                     "CWM Bud bank size", "CWM Plant height", "CWM Leaf dry matter content",
                                     "CWM Leaf surface", "CWM Fine roots diameter", "CWM Root tissue density",
                                     "CWM Specific leaf area", "CWM Specific root length", "CWM Maximum rooting depth",
                                     "CWM Seed mass", "Phylogenetic Richness Faith Index" ,"Litter decomposition", "Land-use intensity",
                                     "Functional MPD Aboveground", "Functional MPD Belowground",
                                     "Functional w_MPD Aboveground", "Functional w_MPD Belowground",
                                     "Phylogenetic MPD ", "Phylogenetic w_MPD", "Leaves fungal infection",
                                     "Soil potential nitrification",
                                     "Species Richness", "Root productivity", "Root decomposition", "Shannon Index")
levels(SEM.results_env_proc$rhs)
levels(SEM.results_env_proc$rhs)<- c("Aboveground traits PC1", "Aboveground traits PC2", "Aboveground traits PC3", "Aboveground traits PC4",
                                     "Absolute environmental indirect effect",
                                     "Aboveground Functional Richness","Belowground Functional Richness", "Soils PC1", "Soils PC2",
                                     "Cumulative environmental indirect effect",
                                     "Belowground traits PC1", "Belowground traits PC2", "Belowground traits PC3", "Belowground traits PC4",
                                     "Aboveground productivity",
                                     "CWM Bud bank size", "CWM Plant height", "CWM Leaf dry matter content",
                                     "CWM Leaf surface", "CWM Fine roots diameter", "CWM Root tissue density",
                                     "CWM Specific leaf area", "CWM Specific root length", "CWM Maximum rooting depth",
                                     "CWM Seed mass", "Phylogenetic Richness Faith Index" ,"Litter decomposition", "Land-use intensity",
                                     "Functional MPD Aboveground", "Functional MPD Belowground",
                                     "Functional w_MPD Aboveground", "Functional w_MPD Belowground",
                                     "Phylogenetic MPD ", "Phylogenetic w_MPD", "Leaves fungal infection",
                                     "Soil potential nitrification",
                                     "Species Richness", "Root productivity", "Root decomposition", "Shannon Index")
#
Exploratories_names <- c("ALB", "HAI", "SCH")
# Change est to Estimate
names(SEM.results_env_proc)[names(SEM.results_env_proc) == 'est'] <- 'Estimate'

# Loop across Exploratories
require(magicfor)
magic_for(silent = TRUE)
for (i in Exploratories_names){

  # Select each exploratory individually
  fitting.SEM.to <- dplyr::filter(SEM.results_env_proc, Exploratory == i)
  # get models estimates with p > 0.05
  fitting.SEM.sig  <- dplyr::filter(fitting.SEM.to, pvalue  > 0.05)
  length(levels(fitting.SEM.sig$Full_model))
  # 176 out of 462 combinations
  fitting.SEM.sig$Full_model <- droplevels(fitting.SEM.sig$Full_model)
  levels(fitting.SEM.sig$Full_model)
  
  #### 1) CWM and Rich ####
  CWM_plantdiv.names<-as.factor(c(
    # Richness
    "Species Richness",
    # PD
    "Phylogenetic MPD",
    "Phylogenetic w_MPD",
    # CWM Above
    "CWM Specific leaf area",
    "CWM Leaf dry matter content",
    "CWM Leaf surface",
    "CWM Plant height",
    "CWM Seed mass",
    # CWM Below
    "CWM Root tissue density",
    "CWM Specific root length",
    "CWM Fine roots diameter",
    "CWM Maximum rooting depth",
    "CWM Bud bank size"
    
  ))
  #
  CWM_fitting.SEM <- dplyr::filter(fitting.SEM.sig, rhs %in% CWM_plantdiv.names)
  
  #Rearrange the order of plant div with SR and PD at top and traits at the bottom
  CWM_fitting.SEM$rhs <- factor(CWM_fitting.SEM$rhs, levels= CWM_plantdiv.names)
  
  # Keep env process as response
  # ecosystem processes
  levels(CWM_fitting.SEM$lhs)
  response.names <-as.factor(c(
    # Aboveground
    "Aboveground productivity",
    "Leaves fungal infection",
    # Belowground
    "Litter decomposition",
    "Root productivity",
    "Root decomposition",
    "Soil potential nitrification"))
  #
  CWM_fitting.SEM     <- dplyr::filter(CWM_fitting.SEM, lhs %in% response.names)
  CWM_fitting.SEM$lhs <- factor(CWM_fitting.SEM$lhs, levels = response.names)
  # keep paths p<0.05
  CWM_fitting.SEM <- dplyr::filter(CWM_fitting.SEM, pvalue.param < 0.05)
  
  # find lowest AIC for each process
  CWM_lowest_AIC <- CWM_fitting.SEM %>% group_by(Process) %>% slice(which.min(aic))
  CWM_lowest_AIC$Best_AIC <- "Lowest AIC"
  # find highest plant div path to each process
  CWM_highest_path <- CWM_fitting.SEM %>% group_by(Process) %>% slice(which.max(abs(Estimate)))
  CWM_highest_path$Best_path <- "Strongest Predictor"
  # merge with estimates dataframe
  CWM_lowest_AIC_tomerge <- dplyr::select(CWM_lowest_AIC, Full_model, Best_AIC)
  CWM_fitting.SEM.Plot <- merge(CWM_fitting.SEM, CWM_lowest_AIC_tomerge[,-1], by ="Full_model", all.x=T)
  CWM_highest_path_tomerge <- dplyr::select(CWM_highest_path, Full_model, Best_path)
  CWM_fitting.SEM.Plot <- merge(CWM_fitting.SEM.Plot, CWM_highest_path_tomerge[,-1], by ="Full_model", all.x=T)
  CWM_fitting.SEM.Plot$Char <- paste(CWM_fitting.SEM.Plot$Best_AIC, CWM_fitting.SEM.Plot$Best_path)
  CWM_fitting.SEM.Plot$Best_AIC[is.na(CWM_fitting.SEM.Plot$Best_AIC)]  <-""
  CWM_fitting.SEM.Plot$Best_path[is.na(CWM_fitting.SEM.Plot$Best_path)]<-""
  
  # plot
  CWM_plot<-ggplot(CWM_fitting.SEM.Plot, aes(lhs, rhs, fill= Estimate)) + 
    geom_tile() +
    geom_text(aes(label = Best_AIC), nudge_y  = 0.1) +
    geom_text(aes(label = Best_path), nudge_y = -0.1) +
    scale_fill_gradient2(midpoint=0, low="#B2182B", high="#2166AC", limits = c(-0.75,0.75))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # loadings on processes for lowest AIC model
  #response.names2<-c(paste(response.names), "Indirect effect")
  CWM_loadings_AIC_process <- dplyr::filter(fitting.SEM.sig, Full_model %in% CWM_lowest_AIC$Full_model)
  CWM_loadings_AIC_process2 <- dplyr::filter(CWM_loadings_AIC_process, lhs %in% response.names)
  CWM_loadings_AIC_process3 <- dplyr::filter(CWM_loadings_AIC_process2, op == "~" | op == ":=")
  
  CWM_Load_proc.plot <- ggplot(CWM_loadings_AIC_process3, aes(x=rhs, y=Estimate)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(~lhs, ncol=3, scales="free_x")+
    ylab("Standardized coefficients")+
    labs(fill = "Coefficient type")+
    theme(
      strip.text.x = element_text(size = 10, colour = "black"),
      axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
      axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
      legend.position = c(0.87, 0.82),
      legend.title = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black", size = 10),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=10),
      #axis.ticks = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white"))
  
  # loadings on processes for highest div path model
  CWM_loadings_DIV_process  <- dplyr::filter(fitting.SEM.sig, Full_model %in% CWM_highest_path$Full_model)
  CWM_loadings_DIV_process2 <- dplyr::filter(CWM_loadings_DIV_process, lhs %in% response.names)
  CWM_loadings_DIV_process3 <- dplyr::filter(CWM_loadings_DIV_process2, op == "~" | op == ":=")
  
  CWM_Load_DIV_proc.plot <- ggplot(CWM_loadings_DIV_process3, aes(x=rhs, y=Estimate)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(~lhs, ncol=3, scales="free_x")+
    ylab("Standardized coefficients")+
    labs(fill = "Coefficient type")+
    theme(
      strip.text.x = element_text(size = 10, colour = "black"),
      axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
      axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
      legend.position = c(0.87, 0.82),
      legend.title = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black", size = 10),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=10),
      #axis.ticks = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white"))
  
  #### 2) Diversities Above and Below ####
  
  # keep CWM and Rich
  Div_plantdiv.names<-as.factor(c(
    # TD Species richness and Shannon
    "Species Richness",
    "Shannon Index",
    # PD Faith and phylo MPD
    "Phylogenetic Richness Faith Index",
    "Phylogenetic MPD",
    # FD: PCA Scores
    "Aboveground traits PC1", "Aboveground traits PC2", "Aboveground traits PC3", "Aboveground traits PC4",                                        
    "Belowground traits PC1", "Belowground traits PC2", "Belowground traits PC3", "Belowground traits PC4",
    # FD: MPD
    "Functional MPD Aboveground", "Functional MPD Belowground",
    "Functional w_MPD Aboveground", "Functional w_MPD Belowground"
  ))
  #
  Div_fitting.SEM <- dplyr::filter(fitting.SEM.sig, rhs %in% Div_plantdiv.names)
  #Rearrange the order of plant div with SR and PD at top and traits at the bottom
  Div_fitting.SEM$rhs <- factor(Div_fitting.SEM$rhs,
                               levels= Div_plantdiv.names)
  
  # keep env process as response
  response.names <-as.factor(c(
    # Aboveground
    "Aboveground productivity",
    "Leaves fungal infection",
    # Belowground
    "Litter decomposition",
    "Root productivity",
    "Root decomposition",
    "Soil potential nitrification"))
  #
  Div_fitting.SEM <- dplyr::filter(Div_fitting.SEM, lhs %in% response.names)
  Div_fitting.SEM$lhs <- factor(Div_fitting.SEM$lhs, levels = response.names)
  # keep paths p<0.05
  Div_fitting.SEM <- dplyr::filter(Div_fitting.SEM, pvalue.param < 0.05)
  
  # find lowest AIC for each process
  Div_lowest_AIC <- Div_fitting.SEM %>% group_by(Process) %>% slice(which.min(aic))
  Div_lowest_AIC$Best_AIC <- "Lowest AIC"
  # find highest plant div path to each process
  Div_highest_path <- Div_fitting.SEM %>% group_by(Process) %>% slice(which.max(abs(Estimate)))
  Div_highest_path$Best_path <- "Strongest Predictor"
  # merge with estimates dataframe
  Div_lowest_AIC_tomerge <- dplyr::select(Div_lowest_AIC, Full_model, Best_AIC)
  Div_fitting.SEM.Plot <- merge(Div_fitting.SEM, Div_lowest_AIC_tomerge[,-1], by ="Full_model", all.x=T)
  Div_highest_path_tomerge <- dplyr::select(Div_highest_path, Full_model, Best_path)
  Div_fitting.SEM.Plot <- merge(Div_fitting.SEM.Plot, Div_highest_path_tomerge[,-1], by ="Full_model", all.x=T)
  Div_fitting.SEM.Plot$Char <- paste(Div_fitting.SEM.Plot$Best_AIC, Div_fitting.SEM.Plot$Best_path)
  Div_fitting.SEM.Plot$Best_AIC[is.na(Div_fitting.SEM.Plot$Best_AIC)]<-""
  Div_fitting.SEM.Plot$Best_path[is.na(Div_fitting.SEM.Plot$Best_path)]<-""
  
  # plot
  Div_plot<-ggplot(Div_fitting.SEM.Plot, aes(lhs, rhs, fill= Estimate)) + 
    geom_tile() +
    geom_text(aes(label = Best_AIC), nudge_y  = 0.1) +
    geom_text(aes(label = Best_path), nudge_y = -0.1) +
    scale_fill_gradient2(midpoint=0, low="#B2182B", high="#2166AC", limits = c(-0.70,0.70))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # loadings on processes for lowest AIC model
  
  Div_loadings_AIC_process  <- dplyr::filter(fitting.SEM.sig, Full_model %in% Div_lowest_AIC$Full_model)
  Div_loadings_AIC_process2 <- dplyr::filter(Div_loadings_AIC_process, lhs %in% response.names)
  Div_loadings_AIC_process3 <- dplyr::filter(Div_loadings_AIC_process2, op == "~" | op == ":=")
  
  Div_Load_proc.plot <- ggplot(Div_loadings_AIC_process3, aes(x=rhs, y=Estimate)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(~lhs, ncol=3, scales="free_x")+
    ylab("Standardized coefficients")+
    labs(fill = "Coefficient type")+
    theme(
      strip.text.x = element_text(size = 10, colour = "black"),
      axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
      axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
      legend.position = c(0.87, 0.82),
      legend.title = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black", size = 10),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=10),
      #axis.ticks = element_blank(),
      panel.background = element_rect(fill = "white",
                                      colour = "white"))
  
  # loadings on processes for highest div path model
  Div_loadings_DIV_process  <- dplyr::filter(fitting.SEM.sig, Full_model %in% Div_highest_path$Full_model)
  Div_loadings_DIV_process2 <- dplyr::filter(Div_loadings_DIV_process, lhs %in% response.names)
  Div_loadings_DIV_process3 <- dplyr::filter(Div_loadings_DIV_process2, op == "~" | op == ":=")
  
  Div_Load_DIV_proc.plot <- ggplot(Div_loadings_DIV_process3, aes(x=rhs, y=Estimate)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(~lhs, ncol=3, scales="free_x")+
    ylab("Standardized coefficients")+
    labs(fill = "Coefficient type")+
    theme(
      strip.text.x = element_text(size = 10, colour = "black"),
      axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
      axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
      legend.position = c(0.87, 0.82),
      legend.title = element_text(color = "black", size = 10),
      legend.text = element_text(color = "black", size = 10),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=10),
      #axis.ticks = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white"))
  
  put( # CWM
      CWM_plot, CWM_Load_proc.plot,
      CWM_fitting.SEM, CWM_loadings_AIC_process3,
      CWM_loadings_DIV_process3,
      # Diversities
      Div_plot,  Div_Load_proc.plot,
      Div_fitting.SEM,
      Div_loadings_AIC_process3,
      Div_loadings_DIV_process3
  )
  
  
}

#
Estimates.plots <- magic_result()

#### II) Plot loadings all Explo combined for AIC and best DIV predict  ####

#### II.1) Div ####

Div_Loadings_AIC_Explo <- bind_rows(Estimates.plots[["Div_loadings_AIC_process3"]])
#

Div_Loadings_AIC_Explo.plot <- ggplot(Div_Loadings_AIC_Explo, aes(x=rhs, y=Estimate, fill = Exploratory)) + 
  geom_bar(position=position_dodge2(preserve="single"), stat="identity")+
  geom_hline(yintercept = 0)+
  facet_wrap(~lhs, ncol=3, scales="free_x")+
  ylab("Standardized coefficients")+
  labs(fill = "Exploratory")+
  theme(
    strip.text.x = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
    legend.position = c(0.95, 0.25),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))

Div_Loadings_AIC_Explo.plot

# All predictors (consistency)

Div_Loadings_All_Explo <- bind_rows(Estimates.plots[["Div_fitting.SEM"]])

#
Div_Loadings_All_Explo.plot <- ggplot(Div_Loadings_All_Explo, aes(x=rhs, y=Estimate, fill = Exploratory)) + 
  geom_bar(position="stack", stat="identity")+
  geom_hline(yintercept=0)+
  facet_wrap(~lhs, ncol=3, scales="free_x")+
  ylab("Standardized coefficients")+
  labs(fill = "Exploratory")+
  theme(
    strip.text.x = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
    legend.position = c(0.87, 0.12),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))

Div_Loadings_All_Explo.plot



#### II.2) CWM ###
# AIC
CWM_Loadings_AIC_Explo <- bind_rows(Estimates.plots[["CWM_loadings_AIC_process3"]])
#

CWM_Loadings_AIC_Explo.plot <- ggplot(CWM_Loadings_AIC_Explo, aes(x=rhs, y=Estimate, fill = Exploratory)) + 
  geom_bar(position=position_dodge2(preserve="single"), stat="identity")+
  geom_hline(yintercept = 0)+
  facet_wrap(~lhs, ncol=3, scales="free_x")+
  ylab("Standardized coefficients")+
  labs(fill = "Exploratory")+
  theme(
    strip.text.x = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
    legend.position = c(0.87, 0.12),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))

CWM_Loadings_AIC_Explo.plot

# All predictors (consistency)

CWM_Loadings_All_Explo <- bind_rows(Estimates.plots[["CWM_fitting.SEM"]])
#
CWM_Loadings_All_Explo.plot <- ggplot(CWM_Loadings_All_Explo, aes(x=rhs, y=Estimate, fill = Exploratory)) + 
  geom_bar(position="stack", stat="identity")+
  geom_hline(yintercept=0)+
  facet_wrap(~lhs, ncol=3, scales="free_x")+
  ylab("Standardized coefficients")+
  labs(fill = "Exploratory")+
  theme(
    strip.text.x = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
    legend.position = c(0.87, 0.12),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))

CWM_Loadings_All_Explo.plot



#### II.3) Sum of estimates Above and Belowground for processes ####

# With PC #
fitting.SEM.sig  <- dplyr::filter(SEM.results_env_proc, pvalue  > 0.05)
length(levels(fitting.SEM.sig$Full_model))
# 174 out of 462 combinations
fitting.SEM.sig$Full_model <- droplevels(fitting.SEM.sig$Full_model)
levels(fitting.SEM.sig$Full_model)
levels(fitting.SEM.sig$rhs)
PC_plantdiv.names<-as.factor(c(
  # TD Species richness and Shannon
  "Species Richness",
  "Shannon Index",
  # PD Faith and phylo MPD
  "Phylogenetic Richness Faith Index",
  "Phylogenetic MPD",
  # FD: PCA Scores
  "Aboveground traits PC1", "Aboveground traits PC2", "Aboveground traits PC3", "Aboveground traits PC4",                                        
  "Belowground traits PC1", "Belowground traits PC2", "Belowground traits PC3", "Belowground traits PC4",
  # FD: MPD
  "Functional MPD Aboveground", "Functional MPD Belowground",
  "Functional w_MPD Aboveground", "Functional w_MPD Belowground"
    ))
#
PC_fitting.SEM <- dplyr::filter(fitting.SEM.sig, rhs %in% PC_plantdiv.names)

#Rearrange the order of plant div with SR and PD at top and traits at the bottom
PC_fitting.SEM$rhs <- factor(PC_fitting.SEM$rhs, levels= PC_plantdiv.names)

# Keep env process as response
# ecosystem processes
levels(PC_fitting.SEM$lhs)
response.names <-as.factor(c(
  # Aboveground
  "Aboveground productivity",
  "Leaves fungal infection",
  # Belowground
  "Litter decomposition",
  "Root productivity",
  "Root decomposition",
  "Soil potential nitrification"))
#
PC_fitting.SEM     <- dplyr::filter(PC_fitting.SEM, lhs %in% response.names)
PC_fitting.SEM$lhs <- factor(PC_fitting.SEM$lhs, levels = response.names)
#
# sum
PC_fitting.SEM$rhs <- gsub(pattern = " PC1", replacement = "", x = PC_fitting.SEM$rhs)
PC_fitting.SEM$rhs <- gsub(pattern = " PC2", replacement = "", x = PC_fitting.SEM$rhs)
PC_fitting.SEM$rhs <- gsub(pattern = " PC3", replacement = "", x = PC_fitting.SEM$rhs)
PC_fitting.SEM$rhs <- gsub(pattern = " PC4", replacement = "", x = PC_fitting.SEM$rhs)
PC_fitting.SEM$rhs <- gsub(pattern = "Functional MPD Aboveground", replacement = "Aboveground traits", x = PC_fitting.SEM$rhs)
PC_fitting.SEM$rhs <- gsub(pattern = "Functional w_MPD Aboveground", replacement = "Aboveground traits", x = PC_fitting.SEM$rhs)
PC_fitting.SEM$rhs <- gsub(pattern = "Functional MPD Belowground", replacement = "Belowground traits", x = PC_fitting.SEM$rhs)
PC_fitting.SEM$rhs <- gsub(pattern = "Functional w_MPD Belowground", replacement = "Belowground traits", x = PC_fitting.SEM$rhs)

# significant only
PC_fitting.SEM <- dplyr::filter(PC_fitting.SEM, pvalue.param < 0.05)

# stacked bars

PC_Loadings_STACKED_Explo.plot <- ggplot(PC_fitting.SEM, aes(x=rhs, y=abs(Estimate), fill = Exploratory)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~lhs, ncol=3, scales="free_x")+
  ylab("Standardized coefficients")+
  labs(fill = "Exploratory")+
  theme(
    strip.text.x = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 45, size = 10, hjust = 0.8 , colour = "black"),
    axis.text.y = element_text(angle = 0, size = 10, hjust = 0, colour = "black"),
    legend.position = c(0.95, 0.25),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 10),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    #axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "white"))

PC_Loadings_STACKED_Explo.plot



#### III) Exploratory plots ####
#### III.1) Raincloud plots SEM variables predictors ####
names(Div_EPs)
Div_proc.subset <- dplyr::select(Div_EPs, Exploratory,
                                 "Mowing and Fertilization" = MF_std,
                                 "Soils PC1" = Axis1_Soils_Explo,
                                 "Soils PC2" = Axis2_Soils_Explo,
                                 "Aboveground traits PC1" = Above_PC1,
                                 "Aboveground traits PC2" = Above_PC2,
                                 "Belowground traits PC1" = Below_PC1,
                                 "Belowground traits PC2" = Below_PC2,
                                 "Belowground traits PC3" = Below_PC3,
                                 "Belowground traits PC4" = Below_PC4,
                                 "Functional Dispersion Aboveground" = Across_Above_FDis_global,
                                 "Functional Dispersion Belowground" = Across_Below_FDis_global,
                                 "Species Richness" = Rich, 
                                 "MPD Phylogeny" = MPDp_abd)

summary(Div_proc.subset)
str(Div_proc.subset)
#
png("p_PCA_axis_Dispersion_Mypairs.png", width = 10000, height = 5000, units = 'px', res = 300)
Mypairs(Div_proc.subset[,-1])
dev.off()
#
Div_proc.subset_melt <- Div_proc.subset %>% gather("Variable", "Value", - Exploratory)

# split into a list for each Variable (necessary for density plot)
df_spl <- split(Div_proc.subset_melt, f = Div_proc.subset_melt$Variable)

# make a raincloud plot for each data frame 

plist <- lapply(
  seq_along(df_spl), 
  function(x) 
  {
    ggplot(df_spl[[x]], aes(x = Exploratory, y = Value, fill =Exploratory)) +
      geom_point(position = position_jitter(width = .05), size = 1, shape = 20)+
      geom_boxplot(outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
      geom_flat_violin(position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = "black")+
      ggtitle(names(df_spl)[x])+
      theme_cowplot()+
      theme(
        axis.title.x=element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")+
      scale_colour_manual(values = c("#F8766D", "#00BA38" , "#619CFF"))+
      scale_fill_manual(values = c("#F8766D", "#00BA38" , "#619CFF"))
  }
)

png("p_Div_predict.Explo.png", width = 3500, height = 2250, units = 'px', res = 300)

plot_grid(plotlist = plist, ncol = 3)

dev.off()

#### III.2) Raincloud plots SEM variables processes ####
names(Div_EPs)
Div_proc.subset <- dplyr::select(Div_EPs, Exploratory,
                                 "Aboveground productivity"  = Biomass,
                                 "Leaves pathogen infection" = pathogen_infection,
                                 # Belowground
                                 "Litter decomposition"      = Litter_decomposition,
                                 "Root productivity"         = Root_biomass,
                                 "Root decomposition"        = Root_decomposition,
                                 "Potential nitrification"   = Potential_nitrification,
                                 "Mycorrhizal hyphae"        = mAMFhyphae,
                                 "Soil aggregation"          = Aggregation)

summary(Div_proc.subset)
str(Div_proc.subset)
# Count NAs per Exploratory #

DT <- data.table(Div_proc.subset)
NAS<-DT[, lapply(.SD, function(x) sum(is.na(x))) , by = list(Exploratory)]
#
Div_proc.subset_melt <- Div_proc.subset %>% gather("Variable", "Value", - Exploratory)

# split into a list for each Variable (necessary for density plot)
df_spl <- split(Div_proc.subset_melt, f = Div_proc.subset_melt$Variable)

# make a raincloud plot for each data frame 

plist <- lapply(
  seq_along(df_spl), 
  function(x) 
  {
    ggplot(df_spl[[x]], aes(x = Exploratory, y = Value, fill =Exploratory)) +
      geom_point(position = position_jitter(width = .05), size = 1, shape = 20)+
      geom_boxplot(outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
      geom_flat_violin(position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = "black")+
      ggtitle(names(df_spl)[x])+
      theme_cowplot()+
      theme(
        axis.title.x=element_blank(),
        plot.title = element_text(size=12),
        legend.position = "none")+
      scale_colour_manual(values = c("#F8766D", "#00BA38" , "#619CFF"))+
      scale_fill_manual(values = c("#F8766D", "#00BA38" , "#619CFF"))
  }
)

png("p_Div_predict.Explo.png", width = 3500, height = 2250, units = 'px', res = 300)

plot_grid(plotlist = plist, ncol = 3)

dev.off()