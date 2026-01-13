library(fdaPDE)
library(fda)
library(tidyverse)
library(reshape)
library(viridis)

pffr20 = readRDS("./Simulation Results/Cleaned/pffr 20 Processed.rds")
pffr50 = readRDS("./Simulation Results/Cleaned/pffr 50 Processed.rds")
pffr100 = readRDS("./Simulation Results/Cleaned/pffr 100 Processed.rds")

FDboost20 = readRDS("./Simulation Results/Cleaned/FDboost 20 Processed.rds")
FDboost50 = readRDS("./Simulation Results/Cleaned/FDboost 50 Processed.rds")
FDboost100 = readRDS("./Simulation Results/Cleaned/FDboost 100 Processed.rds")

OLS20 = readRDS("./Simulation Results/Cleaned/OLS 20 Processed.rds")
OLS50 = readRDS("./Simulation Results/Cleaned/OLS 50 Processed.rds")
OLS100 = readRDS("./Simulation Results/Cleaned/OLS 100 Processed.rds")

Pen20 = readRDS("./Simulation Results/Cleaned/Pen 20 Processed.rds")
Pen50 = readRDS("./Simulation Results/Cleaned/Pen 50 Processed.rds")
Pen100 = readRDS("./Simulation Results/Cleaned/Pen 100 Processed.rds")

Master.Results = list(pffr20, pffr50, pffr100,
                      FDboost20, FDboost50, FDboost100,
                      OLS20, OLS50, OLS100,
                      Pen20, Pen50, Pen100)

nfine = 101
tfine = seq(0,1, length = nfine)
obspts = cbind(rep(tfine, nfine), rep(tfine, each = nfine))

Betas = lapply(Master.Results, function(X){
  
  melted = data.frame(cbind(obspts, melt(X$beta)$value))
  
  return(melted)
  
})

Models = c(rep("pffr",3), rep("FDboost",3), rep("OLS",3), rep("Pen",3))
Sample = rep(c(20, 50, 100),4)

plotting = c()

for(i in 1:12){
  
  mod_frame = Betas[[i]]
  
  mod_frame$Model = Models[i]
  mod_frame$Sample = Sample[i]
  colnames(mod_frame) = c("s", "t", "Beta", "Model", "Sample")
  
  plotting = rbind(plotting, mod_frame)
  
}

plotting$Sample = factor(plotting$Sample, levels = c(20, 50, 100))


Tiles = ggplot(plotting, aes( x = s, y = t, fill = Beta))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis(20), limits = c(-65, 80), breaks = c(-50, 0 ,50))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Value")+
  theme(strip.background = element_rect(fill = "white", colour = "black"), 
        strip.text = element_text(colour = "black", face = 'bold'))+
  facet_grid(Sample ~ Model)


Tiles

ggsave("HFLM Beta Results.jpg", Tiles, dpi = 300, width = 2200, height = 1500, units = 'px')
