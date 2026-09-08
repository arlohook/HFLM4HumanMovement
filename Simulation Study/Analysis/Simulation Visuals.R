library(fdaPDE)
library(fda)
library(tidyverse)
library(reshape)
library(viridis)

setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement/Github/Simulation Study")

pffr20 = readRDS("./Results/pffr 20 Processed.rds")
pffr50 = readRDS("./Results/pffr 50 Processed.rds")
pffr100 = readRDS("./Results/pffr 100 Processed.rds")

PW20 = readRDS("./Results/PW 20 Processed.rds")
PW50 = readRDS("./Results/PW 50 Processed.rds")
PW100 = readRDS("./Results/PW 100 Processed.rds")

FEM20 = readRDS("./Results/FEM 20 Processed.rds")
FEM50 = readRDS("./Results/FEM 50 Processed.rds")
FEM100 = readRDS("./Results/FEM 100 Processed.rds")



Master.Results = list(pffr20, pffr50, pffr100,
                      FEM20, FEM50, FEM100,
                      PW20, PW50, PW100)

nfine = 101
tfine = seq(0,1, length = nfine)
obspts = cbind(rep(tfine, nfine), rep(tfine, each = nfine))

Betas = lapply(Master.Results, function(X){
  
  melted = data.frame(cbind(obspts, melt(X$beta)$value))
  
  return(melted)
  
})

Models = c(rep("pffr",3), rep("Finite Elements",3),rep("Point Wise",3))
Sample = rep(c(20, 50, 100),3)

plotting = c()

for(i in 1:9){
  
  mod_frame = Betas[[i]]
  
  mod_frame$Model = Models[i]
  mod_frame$Sample = Sample[i]
  colnames(mod_frame) = c("s", "t", "Beta", "Model", "Sample")
  
  plotting = rbind(plotting, mod_frame)
  
}

plotting$Sample = factor(plotting$Sample, levels = c(20, 50, 100))
plotting$Beta = ifelse(plotting$s>plotting$t, NA, plotting$Beta)

Tiles = ggplot(plotting, aes( x = s, y = t, fill = Beta))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis(20),limits = c(min(plotting$Beta), max(plotting$Beta)))+
  scale_x_continuous(expand = c(0,0), labels = c("0", "0.25", "0.5", "0.75", "1"), breaks = c(0,0.25, 0.5, 0.75, 1))+
  scale_y_continuous(expand = c(0,0), , labels = c("0", "0.25", "0.5", "0.75", "1"), breaks = c(0,0.25, 0.5, 0.75, 1))+
  labs(fill = "Value")+
  theme(strip.background = element_rect(fill = "white", colour = "black"), 
        strip.text = element_text(colour = "black", face = 'bold'))+
  facet_grid(Sample ~ Model)


Tiles

True.Beta  = readRDS("./Data/Simulations Data 20.rds")$Model$beta

B = data.frame(s = obspts[,1], t = obspts[,2], Value = eval.FEM(True.Beta, obspts))

TrB = ggplot(B, aes( x = s, y = t, fill = Value))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis(20), limits = c(min(plotting$Beta), max(plotting$Beta)))+
  scale_x_continuous(expand = c(0,0), labels = c("0", "0.25", "0.5", "0.75", "1"), breaks = c(0,0.25, 0.5, 0.75, 1))+
  scale_y_continuous(expand = c(0,0), labels = c("0", "0.25", "0.5", "0.75", "1"), breaks = c(0,0.25, 0.5, 0.75, 1))+
  labs(fill = "Value", title = "True")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

TrB



library(patchwork)


layout = c("AAABBBBBB
            AAABBBBBB
            CCCBBBBBB
            CCCBBBBBB")

Fig2 = TrB + Tiles + guide_area() + plot_layout(design = layout, guides = 'collect')
Fig2

ggsave("HFLMHM - Figure 5 PW2.jpg", Fig2, dpi = 300, width = 2200, height = 1500, units = 'px')



plotting2 = plotting %>% group_by(Model, Sample) %>% mutate(Dif = Beta - B$Value)

Tiles2 = ggplot(plotting2, aes( x = s, y = t, fill = Dif))+
  geom_tile()+
  theme_light()+
  scale_fill_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0)+
  scale_x_continuous(expand = c(0,0), labels = c("0", "0.25", "0.5", "0.75", "1"), breaks = c(0,0.25, 0.5, 0.75, 1))+
  scale_y_continuous(expand = c(0,0), labels = c("0", "0.25", "0.5", "0.75", "1"), breaks = c(0,0.25, 0.5, 0.75, 1))+
  labs(fill = "Value")+
  theme(strip.background = element_rect(fill = "white", colour = "black"), 
        strip.text = element_text(colour = "black", face = 'bold'))+
  facet_grid(Sample ~ Model)

Tiles2

ggsave("HFLMHM - Figure 6 PW2.jpg", Tiles2, dpi = 300, width = 2200, height = 1500, units = 'px')


TB = data.frame("Model" = Models, "Sample" = Sample,
                "Alpha_mu" = round(sapply(Master.Results, function(x){x$AMSE}),1),
                "Alpha_sd" = round(sapply(Master.Results, function(x){x$ASDSE}),1),
                "Beta_mu" = round(sapply(Master.Results, function(x){x$BMSE}),1),
                "Beta_sd" = round(sapply(Master.Results, function(x){x$BSDSE}),1))


RTpffr20 = readRDS("./Results/pffr Simulation Results 20.rds")$Results %>% sapply(function(X){X$`Run Time`})
RTpffr50 = readRDS("./Results/pffr Simulation Results 50.rds")$Results %>% sapply(function(X){X$`Run Time`})
RTpffr100 = readRDS("./Results/pffr Simulation Results 100.rds")$Results %>% sapply(function(X){X$`Run Time`})

# correct if over a minute 
RTpffr100 = ifelse(RTpffr100 < 5, RTpffr100*60, RTpffr100)

RTFEM20 = readRDS("./Results/FEM Simulation Results 20.rds")$Results %>% sapply(function(X){X$`Run Time`})
RTFEM50 = readRDS("./Results/FEM Simulation Results 50.rds")$Results %>% sapply(function(X){X$`Run Time`})
RTFEM100 = readRDS("./Results/FEM Simulation Results 100.rds")$Results %>% sapply(function(X){X$`Run Time`})

RTPW20 = readRDS("./Results/PW2 Simulation Results 20.rds")$Results %>% sapply(function(X){X$`Run Time`})
RTPW50 = readRDS("./Results/PW2 Simulation Results 50.rds")$Results %>% sapply(function(X){X$`Run Time`})
RTPW100 = readRDS("./Results/PW2 Simulation Results 100.rds")$Results %>% sapply(function(X){X$`Run Time`})

RTS = list(RTpffr20, RTpffr50, RTpffr100,
           RTFEM20, RTFEM50, RTFEM100,
           RTPW20, RTPW50, RTPW100)

TB$`Run_Time_mu` = round(sapply(RTS, mean),1)

t = seq(0,1,0.01)
alphas = lapply(Master.Results, function(X){
  
  A = eval.fd(t, X$alpha)
  
  })
alphas = do.call(cbind, alphas)

colnames(alphas) = c("pffr_20", "pffr_50", "pffr_100",
                     "Finite Elements_20", "Finite Elements_50", "Finite Elements_100",
                     "Pointwise_20", "Pointwise_50", "Pointwise_100")

alphas = data.frame(melt(alphas))
alphas$Var1 = rep(t, 9)
alphas$Mu = rep(c(eval.fd(t, readRDS("./Data/Simulations Data 20.rds")$Model$alpha)), 9)

colnames(alphas) = c("t", "Mod_Sam", "Hat", "True")

alphas = alphas %>% separate(Mod_Sam, into = c("Model", "Sample"), sep = "_") %>% 
          mutate(Sample = factor(Sample, levels = c("20", "50", "100"), ordered = T))

A.fig = ggplot(alphas, aes(x = t, colour  = Model))+
  geom_line(aes(y = Hat))+
  geom_line(aes(y = True), colour = "black")+
  theme_light()+
  labs(y = "Value")+
  geom_hline(yintercept = 0, lty = 2, colour = "black")+
  facet_grid(Sample ~ Model)+
  theme(legend.position = 'none')

A.fig
ggsave("HFLMHM - Figure 4 PW2.jpg", A.fig, dpi = 300, width = 2200, height = 1500, units = 'px')



TrB2 = ggplot(B, aes( x = s, y = t, fill = Value))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis(20), limits = c(min(plotting$Beta), max(plotting$Beta)))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Value")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

TrA = ggplot(filter(alphas, Model == "pffr" & Sample == 20), aes(x = t, y = True))+
        geom_line()+
        labs(y = "Value")+
        geom_hline(yintercept = 0, lty = 2, colour = "black")+
        theme_light()
TrA
TrB2

Fig = ggarrange(TrA, TrB2, widths = c(0.45, 0.55))

ggsave("HFLMHM - Figure 3.jpg", plot = Fig, dpi = 300, units = 'cm', height = 9, width = 20)





