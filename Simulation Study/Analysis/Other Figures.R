library(fdaPDE)
library(fda)
library(tidyverse)
library(reshape)
library(viridis)
library(ggforce)


source("./Functions/Tri_Mesh.R")


m = 5

trial <- Tri_Mesh(m) 

Nodes = trial$Nodes/(m-1) # make nodes along 0-1

mesh <- create.mesh.2D(Nodes, segments = trial$Edges, triangles = trial$Triangles, order = 1)

FEMbasis <- create.FEM.basis(mesh)

FEMco <- c(0, 0, 0, 0, 0, 0, 0, 0 , 2, 0, 0, 0, 0 , 0, 0)

FEMobj <- FEM(coeff = FEMco, FEMbasis = FEMbasis)

t = seq(0,1,0.001)
obspts = cbind(rep(t, 1001), rep(t, each = 1001))

plt = data.frame("s" = rep(t, 1001), "t" = rep(t, each = 1001), "Value" = eval.FEM(FEMobj, obspts))

ggplot(plt, aes(x = s, y = t, fill = Value))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = turbo(20))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_circle(
    aes(x0 = 0.5, y0 = 0.75, r = 0.125),
    lty = 2,
    fill = NA,        # transparent fill
    colour = "black"
  )

ggsave("HFLMHM - Figure 1.jpg", dpi = 300, height = 8, width = 10, units = 'cm')
