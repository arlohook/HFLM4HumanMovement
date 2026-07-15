library(fda)
library(fdaPDE)
library(reshape)
library(tidyverse)
library(tidyfun)
library(ggpubr)
library(refund)

setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement")

source("./Functions/Tri_Mesh.r") 
source("./Functions/HFLM.smidge.r")
source("./Functions/FEM.diff.R")

# create fdobj
hip.mat = gait[,,1]
knee.mat = gait[,,2]

n = ncol(hip.mat)

time = seq(0,1, length = 20)

xybasis = create.bspline.basis(c(0,1), 10, 4)
xyPar = fdPar(xybasis, 2, 1e-9)
xyPar1 = fdPar(xybasis, 2, 1e-6)
hip.fd = smooth.basis(time, hip.mat, xyPar)$fd
knee.fd = smooth.basis(time, knee.mat, xyPar1)$fd

#plot variables
par(mfrow = c(1,2))
plot(hip.fd)
plot(knee.fd)



# create FEM basis
m = which(cumsum(pca.fd(hip.fd, nharm = 10)$varprop)>0.99)[1]
trial <- Tri_Mesh(m) 

Nodes = trial$Nodes/(m-1) # make nodes along 0-1

mesh <- create.mesh.2D(nodes = Nodes, segments = trial$Edges, triangles = trial$Triangles, order = 1)

FEMbasis <- create.FEM.basis(mesh)

basis_co <- diag(1, FEMbasis$nbasis)

FEMbasis1 <- FEM(basis_co, FEMbasis)

nfine = 101

# set up lists for concurrent regression (center x and y)


xfdlist = list(hip.fd)


betalist = list(xybasis, FEMbasis1)

fit = HFLM.smidge(y.fd = knee.fd, xfdlist = xfdlist, betalist = betalist, nfine = 101, lambda = c(0,0,0), fit.int = T)


 
# calculate R2
SSE = sum(inprod((fit$yhat.fd-knee.fd)^2))
TSE = sum(inprod(center.fd(knee.fd)^2))
R2 = round(1 - SSE/TSE, 3)

# plot B(s,t) surfaces
nfine = 101
tfine <- seq(0,1, length.out = nfine)
obspts = cbind(rep(tfine,nfine),rep(tfine,1,each=nfine))

model = matrix(eval.FEM(fit$beta.est.list[[1]], obspts), nfine, nfine) %>% melt()

F1 = ggplot(model, aes(x = X1/100, y = X2/100, fill = value))+
  geom_tile()+
  theme_light()+
  labs(fill = 'Beta', x = "s", y = "t", 
       title = "Historical Effect",
       subtitle = paste0("Variance Explained = ", 100*R2, "%"))+
  scale_fill_gradientn(colors = viridis::viridis(100))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))


F1

# plot a(t)

at = eval.fd(tfine, fit$intercept.fd)

C1 = ggplot()+
  geom_line(aes(x = tfine, y = at), colour = "navy", lwd = 1.5)+
  theme_light()+
  labs(x = "t", y = "Beta", 
       title = "Intercept",
       subtitle = " ")

C1

mean.fd = fd(matrix(rep(mean.fd(knee.fd)$coefs, n), knee.fd$basis$nbasis, n), knee.fd$basis)

tidy.data = data.frame("ID" = colnames(knee.mat))
tidy.data$True = tfd(t(eval.fd(tfine, knee.fd)), arg = tfine)
tidy.data$Full = tfd(t(eval.fd(tfine, fit$yhat.fd)), arg = tfine)



True = ggplot(tidy.data, aes(y = True, colour = ID))+
  geom_spaghetti(lwd = 1)+
  scale_colour_manual(values = viridis::turbo(39))+
  theme_light()+
  theme(legend.position = 'none')+
  labs(x = 't', y = "Knee Flexion", title = "True Knee")



Full = ggplot(tidy.data, aes(y = Full, colour = ID))+
  geom_spaghetti(lwd = 1)+
  scale_colour_manual(values = viridis::turbo(39))+
  theme_light()+
  theme(legend.position = 'none')+
  labs(x = 't', y = "Knee Flexion", title = "Estimated Knee")





ggarrange(C1, F1, True, Full, ncol = 2, nrow = 2)

ggsave("HFLMHM - Figure 9.1.jpg", plot = last_plot(), dpi = 300, units = 'cm', height = 20, width = 20)


