library(fda)
library(fdaPDE)
library(reshape)
library(tidyverse)
library(tidyfun)
library(ggpubr)
source("Tri_Mesh.r") 
source("HFLM.r")
source("FEM.diff.R")


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
m = 6
trial <- Tri_Mesh(m) 

Nodes = trial$Nodes/(m-1) # make nodes along 0-1

mesh <- create.mesh.2D(nodes = Nodes, segments = trial$Edges, triangles = trial$Triangles, order = 1)

FEMbasis <- create.FEM.basis(mesh)

basis_co <- diag(1, FEMbasis$nbasis)

FEMbasis1 <- FEM(basis_co, FEMbasis)

nfine = 101

# set up lists for concurrent regression (center x and y)

x.fd = center.fd(hip.fd)
xfdlist = list(x.fd)


betalist = list(xybasis)
betalist1 = list(FEMbasis1)

y.fd = center.fd(knee.fd)

##########################################################
#### Run 2 step estimation process########################
##########################################################

# fit concurrent model
tfine = seq(0,1, length = nfine)
Xs = eval.fd(tfine, x.fd)
Ys = eval.fd(tfine, y.fd)

sigX = (1/n)*Ys%*%t(Ys)
sigXY = (1/n)*Ys%*%t(Xs)


alpha.hat = diag(sigXY/sigX)

plot(tfine, alpha.hat, type = 'l')


# extract residuals

conc.hat = apply(Xs, 2, function(p){p*alpha.hat})
conc.hat.fd = smooth.basis(tfine, conc.hat, xybasis)$fd
matplot(conc.hat, type = 'l')
ey = Ys - conc.hat
ey.fd = smooth.basis(tfine, ey, xybasis)$fd

par(mfrow = c(1,1))
plot(ey.fd)

# remove concurrent component

vs = sigX/diag(sigX)

delta = Xs - vs%*%Xs

delta.fd = smooth.basis(tfine, delta, xyPar)$fd
plot(delta.fd)


# fit the HFLM element
betalist1 = list(FEMbasis1)
xfdlist1 = list(delta.fd)

fit = HFLM(y.fd = ey.fd, xfdlist = xfdlist1, betalist = betalist1, nfine = 101, lambda = 1e-10, fit.int = F)



# plot residuals
par(mfrow = c(1,2))
plot(fit$y.fd-fit$yhat.fd, main = "Full")
plot(hist.fit$y.fd-hist.fit$yhat.fd, main = "Historical")

# calculate R2

Cr2 = round(1 - (sum(inprod((y.fd-conc.hat.fd)^2))/sum(inprod(y.fd^2))),2)

yhatfull.fd = (conc.hat.fd+fit$yhat.fd)

Fr2 = round(1 - (sum(inprod((y.fd-yhatfull.fd)^2))/sum(inprod(y.fd^2))),2)


# plot B(s,t) surfaces
nfine = 101
tfine <- seq(0,1, length.out = nfine)
obspts = cbind(rep(tfine,nfine),rep(tfine,1,each=nfine))

model = matrix(eval.FEM(fit$beta.est.list[[1]], obspts), nfine, nfine) %>% melt()
model2 = matrix(eval.FEM(hist.fit$beta.est.list[[1]], obspts), nfine, nfine) %>% melt()

F1 = ggplot(model, aes(x = X1, y = X2, fill = value))+
  geom_tile()+
  theme_light()+
  labs(fill = 'Beta', x = "s", y = "t", 
       title = "Historical Effect",
       subtitle = paste0("Variance Explained = ", 100*(Fr2-Cr2), "%"))+
  scale_fill_gradientn(colors = viridis::viridis(100))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))


H1 = ggplot(model2, aes(x = X1, y = X2, fill = value))+
  geom_tile()+
  theme_light()+
  labs(fill = 'Beta', x = "s", y = "t", 
       title = "Historical Model",
       subtitle = paste0("R-squared = ", Hr2))+
  scale_fill_gradientn(colors = viridis::magma(100), limits = c(-20,20))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))

# plot B(t)

Bt = alpha.hat

C1 = ggplot()+
  geom_line(aes(x = tfine*100, y = Bt), colour = "navy", lwd = 1.5)+
  theme_light()+
  labs(x = "t", y = "Beta", 
       title = "Concurrent Effect",
       subtitle = paste0("Variance Explained = ", 100*Cr2, "%"))



mean.fd = fd(matrix(rep(mean.fd(knee.fd)$coefs, n), y.fd$basis$nbasis, n), y.fd$basis)

tidy.data = data.frame("ID" = colnames(knee.mat))
tidy.data$True = tfd(t(eval.fd(tfine, knee.fd)))
tidy.data$Full = tfd(t(eval.fd(tfine, mean.fd+conc.hat.fd+fit$yhat.fd)))



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

ggsave("HFLM Figure.jpg", plot = last_plot(), dpi = 300, units = 'cm', height = 20, width = 20)

