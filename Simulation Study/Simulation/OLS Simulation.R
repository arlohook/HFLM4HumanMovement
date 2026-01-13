## simulation 

library(fda)
library(fdaPDE)
library(tidyverse)
library(reshape2)
source("Tri_Mesh.r")
source("HFLM.r")


# read in data set

Data = readRDS("Simulations Data 20.rds")


# create the basis for the coefficients 


m = 9
tri_dim <- Tri_Mesh(m)
nodes = tri_dim$Nodes/(m-1)
triangles = tri_dim$Triangles
edges = tri_dim$Edges
mesh = create.mesh.2D(nodes = nodes, segments = edges, triangles = triangles, order = 1)
plot(mesh)
b.basis = create.FEM.basis(mesh)

beta.basis = FEM(coeff = diag(1, b.basis$nbasis), FEMbasis = b.basis)
alpha.basis = create.bspline.basis(rangeval = c(0,1), nbasis = 50, norder = 4)


betalist = list(alpha.basis, beta.basis)

nfine = 201

Fit.List <- vector('list', 1000)

for (b in 1:1000) {
  


print(paste0("Simulation ", b, " of 1000"))


# read in data
  x.fd = Data[[3]][[b]]$x.fd
  y.fd = Data[[3]][[b]]$yplusE

## fit the model on sim data

xfdlist = list(x.fd)

Fit.List[[b]] = HFLM(y.fd = y.fd, xfdlist = xfdlist, betalist = betalist, nfine = nfine, lambda = 0)



}


plot(Fit.List[[1]]$yhat.fd)
plot(Fit.List[[1]]$xfdlist[[1]])
# save sim data to save time 

details = "OLS fit, order = 1, m = 9"


sim_results <- list("Details" = details, "Results" = Fit.List)
saveRDS(sim_results, "OLS Simulation Results 20.rds")

