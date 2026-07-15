## simulation 

library(fda)
library(fdaPDE)
library(tidyverse)
library(reshape2)
library(parallel)
library(foreach)
library(doRNG)
library(doParallel)


setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement")


source("./Functions/Tri_Mesh.r")
source("./Functions/HFLM.smidge.r")
source("./Functions/FEM.diff.r")

# read in data set

Data = readRDS("./Simulation Study/Data/Simulations Data 100.rds")


# create the basis for the coefficients 


m = 6
tri_dim <- Tri_Mesh(m)
nodes = tri_dim$Nodes/(m-1)
triangles = tri_dim$Triangles
edges = tri_dim$Edges
mesh = create.mesh.2D(nodes = nodes, segments = edges, triangles = triangles, order = 1)

b.basis = create.FEM.basis(mesh)

beta.basis = FEM(coeff = diag(1, b.basis$nbasis), FEMbasis = b.basis)
alpha.basis = create.bspline.basis(rangeval = c(0,1), nbasis = 20, norder = 4)


betalist = list(alpha.basis, beta.basis)

nfine = 301

fit.sim = function(b){

  xfdlist =list(Data[[3]][[b]]$x.fd)
  y.fd = Data[[3]][[b]]$yplusE
  
  A = Sys.time()
   
  out = HFLM.smidge(y.fd = y.fd, xfdlist = xfdlist, betalist = betalist, nfine = nfine, lambda = c(0,0,0), fit.int = T)
  out = append(out, Sys.time()-A)
  names(out)[length(out)] = "Run Time"
  
  return(out)
  
}




# set up cluster
n_cores <- 15
cl <- makeCluster(n_cores)
registerDoParallel(cl)


clusterExport(
  cl,
  c("Data", "nfine", "betalist", "fit.sim", "HFLM.smidge", "FEM.diff"),
  envir = environment()
)

# run simulation

RESULTS  = foreach(L = 1:1000,
                   .packages = c("fda", "fdaPDE", "reshape2", "tidyr", "dplyr"),
                   .combine = function(x, y) {
                     
                     append(x, list(y))},
                   
                   .init = list()) %dopar% {fit.sim(L)}

stopCluster(cl)



# save results

details = "OLS fit, order = 1, m = 6"
sim_results <- list("Details" = details, "Results" = RESULTS)


saveRDS(sim_results, "./Simulation Study/Results/FEM Simulation Results 100.rds")



