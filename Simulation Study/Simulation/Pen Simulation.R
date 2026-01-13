## simulation 

library(fda)
library(fdaPDE)
library(tidyverse)
library(reshape2)
source("Tri_Mesh.r")
source("HFLM.r")
source("FEM.diff.r")



# read in data set

Data = readRDS("Simulations Data 20.rds")


# create the basis for the coefficients 


m = 9
tri_dim <- Tri_Mesh(m)
nodes = tri_dim$Nodes/(m-1)
triangles = tri_dim$Triangles
edges = tri_dim$Edges
mesh = create.mesh.2D(nodes = nodes, segments = edges, triangles = triangles, order = 2)
#plot(mesh)
b.basis = create.FEM.basis(mesh)

beta.basis = FEM(coeff = diag(1, b.basis$nbasis), FEMbasis = b.basis)
alpha.basis = create.bspline.basis(rangeval = c(0,1), nbasis = 50, norder = 4)


betalist = list(alpha.basis, beta.basis)

nfine = 2*b.basis$nbasis+1

Fit.List <- vector('list', 1000)
Fit.List = readRDS("Pen Simulation Results 20 1-167.rds")[[2]]
lambdas = c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e+0)
            #1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7, 1e+8, 1e+9, 1e+10)

START = Sys.time()

for(b in 1:1000) {
  
  
  
  print(paste0("Simulation ", b, " of 1000"))
  
  if(!is.null(Fit.List[[b]])){next}
  # read in data
  xfdlist =list(Data[[3]][[b]]$x.fd)
  y.fd = Data[[3]][[b]]$yplusE
  
  GCV = c()
  
    for(i in 1:length(lambdas)){
  
      print(paste0("  Testing lambda = ", lambdas[i]))

        train.fit <- HFLM(y.fd = y.fd, 
                          xfdlist = xfdlist, 
                          betalist = betalist, 
                          nfine = nfine, 
                          lambda = lambdas[i], fit.int = T)
    
        GCV = append(GCV, train.fit$GCV) 
  
  }
  
 
  l = which.min(GCV)
  
  Fit.List[[b]] = HFLM(y.fd = y.fd, xfdlist = xfdlist, betalist = betalist, nfine = nfine, lambda = 0, fit.int = T)
  
}

plot(Fit.List[[1]]$yhat.fd)


details = "Penalised fit, order = 2, m = 9, lambda selected from a fine grid by minimising GCV"


sim_results <- list("Details" = details, "Results" = Fit.List)

saveRDS(sim_results, "Pen Simulation Results 20.rds")

Sys.time() - START


