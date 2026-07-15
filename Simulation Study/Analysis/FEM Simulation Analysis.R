library(fda)
library(fdaPDE)
library(tidyverse)
library(reshape)
library(fields)

setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement")

sim.results <- readRDS("./Simulation Study/Results/FEM Simulation Results 100.rds")
sim.data <- readRDS("./Simulation Study/Data/Simulations Data 100.rds")

true.beta <- sim.data$Model$beta
true.alpha <- sim.data$Model$alpha

 # set up parameters for evaluating the surface
nfine <- 101
tfine <- seq(0, 1, length = nfine)
obspts <- cbind(rep(tfine,nfine),rep(tfine,1,each=nfine))

## get betas

bmat.list <- lapply(sim.results$Results, function(x){
  
  bm <- eval.FEM(FEM = x$beta.est.list[[1]], locations = obspts)
  
})



beta.mean <- bmat.list[[1]]

for (i in 2:1000) {
  
  beta.mean <- beta.mean + bmat.list[[i]]
  
}


beta.mean = matrix(beta.mean/1000, nfine, nfine)



true.b.mat <- matrix(eval.FEM(FEM = true.beta, obspts), nfine, nfine)



# get mean and sd SSIE

integrate_surface_error <- function(B, Bhat, s, t) {
  stopifnot(all(dim(B) == dim(Bhat)))
  n1 <- length(s)
  n2 <- length(t)
  stopifnot(n1 == nrow(B), n2 == ncol(B))
  
  # squared error matrix
  SE <- (B - Bhat)^2
  
  # trapezoidal weights for each axis
  ws <- c(0.5, rep(1, n1 - 2), 0.5) * (s[2] - s[1])
  wt <- c(0.5, rep(1, n2 - 2), 0.5) * (t[2] - t[1])
  
  # tensor product weights
  W <- outer(ws, wt, "*")
  
  # approximate integral
  sum(SE * W)
}

B = true.b.mat
B[is.na(B)] = 0


bSIE.list <- lapply(bmat.list, function(x){
  
  Bhat = matrix(x, nfine, nfine)
  Bhat[is.na(Bhat)] = 0
  
  integrate_surface_error(B = B, Bhat = Bhat, s = tfine, t = tfine)

  
})

bSIE.vec <- unlist(bSIE.list)


bMSIE <- mean(bSIE.vec)
bSDSIE = sd(bSIE.vec)


# get alpha

alpha.list <- lapply(sim.results$Results, function(x){
  x$intercept.fd
})

a.co <- c()

for (i in 1:1000) {
  
  a.co = cbind(a.co, alpha.list[[i]]$coefs[,1])
  
}

alpha.mean.fd <- mean.fd(fd(a.co, alpha.list[[1]]$basis))


aSIE.list <- lapply(alpha.list, function(x){
  
  s.res <- (true.alpha-x)^2
  cat(".")
  sie <- inprod(s.res)
  
  return(sie)
})




aSIE.vec <- unlist(aSIE.list)
aMSIE <- mean(aSIE.vec)
aSDSIE <- sd(aSIE.vec)


result.list = list("alpha" = alpha.mean.fd, "beta" = beta.mean,
                   "BMSE" = bMSIE, "BSDSE" = bSDSIE,
                   "AMSE" = aMSIE, "ASDSE" = aSDSIE)
saveRDS(result.list,file = "./Simulation Study/Results/FEM 100 Processed.rds")

