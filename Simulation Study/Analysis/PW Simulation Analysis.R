library(fda)
library(fdaPDE)
library(tidyverse)
library(reshape)
library(fields)

setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement")

sim.results <- readRDS("./Simulation Study/Results/PW2 Simulation Results 20.rds")
sim.data <- readRDS("./Simulation Study/Data/Simulations Data 20.rds")

plot(sim.results$Results[[1]]$alpha$A)

true.beta <- sim.data$Model$beta
true.alpha <- sim.data$Model$alpha

 # set up parameters for evaluating the surface
nfine <- 101
tfine <- seq(0, 1, length = nfine)
obspts <- cbind(rep(tfine,nfine),rep(tfine,1,each=nfine))

## get betas

bmat.list <- lapply(sim.results$Results, function(x){
  
  bm = matrix(NA, 101, 101)
  bm[upper.tri(bm)] = x$beta$B
  
  bm
  
})


beta.mean <- bmat.list[[1]]

for (i in 2:1000) {
  
  beta.mean <- beta.mean + bmat.list[[i]]
  
}


beta.mean = matrix(beta.mean/1000, nfine, nfine)

fields::image.plot(beta.mean)

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
  x$alpha$A
})

a.co <- do.call(cbind, alpha.list)

bs = create.bspline.basis(c(0,1), 101, 4)

alpha.mean.fd <- mean.fd(fd(a.co, bs))


aSIE.list <- lapply(alpha.list, function(x){
  
  s.res <- (true.alpha-fd(x, bs))^2
  cat(".")
  sie <- inprod(s.res)
  
  return(sie)
})


plot(alpha.mean.fd)

aSIE.vec <- unlist(aSIE.list)
aMSIE <- mean(aSIE.vec)
aSDSIE <- sd(aSIE.vec)


result.list = list("alpha" = alpha.mean.fd, "beta" = beta.mean,
                   "BMSE" = bMSIE, "BSDSE" = bSDSIE,
                   "AMSE" = aMSIE, "ASDSE" = aSDSIE)
saveRDS(result.list,file = "./Simulation Study/Results/PW2 20 Processed.rds")

