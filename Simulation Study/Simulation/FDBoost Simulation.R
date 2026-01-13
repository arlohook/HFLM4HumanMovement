## simulation 
library(FDboost)
library(fda)
library(tidyverse)
library(reshape2)



# read in data set
n = 20
Data = readRDS("Simulations Data 20.rds")

nfine = 101
tfine = seq(0, 1, length = nfine)

Fit.List <- vector('list', 1000)

get.info <- function(model){
  
  alpha <- rowMeans(matrix(modelF$offset, 101, 100, T))
  D = sqrt(length(model$coef()[[2]]))
  beta = matrix(model$coef()[[2]], D, D)
  residuals <- matrix(model$resid(), n, nfine)  
  
  return(list("alpha" = alpha, "beta" = beta, "residuals" = residuals))
  
}

A = Sys.time()

for (b in 1:1000) {
  

print(paste0("Simulation ", b, " of 1000"))


# read in data
  
df = list(Y = t(eval.fd(tfine, Data[[3]][[b]]$yplusE)), 
          X = t(eval.fd(tfine, Data[[3]][[b]]$x.fd)),
          t = tfine,
          s = tfine)
## fit the model on sim data

modelF = FDboost(Y ~ 1 + bhist(x = X, s = s, time = t, limits = "s<=t", knots = 20),
              timeformula = ~ bbs(t, knots = 20), data = df)



folds     <- cv(rep(1, length(unique(modelF$id))), B = 10)
boostIt   <- applyFolds(modelF, folds = folds, grid = seq(100, 2500, by = 400))

modelF <- modelF[mstop(boostIt)]



Fit.List[[b]] = get.info(model = modelF)



}


# save sim data to save time 

details = "FDboost fit, knots = 20"


sim_results <- list("Details" = details, "Results" = Fit.List)
saveRDS(sim_results, "FDboost Simulation 20 Results.rds")

A - Sys.time()


