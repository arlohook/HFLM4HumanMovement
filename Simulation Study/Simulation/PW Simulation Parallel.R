## simulation 

library(fda)
library(fdaPDE)
library(mgcv)
library(refund)
library(tidyverse)
library(reshape2)
library(parallel)
library(foreach)
library(doRNG)
library(doParallel)


setwd("C:/Users/arloh/OneDrive/Documents/PhD/Projects/HFLM for Human Movement")




source("./Functions/HFLM.PW.r")


# read in data set

Data = readRDS("./Simulation Study/Data/Simulations Data 20.rds")


fit.sim = function(b){

  x.fd = Data[[3]][[b]]$x.fd
  y.fd = Data[[3]][[b]]$yplusE
  
  A = Sys.time()
   
  out = HFLM.PW(y.fd = y.fd, x.fd = x.fd, Ak = 20, Bk = c(20,20), tn = 101)
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
  c("Data", "fit.sim", "HFLM.PW2"),
  envir = environment()
)

# run simulation

RESULTS  = foreach(L = 1:1000,
                   .packages = c("fda", "mgcv", "reshape2", "tidyr", "dplyr", "refund"),
                   .combine = function(x, y) {
                     
                     append(x, list(y))},
                   
                   .init = list()) %dopar% {fit.sim(L)}

stopCluster(cl)


# save results

details = "PW fit with k=20 smoothing"
sim_results <- list("Details" = details, "Results" = RESULTS)


saveRDS(sim_results, "./Simulation Study/Results/PW2 Simulation Results 20.rds")


ggplot(RESULTS[[1000]]$beta, aes(x = s, y = t, fill = B))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis::viridis(20), limits = c(min(Data$Model$beta$coeff),max(Data$Model$beta$coeff)))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Value")



hmm = fit.sim(1)
hmm2 = fit.sim(1)
ggplot(hmm$beta, aes(x = s, y = t, fill = B))+
  geom_tile()+
  theme_light()+
  scale_fill_gradientn(colours = viridis::viridis(20), limits = c(min(Data$Model$beta$coeff),max(Data$Model$beta$coeff)))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Value")

t1hat = filter(hmm$beta, t == 1)
t1hatpc = filter(hmm2$beta, t == 1)
t1true = eval.FEM(Data$Model$beta, cbind(seq(0,1, 0.01), rep(1,101)))


plot(x = seq(0, 1, 0.01), y = t1true, type = 'l', main = "B(s,t) at t = 1, n = 500")+
  lines(x = t1hat$s, t1hat$B, lty = 2)+
  lines(x = t1hatpc$s, t1hatpc$B, lty = 3)+
  legend(
    "topleft",                 # position
    legend = c("True", "Splines", "PCR"),
    lty = c(1,2,3),
    bty = "n"                   # remove box (optional)
  )
