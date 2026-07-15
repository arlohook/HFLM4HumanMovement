library(fda)
library(fdaPDE)
library(reshape)
library(tidyverse)


source("./Functions/predict.FEM.HFLM.r")
source("./Functions/make_data.R")

setwd("./Simulation Study/Data")

dgobj = readRDS("Data Gen.rds")


seeds = dgobj$seeds$`20`

Data = make_data(n = 20, nsets = 10, seeds = seeds)

saveRDS("Simulations Data 20.rds")
