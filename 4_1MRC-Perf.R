# Analyse results of Monte Carlo simulation
# Calculate models performances
rm(list = ls())

# set work dir
wd <- "E:/Olarinoye/Project 4 - Spatial recession analysis/"
setwd(wd)

# load packages
source("./Source modules/fileMRC.R")
library(hydroGOF)
library(parallel)
library(foreach)
library(doParallel)

# number of logical processors
no.cores <- detectCores()

# load data
load("SpringHydrographs.RData") # treated spring timeseries
rm(spring.data)

# monte carlo simulations output
sims.dir <- paste0("./Results/SIMS/")
sims <- stringr::str_sort(list.files(sims.dir,full.names=T), numeric=T)

# initial matrix fraction output
f.dir <- paste0("./Results/Matrix fraction/")
fs <- stringr::str_sort(list.files(f.dir,full.names=T), numeric=T)

# define start and end of calculation
a <- 51
b <- 314

# # initiate output dir
# dir.create("./Results/SIMS_Merge/")
# dir.create(paste0("./Results/SIMS_Merge/",b))
# out.dir <- paste0("./Results/SIMS_Merge/",b,"/")
# 
# # create variables for child clusters
# cl <- makeCluster(8)
# clusterExport(cl, varlist = c('parLapply',"sims.dir","sims"))
# 
# #--- 1 ---
# # merge recession segments to time series
# # for(i in 1:length(sims)){
# for(i in a:b){
#   cat("Merging recession segments of spring no.",i,"\n")
#   sim <- readRDS(sims[[i]])
# 
#   # round float integers to reduce file size
#   sim.x <- parLapply(cl=cl, sim, function(x) {
#     x <- sapply(x,as.numeric)
#     x <- apply(x,2,round,digits=4)
#     })
#   rm(sim)
# 
#   # merge and save result to disk
#   W <- paste0(out.dir,i,"simTS.rds")
#   con <- file(W,"wb")
#   saveRDS(do.call("cbind",sim.x), file=con)
#   close(con)
#   rm(sim.x)
# }
# stopCluster(cl)


# #--- 2 ---
# # merge initial matrix fractions
# # initiate output dir
# dir.create("./Results/Matrix fraction_Merge/")
# dir.create(paste0("./Results/Matrix fraction_Merge/",b))
# out.dir <- paste0("./Results/Matrix fraction_Merge/",b,"/")
# # for(i in 1:length(fs)){
# for(i in a:b){
#   cat("Merging initial matrix fraction df of spring no.",i,"\n")
#   f <- readRDS(fs[[i]])
# 
#   # merge and save result to disk
#   X <- do.call("cbind",f)
#   X <- X[,which(colnames(X)=="f")]
#   W <- paste0(out.dir,i,"matFracTS.rds")
#   con <- file(W,"wb")
#   saveRDS(X, file=con)
#   close(con)
#   rm(f,X)
# }


#--- 3 ---
# Evaluate model performance
# initiate output dir
dir.create("./Results/Eval")
dir.create(paste0("./Results/Eval/",b))
out.dir <- paste0("./Results/Eval/",b,"/")

# load files
# sims <-  stringr::str_sort(list.files(paste0("./Results/SIMS_Merge/",b),full.names=T), numeric=T)
sims <-  stringr::str_sort(list.files(paste0("./Results/SIMS_Merge/1-314"),full.names=T), numeric=T)
parset <- read.csv(paste0("./Parameters10k.csv"), sep=",", header = T)
n <- nrow(parset)

for(i in a:b){
  cat("Calculating KGE for spring no.",i,"...\n")
  sim <- readRDS(sims[[i]])

  # calculate KGE and save result
  W <- paste0(out.dir,i,"KGE.rds")
  con <- file(W,"wb")

  # saveRDS(perfEval(sim, PAR=parset, n=n, skip=3), file=con)

  X <- perfEval(sim, PAR=parset, n=n, skip=3)
  X$id <- metainfo$id[i]
  X$lon <- round(as.numeric(metainfo$lon[i]), 4)
  X$lat <- round(as.numeric(metainfo$lat[i]), 4)
  saveRDS(X, file=con)
  close(con)

  rm(sim,X)
}


#--- 4 ---
# Rank model peformance by KGE; select best model; estimate SD
# load file from Eval and matrix fraction dir
eval <- stringr::str_sort(list.files("./Results/Eval/1-314", full.names=T), numeric=T)
fs <- stringr::str_sort(list.files("./Results/Matrix fraction_Merge/1-314", full.names=T), numeric=T)

for(i in 1:length(eval)){
  #i=1
  cat("Evaluating result of spring no",i, ":)\n")
  perf <- readRDS(eval[[i]])
  f <- readRDS(fs[[i]])
  
  # select simulation with KGE > 0.5
  perf.ord <- perf[order(perf['KGE'],decreasing=T),]
  perf.ord <- subset(perf.ord, KGE >= 0.5)
  
  # select corresponding f
  f.ord <- f[order(perf['KGE'],decreasing=T),] 
  if(i==302 | i==304){
    f.ord <- f.ord[1:nrow(perf.ord)]
  }else{
    f.ord <- f.ord[1:nrow(perf.ord),]
  }
  
  # compute karstication index
  ki <- perf.ord$a_base/perf.ord$a_cond
  
  ki_v2 <- perf.ord$a_cond/perf.ord$a_base
  
  # compute meand and sd of parameters
  metainfo$mean_alpB[i] <- mean(perf.ord$a_base)
  metainfo$sd_alpB[i] <- sd(perf.ord$a_base)
  metainfo$mean_alpC[i] <- mean(perf.ord$a_cond)
  metainfo$sd_alpC[i] <- sd(perf.ord$a_cond)
  metainfo$mean_f[i] <- mean(as.matrix(f.ord))
  metainfo$sd_f[i] <- sd(as.matrix(f.ord))
  metainfo$mean_ki[i] <- mean(ki)
  metainfo$sd_ki[i] <- sd(ki)
  
  metainfo$mean_ki_v2[i] <- mean(ki_v2)
  metainfo$sd_ki_v2[i] <- sd(ki_v2)
  
  rm(perf,f,per.ord,f.ord,ki)
  
}
save(file="metainfo+param.RData",metainfo)

#--- end ---

# verification
# simTS <- list.files(out.dir, full.names=T)
# simTS_1 <- readRDS(simTS[[1]])
# sim_1 <- readRDS(sims[[1]])
# 
# # matrix fraction
# fTS <- list.files(out.dir, full.names=T)
# fTS_1 <- readRDS(fTS[[1]])
# rm(sim_1,simTS_1)

# Eval
# kge <- list.files(out.dir, full.names=T)
# kge_1 <- readRDS(kge[[3]])
# rm(kge,kge_1)
