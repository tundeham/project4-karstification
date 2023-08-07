# Monte Carlo simulation for recession model
# Double exponential (Maillet) recession model

# rm(list = ls())
# wd <- "D:/Project_2"
# setwd(wd)
# source("./Recession analysis/Recession model/KGE.R")

# check install packages
pkg <- c("hydroGOF","parallel","foreach","doParallel")
for(i in 1:length(pkg)){
  flag <- pkg[[i]] %in% installed.packages()
  if(flag==FALSE){
    install.packages(pkg[[i]])
  }
}
library(hydroGOF)
library(parallel)
library(foreach)
library(doParallel)

#--functions required----
rsq <- function(x, y) summary(lm(y~x))$r.squared
#-----
# set min value to take
setMin <- function(expr){
  pmax(0, eval(expr))
}
# Function to find intercept and slope
intcep <- function(x, y){
  if(is.na(y)){
    return(NA)
  }else{
    abs(summary(lm(y~x))$coefficients[1]) 
  }
} 
slope <- function(x, y){
  abs(summary(lm(y~x))$coefficients[2])
}
#-----
# Function to find index of min
minIndex <- function(df, minValue){
  # df=simLS[[1]]; minValue=m[[1]]
  f_ind <- matrix(data=NA, nrow=1, ncol=ncol(df))
  for(i in 1:ncol(df)){
    f_ind[[i]] <- which(df[,i]==minValue[i])-1
  }
  return(f_ind)
}
#-----
# Function for two linear reservior recession model
resMod <- function(data, frac, a1, a2){
  #data = segs
  # check if Q0 is included in df columns
  if(isFALSE("Q0" %in% colnames(data))){
    Q0 <- data[1,2]
  }
  NR <- nrow(data)
  Qend <- data[NR,2]
  #Qd <- Q0-Qend
  data$t <- data[,6]-1
  
  if(is.null(frac)){
    with(data, ((f*Q0) * exp(-(a1)*t)) + setMin(((1-f)*Q0 * exp(-(a2)*t))))
  }else{
	fm <- frac
    if(Q0*fm < Qend){
      fm=NA
      warnings("Lowest matrix fraction defined, reset value to NA")
    }
    with(data, ((fm*Q0) * exp(-(a1)*t)) + setMin(((1-fm)*Q0 * exp(-(a2)*t))))
  }
  
}
#-----
# Function to run MC simulations for recession model
fRes <- function(data, frac, parset, n, COR=F, OUT=T){
  #data=hydrogr[[41]]; frac=0.1; parset=parset; n=n; COR=T; OUT=T
  Qobs <- data[,2]
  NR <- nrow(data)
  r <- 1
  
  # dataframe to store simulation results 
  Qsim <- as.data.frame(matrix(data = NA,nrow = n+4,ncol = NR))
  Qsim_row <- rep("Qsim_", n); row_i <- seq(1,n,1)
  row_names <- paste0(Qsim_row, row_i)
  row.names(Qsim) <- c("Date","Qobs","f","t",row_names)
  Qsim[1,] <- format(data[,1], format = "%d/%m/%Y")
  Qsim[2,] <- data[,2]
  
  # if initail matrix fraction (f) is supplied with recession df
  if(is.null(frac)){
    Qsim[3,] <- data$f
  # if f value is provided separately
  }else{
    Qsim[3,] <- frac
  }
  #Qsim[4,] <- seq(0,(NR-1),1)
  Qsim[4,] <- data[,6]-1
  
  # run parallel monte carlo sim
  a1 <- parset[,2] #baseflow coefficient
  a2 <- parset[,3] #fastflow coefficient
  Qj <- mapply(resMod, a1=a1, a2=a2, MoreArgs=list(data=data, frac=frac))    
  #Qj <- mapply(parset, function(x) resMod(data=data, frac=frac, a1=as.numeric(x[,2]), a2=as.numeric(x[,3])))  
  
  # get simulation with smallest intercept
  if(COR==T){
    r <- min(apply(Qj,2,intcep,x=Qobs))
  }
  
  # monte carlo simulations output
  if(OUT==T){
    #Qjj <- t(Qj)
    Qsim[5:(n+4),] <- t(Qj)
  }
  if(OUT==T){
    return(Qsim)
  }else{
    return(r)
  }
  rm(Qj)
}
#-----
# Function to run MC simulations for recession model
fRes2 <- function(data, frac, parset, n, COR=F, OUT=T){
  # frac = numeric values of initial matrix fraction between 0 and 1
  # parset = dataframe of priori parameters distribution
  # if OUT == T, simulated Q values are returned
  
  # data=data.frame(resDF_proj31[[2]][7])
  # frac=frac
  # parset=parset[1:500, ]
  # n=500; COR=T; OUT=T
  
  Qobs <- data[,2]
  NR <- nrow(data)
  r <- 1
  
  # dataframe to store simulation results 
  Qsim <- as.data.frame(matrix(data = NA,nrow = n+3,ncol = NR))
  Qsim_row <- rep("Qsim_", n); row_i <- seq(1,n,1)
  row_names <- paste0(Qsim_row, row_i)
  row.names(Qsim) <- c("Date","Qobs","t",row_names)
  Qsim[1,] <- format(data[,1], format = "%d/%m/%Y")
  Qsim[2,] <- data[,2]
  
  #Qsim[3,] <- seq(0,(NR-1),1)
  Qsim[3,] <- data[,6]-1
  
  # create df to save best f for each recession pairs
  f <- data.frame(matrix(data=NA, nrow=n, ncol=2))
  colnames(f) <- c("No","f")
  f[,1] <- seq(1,n)
  
  
  # run parallel monte carlo sim
  # create cluster units for parallel run
  # cl <- makeCluster(no.cores)
  # registerDoParallel(cl)
  
  # initial matrix fraction
  if(COR==T){
    
    f[,2] <- foreach(i=1:nrow(parset), .combine=rbind, .export=c("resMod","setMin","intcep")) %dopar% {
      # i=1
      a1 <- parset[i,2]
      a2 <- parset[i,3]
      Qj <- mapply(resMod, frac=frac, MoreArgs=list(data=data, a1=a1, a2=a2))
      cor <- apply(Qj,2,intcep,x=Qobs)
      r <- which(cor==min(cor, na.rm=T))
      f[i,2] <- (r-1)*0.05 
    }
  }
  
  # recession discharge simulations
  if(OUT==T){
    
    Qsim[4:(n+3),] <- foreach(i=1:nrow(parset), .combine=rbind, .export=c("resMod","setMin","intcep")) %dopar% {
      # i=1
      a1 <- parset[i,2]
      a2 <- parset[i,3]
      Qj <- mapply(resMod, frac=frac, MoreArgs=list(data=data, a1=a1, a2=a2))
      cor <- apply(Qj,2,intcep,x=Qobs)
      r <- which(cor==min(cor, na.rm=T))
      Qsim[i+3,] <- Qj[,r]
    }
  }
  # stop cluster
  # stopCluster(cl)
  
  if(COR==T){
    return(f)
  }
  if(OUT==T){
    return(Qsim)
  }
  
  rm(Qj)
}
#-----
# Function to analysis model performance of Monte Carlo simulation
perfEval <- function(X, PAR, n, skip){
  #X <- sim_segs; Qobs <- as.numeric(X[2,]); PAR <- parset; skip <- 3
  # X = Monte Carlo Simulation object
  # Qobs = observed dicharge ts
  # PAR = parameter combination object data used for Monte Carlo
  # n = number of Monte Carlo run
  # skip = number of headers line to skip in MC object
  
  Qobs <- as.numeric(X[2,])
  X <- X[-(1:skip),]
  X <- apply(X,1,as.numeric)
  
  # output dataframe
  batchPerf <- as.data.frame(matrix(data = NA,nrow = n,ncol = 6))
  colnames(batchPerf) <- c("Parset","a_base","a_cond","KGE","RMSE","NSE")
  batchPerf[,1] <- seq(1,n)
  batchPerf[,2] <- PAR[,2]
  batchPerf[,3] <- PAR[,3]
  
  batchPerf[,4] <- apply(X,2,KGE,obs=Qobs)
  batchPerf[,5] <- apply(X,2,rmse,obs=Qobs)
  batchPerf[,6] <- apply(X,2,NSE,obs=Qobs)
  return(batchPerf)
}


