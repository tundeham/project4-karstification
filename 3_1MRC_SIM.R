# Project 4 - Master recession analysis of extracted 
# recession events
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
load("RecessionEvents.RData") # extracted recessions from all springs

# define recession parameters sample 
n <- 10000
parset<-as.data.frame(matrix(data = NA,nrow = n,ncol = 3))
colnames(parset)<-c("No","alp_b","alp_c")
parset[,1]<-seq(1,n)
set.seed(100)
parset[,2]<-runif(n, 0.00001,	0.05)
parset[,3]<-runif(n, 0.01,	1)
frac <- seq(0,1,0.05)  # matrix fraction
write.csv(parset, file="Parameters10k.csv", sep=",", row.names = FALSE)


# STEP 1: Run Monte Carlo simulation
# initialize output folder
dir.create("Results")
dir.create("Results/SIMS/")
# dir.create("Results/SIMS_")
sim.dir <- paste0("./Results/SIMS/")
# dir.create("Results/Matrix fraction/")
# f.dir <- paste0("./Results/Matrix fraction/")

# function to run parallel recession model 'resMod'
parfRes2 <- function(data,frac,parset,n,COR=F,OUT=T,W=NULL){
  # W = output file dir
  if(is.null(W)){
    mclapply(data, function(x){
      fRes2(data=x, frac=frac, parset=parset, n=n, COR=COR, OUT=OUT)
    })
  }else{
    con <- file(W, "wb")
    saveRDS(mclapply(data, function(x){
      fRes2(data=x, frac=frac, parset=parset, n=n, COR=COR, OUT=OUT)
    }), 
    file=con)
    close(con)
  }
  
}

# another variation of parfRes2 function above using "parLapply"
# this should be faster and for creating multiple cores
parfRes3 <- function(data,frac,parset,n,COR=F,OUT=T,W=NULL){
  # W = output file dir
  if(is.null(W)){
    parLapply(cl=cl, X=data, fun=function(x){
      fRes2(data=x, frac=frac, parset=parset, n=n, COR=COR, OUT=OUT)
    })
  }else{
    con <- file(W, "wb")
    saveRDS(parLapply(cl=cl, X=data, fun=function(x){
      fRes2(data=x, frac=frac, parset=parset, n=n, COR=COR, OUT=OUT)
    }), 
    file=con)
    close(con)
  }
}

# create variables for child clusters
# cl <- makeCluster(no.cores)
cl <- makeCluster(no.cores)
clusterExport(cl, varlist = c('parfRes3','frac','fRes2','n','parset','resMod',
                              'setMin','intcep','parLapply',"%dopar%","foreach"))

# initial matrix fraction of recession segments
# fTime1 <- Sys.time()
# for(i in 1:80){
  # cat("Calculating initial matrix fraction for spring no.",i,"\n")
  # res <- res.list[[i]]
  # W <- paste0(f.dir,i,"matFrac.rds")
  # parfRes3(data=res, W=W, frac=frac, parset=parset, n=n, COR=T, OUT=F)
# }
# fTime2 <- Sys.time()

# simulations of recession events
simTime1 <- Sys.time()
for(i in 103:150){
  cat("running MC simulation for spring no.",i,"\n")
  res <- res.list[[i]]
  W <- paste0(sim.dir,i,"sim.rds")
  parfRes3(data=res, W=W, frac=frac, parset=parset, n=n, COR=F, OUT=T)
}
simTime2 <- Sys.time()



# check codes
# ftest <- readRDS(paste0(f.dir,"5matFrac.rds"))
# simtest <- readRDS(paste0(sim.dir,"5sim.rds"))
# length(res.list[[5]])
# rm(ftest)
# df <- res.list[[151]]

# tryFunction <- function(x){
#   out <- list()
#   out[[1]] <- lapply(x, function(x) x+x)
#   out[[2]] <- lapply(x, function(x) x^x)
#   return(out)
# }
# x <- c(2,4,6,8,10)
# tryResult <- tryFunction(x)
# tryResult[[1]]
# rm(tryFunction,tryResult,x)
