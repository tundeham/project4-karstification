# Porject 4 - Extract recession events of the spring discharge
# observations with daily time step measurements
rm(list = ls())

# set working dir
library("tidyr")
wd <- "E:/Olarinoye/Project 4 - Spatial recession analysis/"
setwd(wd)
source("./Source modules/REM.R")
wokas.dir <- "Data/Spring hydrographs/"
wokas.hydro <- list.files(paste0(wd, wokas.dir))

#--- 1 ---
# initiate output data 
spring.data <- list()
metainfo <- data.frame(matrix(data=NA, nrow=length(wokas.hydro), ncol=3))
colnames(metainfo) <- c("id","lon","lat")

for(i in 1:length(wokas.hydro)){
  # 69,74,82 - index of data with different sep value
  # i=212
  hydro <- wokas.hydro[i]
  id <- unlist(strsplit(hydro, "@"))[1]
  metainfo[i,1] <- id
  
  # aggregate sub-daily obs to daily
  h <- read.csv(paste0(wd,wokas.dir,hydro), sep="\t", dec=".", header=T, skip=8)
  h$date <- as.Date(h$Timestamp, "%d.%m.%Y %H:%M:%S")
  h$date <- format(h$date, "%Y-%m-%d")
  h <- aggregate(discharge ~ date, data = h, FUN = "mean")
  
  # check if dates are complete
  # and fill missing date rows
  h$date <- as.Date(h$date)
  h <- h %>%
    complete(date = seq(min(date), max(date), by = 'day')) 
  
  spring.data[[i]] <- h

}

# complete metainfo dataframe
wokas.meta <- readxl::read_xlsx(paste0(wd,"Data/WoKaS_Hydrograph_Metafile/WoKaS_Hydrograph_Metafile_updated.xlsx"))
for(i in 1:nrow(metainfo)){
  # i=1
  ind <- which(metainfo$id[i]==wokas.meta$WoKaS_ID)
  metainfo$lon[i] <- wokas.meta$Longitude[ind]
  metainfo$lat[i] <- wokas.meta$Latitude[ind]
  metainfo$name[i] <- wokas.meta$Name[ind]
}
# save this variables to a workspace
# save(file="SpringHydrographs.RData",spring.data,metainfo)

#--- 2 ---
# load data
load("SpringHydrographs.RData") # treated spring timeseries

# recession extraction method
models <- "Vogel"
res.list <- list()

for(i in 1:length(spring.data)){
  # 212 - index of data with error message
  # i=212
  h <- data.frame(spring.data[i])
  ex_model <- models
  # nm <- paste0(spring_id, "@", ext_model)
  res <- rec.extract(h, model=ex_model, par1=NULL, par2=NULL, CoV=0.20, fit=0.90, len=15, plot=F)
  res.list[[i]] <- res
}
# save extracted recession events to a workspace
save(file="RecessionEvents.RData", res.list)

# find spring with less than two recession event
n0 <- lapply(res.list, function(x) 
  length(x)<2)
n0 <- which(n0==TRUE)

# extract event using less days
for(i in 1:length(n0)){
  i=1
  cat("...Processing spring no",i,"\n")
  j <- n0[[i]]
  hj <- data.frame(spring.data[[j]])
  ex_model <- models
  res <- rec.extract(hj, model=ex_model, par1=NULL, par2=NULL, CoV=0.20, fit=0.90, len=15, plot=T)
  res.list[[j]] <- res
}
save(file="RecessionEvents.RData", res.list) # save data again

# load data
load("RecessionEvents.RData") # extracted recessions from all springs



# code check
# n <- 212
# hh <- spring.data[[n]]
# plot(hh$date,hh$discharge, type="l")
# metainfo[n,]
# res.list[[212]]
