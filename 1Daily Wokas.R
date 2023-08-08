# Porject 4 - script to select only spring discharge
# with daily or sub-daily observations from the WoKaS database

# set working dir
# rm(list = ls())
wd <- "D:/Project 4 - Spatial recession analysis/"
setwd(wd)

# select spring hydrographs with daily resolutions
# list all WoKaS hydrographs
# firstly, copy data from "WoKaS_Hydrograph_Datasets" to "Spring hydrographs" folder
wokas.dir <- "Data/Spring hydrographs/"
wokas.hydro <- list.files(paste0(wd, wokas.dir))

#initiate output to save hydrographs with dialy time step
wokas.daily <- list()
comp.list <- list()
ii <- 0

for(i in 1:length(wokas.hydro)){
  #i=400
  hydro <- read.csv(paste0(wd, wokas.dir, wokas.hydro[i]), sep="\t", skip=8, header=T)
  time1 <- as.Date(hydro[1,1], format=c("%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M"))
  time2 <- as.Date(hydro[nrow(hydro),1], format=c("%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M"))
  
  # remove invalid time format
  if(length(time1)>1){
    time1 <- time1[which(!is.na(time1))]
  }
  if(length(time2)>1){
    time2 <- time2[which(!is.na(time2))]
  }
  
  # retain one of duplicted time
  if(length(time1)>1){
    time1 <- time1[1]
  }
  if(length(time2)>1){
    time2 <- time2[1]
  }
 
  # format time to daily
  d.time1 <- as.Date(time1,"%d.%m.%Y")
  d.time2 <- as.Date(time2,"%d.%m.%Y")
  obs.span <- seq.Date(d.time1, d.time2,1)
  
  # check completeness of hydrograph
  comp <- nrow(hydro)/length(obs.span)
  comp.list[i] <- comp
  if(comp>0.51){
    ii = ii + 1
    wokas.daily[ii] <- wokas.hydro[i]
  }
}
# files to delete
del.file <- wokas.hydro[!(wokas.hydro %in% wokas.daily)]
# file.remove(paste0(wd,wokas.dir,del.file))







