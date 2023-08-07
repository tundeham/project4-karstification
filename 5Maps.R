# Spatial distribution of recession parameters
# Create map plots of mean and Sd of parameters
# World Robinson map EPSG:54030
rm(list = ls())

# set work dir
wd <- "E:/Olarinoye/Project 4 - Spatial recession analysis/"
setwd(wd)

# load packages
library("rgdal")      # for spTransform() & project()
library(ggplot2)    # for ggplot()
library(ggpubr)
library(sp)
library(raster)
library(mapproj)
library(geosphere)

# load meta info and other data
load("metainfo+param.RData")
wokam <- readOGR(dsn="./Wokam_karst", layer="wokam_final", verbose=T) #wokam map
graticules.30 <- readOGR(dsn="./ne_10m_graticules_all", layer="ne_10m_graticules_30", verbose=T)

# Load ready to use data for creating global map from GitHub
load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
# This will load 6 objects:
#   xbl.X & lbl.Y are two data.frames that contain labels for graticule lines
#       They can be created with the code at this link: 
#       https://gist.github.com/valentinitnelav/8992f09b4c7e206d39d00e813d2bddb1
#   NE_box is a SpatialPolygonsDataFrame object and represents a bounding box for Earth 
#   NE_countries is a SpatialPolygonsDataFrame object representing countries 
#   NE_graticules is a SpatialLinesDataFrame object that represents 10 dg latitude lines and 20 dg longitude lines
#           (for creating graticules check also the graticule package or gridlines fun. from sp package)
#           (or check this gist: https://gist.github.com/valentinitnelav/a7871128d58097e9d227f7a04e00134f)
#   NE_places - SpatialPointsDataFrame with city and town points
#   NOTE: data downloaded from http://www.naturalearthdata.com/
#         here is a sample script how to download, unzip and read such shapefiles:
#         https://gist.github.com/valentinitnelav/a415f3fbfd90f72ea06b5411fb16df16

# Project from long-lat (unprojected data) to Robinson projection
# spTransform() is used for shapefiles and project() in the case of data frames
# for more PROJ.4 strings check the followings
#   http://proj4.org/projections/index.html
#   https://epsg.io/

# Create map in Robinson projection
PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
countries_rob  <- spTransform(NE_countries, CRSobj = PROJ)
graticules_rob <- spTransform(NE_graticules, CRSobj = PROJ)
graticules.30_rob <- spTransform(graticules.30, CRSobj = PROJ)
box_rob <- spTransform(NE_box, CRSobj = PROJ)
wokam_rob <- spTransform(wokam, CRSobj = PROJ)

# project long-lat coordinates for graticule label data frames 
# (two extra columns with projected XY are created)
prj.coord <- project(cbind(lbl.Y$lon, lbl.Y$lat), proj=PROJ)
lbl.Y.prj <- cbind(prj.coord, lbl.Y)
names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj")

prj.coord <- project(cbind(lbl.X$lon, lbl.X$lat), proj=PROJ)
lbl.X.prj <- cbind(prj.coord, lbl.X)
names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj")

# project coordinates of spring to robinson
metainfo$lon.proj <- as.numeric(metainfo$lon)
metainfo$lat.proj <- as.numeric(metainfo$lat)
lonlat.proj <- project(cbind(metainfo$lon.proj, metainfo$lat.proj), proj=PROJ)
metainfo[,13:14] <- lonlat.proj

# Plot layers
mid <- function(x){
  (max(x, na.rm=T)-min(x,na.rm=T))/2
}
#--- 1 ---
# Slow flow recession coefficent
a1 <- ggplot() +
  # add Natural Earth countries projected to Robinson, give black border and fill with gray
  geom_polygon(data=countries_rob, aes(long,lat, group=group), colour="white", fill="lightgrey", alpha = 0.8, size = 0.25) +
  geom_polygon(data=wokam_rob, aes(long,lat, group=group), colour="transparent", fill="powderblue", size = 0.25) +
  # Note: "Regions defined for each Polygons" warning has to do with fortify transformation. Might get deprecated in future!
  # alternatively, use use map_data(NE_countries) to transform to data frame and then use project() to change to desired projection.
  # add Natural Earth box projected to Robinson
  geom_polygon(data=box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  # add graticules projected to Robinson
  geom_path(data=graticules.30_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
  # add graticule labels - latitude and longitude
  # geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_fixed(ratio = 1) +
  # remove the background and default gridlines
  theme_void() +
  # add locations of spring discharge dataset, data = metaFile
  geom_point(data=metainfo, aes(x=lon.proj, y=lat.proj, colour=mean_alpB), size=0.25, na.rm=TRUE)+
  scale_colour_gradient2(low="blue", mid="#299617", high="red", midpoint=mid(metainfo$mean_alpB), guide="colourbar")+
  scale_x_continuous()+
  scale_y_continuous()+
  #ggtitle("(a) Mean estimated slow flow recession coefficient")+
  theme(plot.margin=unit(c(0,0,0,0),"cm"),
        legend.title=element_blank(), legend.key.height=unit(0.1,"cm"), 
        legend.key.width=unit(0.5, "cm"), legend.text=element_text(size=6), 
        legend.direction="horizontal", legend.position=c(0.5,0.15), legend.spacing.y=unit(0.1, "cm"))

# save to pdf and png file
# ggsave("./Results/Maps/slow coefficient.pdf", width=15, height=12, units="cm")
# ggsave(paste0(basePath,"map_plot/maps/map_draft2_NoS.png"), width=28, height=13.5, units="cm", dpi=600)

a2 <- ggplot() +
  # add Natural Earth countries projected to Robinson, give black border and fill with gray
  geom_polygon(data=countries_rob, aes(long,lat, group=group), colour="white", fill="lightgrey", alpha = 0.8, size = 0.25) +
  geom_polygon(data=wokam_rob, aes(long,lat, group=group), colour="transparent", fill="powderblue", size = 0.25) +
  # Note: "Regions defined for each Polygons" warning has to do with fortify transformation. Might get deprecated in future!
  # alternatively, use use map_data(NE_countries) to transform to data frame and then use project() to change to desired projection.
  # add Natural Earth box projected to Robinson
  geom_polygon(data=box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  # add graticules projected to Robinson
  geom_path(data=graticules.30_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
  # add graticule labels - latitude and longitude
  # geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_fixed(ratio = 1) +
  # remove the background and default gridlines
  theme_void() +
  # add locations of spring discharge dataset, data = metaFile
  geom_point(data=metainfo, aes(x=lon.proj, y=lat.proj, colour=sd_alpB), size=0.25, na.rm=TRUE)+
  scale_colour_gradient2(low="blue", mid="#299617", high="red", midpoint=mid(metainfo$sd_alpB), guide="colourbar")+
  scale_x_continuous()+
  scale_y_continuous()+
  # ggtitle("(a) Uncertainty of slow flow recession coefficient")+
  theme(plot.margin=unit(c(0,0,0,0),"cm"),
        legend.title=element_blank(), legend.key.height=unit(0.1,"cm"), 
        legend.key.width=unit(0.5, "cm"), legend.text=element_text(size=6), 
        legend.direction="horizontal", legend.position=c(0.5,0.15), legend.spacing.y=unit(0.1, "cm"))

# save to pdf and png file
# ggsave("./Results/Maps/slow coefficient.pdf", width=15, height=12, units="cm")

#--- 2 ---
# Quick flow recession coefficient
c1 <- ggplot() +
  # add Natural Earth countries projected to Robinson, give black border and fill with gray
  geom_polygon(data=countries_rob, aes(long,lat, group=group), colour="white", fill="lightgrey", alpha = 0.8, size = 0.25) +
  geom_polygon(data=wokam_rob, aes(long,lat, group=group), colour="transparent", fill="powderblue", size = 0.25) +
  # Note: "Regions defined for each Polygons" warning has to do with fortify transformation. Might get deprecated in future!
  # alternatively, use use map_data(NE_countries) to transform to data frame and then use project() to change to desired projection.
  # add Natural Earth box projected to Robinson
  geom_polygon(data=box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  # add graticules projected to Robinson
  geom_path(data=graticules.30_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
  # add graticule labels - latitude and longitude
  # geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_fixed(ratio = 1) +
  # remove the background and default gridlines
  theme_void() +
  # add locations of spring discharge dataset, data = metaFile
  geom_point(data=metainfo, aes(x=lon.proj, y=lat.proj, colour=mean_alpC), size=0.25, na.rm=TRUE)+
  scale_colour_gradient2(low="blue", mid="#299617", high="red", midpoint=mid(metainfo$mean_alpC), guide="colourbar")+
  scale_x_continuous()+
  scale_y_continuous()+
  #ggtitle("(a) Mean estimated quick flow recession coefficient")+
  theme(plot.margin=unit(c(0,0,0,0),"cm"),
        legend.title=element_blank(), legend.key.height=unit(0.1,"cm"), 
        legend.key.width=unit(0.5, "cm"), legend.text=element_text(size=6), 
        legend.direction="horizontal", legend.position=c(0.5,0.15), legend.spacing.y=unit(0.1, "cm"))
  
  # save to pdf and png file
  # ggsave("./Results/Maps/slow coefficient.pdf", width=15, height=12, units="cm")
  # ggsave(paste0(basePath,"map_plot/maps/map_draft2_NoS.png"), width=28, height=13.5, units="cm", dpi=600)
  
c2 <- ggplot() +
  # add Natural Earth countries projected to Robinson, give black border and fill with gray
  geom_polygon(data=countries_rob, aes(long,lat, group=group), colour="white", fill="lightgrey", alpha = 0.8, size = 0.25) +
  geom_polygon(data=wokam_rob, aes(long,lat, group=group), colour="transparent", fill="powderblue", size = 0.25) +
  # Note: "Regions defined for each Polygons" warning has to do with fortify transformation. Might get deprecated in future!
  # alternatively, use use map_data(NE_countries) to transform to data frame and then use project() to change to desired projection.
  # add Natural Earth box projected to Robinson
  geom_polygon(data=box_rob, aes(x=long, y=lat), colour="black", fill="transparent", size = 0.25) +
  # add graticules projected to Robinson
  geom_path(data=graticules.30_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
  # add graticule labels - latitude and longitude
  # geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size=2) +
  # the default, ratio = 1 in coord_fixed ensures that one unit on the x-axis is the same length as one unit on the y-axis
  coord_fixed(ratio = 1) +
  # remove the background and default gridlines
  theme_void() +
  # add locations of spring discharge dataset, data = metaFile
  geom_point(data=metainfo, aes(x=lon.proj, y=lat.proj, colour=sd_alpC), size=0.25, na.rm=TRUE)+
  scale_colour_gradient2(low="blue", mid="#299617", high="red", midpoint=mid(metainfo$sd_alpC), guide="colourbar")+
  scale_x_continuous()+
  scale_y_continuous()+
  # ggtitle("(a) Uncertainty of quick flow recession coefficient")+
  theme(plot.margin=unit(c(0,0,0,0),"cm"), 
        legend.title=element_blank(), legend.key.height=unit(0.1,"cm"), 
        legend.key.width=unit(0.5, "cm"), legend.text=element_text(size=6), 
        legend.direction="horizontal", legend.position=c(0.5,0.15), legend.spacing.y=unit(0.1, "cm"))

# save to pdf and png file
# ggsave("./Results/Maps/slow coefficient.pdf", width=15, height=12, units="cm")
cowplot::plot_grid(a1,c1,a2,c2, nrow=2, ncol=2, align="hv",
                   labels=c("(a) Mean of estimated slow flow recession coefficient",
                             "(b) Mean of estimated quick flow recession coefficent",
                            "(c) Uncertainty of estimated slow flow \nrecession coefficient",
                            "(d) Uncertainty of estimated quick flow \nrecession coefficent"),
                   label_size=6, hjust=0, vjust=c(1.5,1.5,1,1))
ggsave("./Results/Maps/slow coefficient.pdf", width=12, height=8, units="cm")

