# Random forest model for prediction recession characteritics and 
# degree of karstification using catchment attributes
rm(list = ls())

# set work dir
wd <- "E:/Olarinoye/Project 4 - Spatial recession analysis/"
setwd(wd)

# load packages
pkg <- list("randomForest","rsample","ggplot2","ranger","caret","maps")
lapply(pkg, function(x)
  if(!(x %in% installed.packages())){
    install.packages(x)
  })
library(hydroGOF)
library(randomForest)
library(rsample)
library(ggplot2)
library(ranger)
library(caret)
library(maps)

# load meta info and other data
load("metainfo+param.RData")

#--- Pre-analysis: data cleaning and sorting ---

# load catchment attributes
attr.dir <- "./Data/Catchment_attributes/"
catch.attr <- list.files(attr.dir)  # folder contains 6 files for different spring areas
out.dist <- catch.attr[which(catch.attr == "out_dist.txt")]
out.attr <- read.table(paste0(attr.dir,out.dist), sep=",", header=T)
out.attr$id <-  gsub("-a", "(A)", out.attr$id); gsub("-b", "(B)", out.attr$id)

# merge df to get springs with attributes info
metainfo2 <- merge(metainfo, out.attr, by="id")
no.attr <- metainfo$id[which(!(metainfo$id %in% out.attr$id))]
write.table(no.attr, file="./noAttr_id.txt", row.names=F)
rm(no.attr, out.dist)

# plot spring for spatial overview
pdf("./Results/Maps/springs_with_attr.pdf", width=7.5, height=5)
map(wrap=T, col="grey")
points(metainfo2$lon,metainfo2$lat,pch=19,cex=0.7,col=1)
title(main= sprintf("springs (%s) with catchment attributes",nrow(metainfo2)))
box()
dev.off()

#--- 1. ANALYSIS: Basic RF model ---
# names of predictors: catchment attributes
attr.list <- c("HI","P","Psi","PET","PETsi","CORR","TA","ELEV")

set.seed(600)
attr.split <- initial_split(metainfo2, prop=0.7)
attr.train <- training(attr.split)
attr.test <- testing(attr.split)

no.train <- nrow(attr.train)
no.test <- nrow(attr.test)

# RF model - Basic model ----
set.seed(600)

# sample colums of df with predictor attributes and outcome "mean_ki"
attr.train <- attr.train[,c("mean_ki",attr.list)]
attr.train <- attr.train[!(attr.train$mean_ki %in% c(NaN)),]

attr.test <- attr.test[,c("mean_ki",attr.list)]
attr.test <- attr.test[!(attr.test$mean_ki %in% c(NaN)),]

# default RF
m1 <- randomForest(formula=mean_ki ~ ., data=attr.train, 
                   importance=T, proximity=T, oob.prox=F)
pred.m1 <- predict(m1, attr.test)

# make useful plots
plot(m1)

reg.train <- lm(attr.train$mean_ki ~ m1$predicted)
reg.test <- lm(attr.test$mean_ki ~ pred.m1)

# calculate error metrics between obs and predicted
r2.train <- round(summary(reg.train)$r.squared,2)
r2.test <- round(summary(reg.test)$r.squared,2)

rmse.train <- round(rmse(attr.train$mean_ki, m1$predicted),2)
rmse.test <- round(rmse(attr.test$mean_ki, pred.m1),2)

# find max value for plot axis
lim <- round(max(c(attr.train$mean_ki, m1$predicted)),2)


# save plot
pdf("./Results/Regression_analysis/obsVpre-k-defaultRF.pdf", width=5, height=5)
par(pty="s", las=1, cex.axis=0.8)
plot(attr.train$mean_ki, m1$predicted, xlim=c(0,lim),ylim=c(0,lim), pch=19, cex=0.75,
     main=sprintf("Validation of %s data point with %s training data sample", no.test, no.train), 
     cex.main=0.75, xlab="Observed mean ki", ylab="Predicted ki")
points(attr.test$mean_ki, pred.m1, pch=19, cex=0.75, col="grey")
text(x=lim-0.25, y=lim-0.01, labels=paste0("R2 Train = ", r2.train, "\nRMSE Train = ", rmse.train), cex=0.75)
text(x=lim-0.1, y=lim-0.01, labels=paste0("R2 Val = ", r2.test, "\nRMSE Val = ", rmse.test), cex=0.75)
abline(reg.train)
abline(reg.test, col="grey")
abline(0,1,lty=2,col="grey")
box()
dev.off()

# RMSE of trres with lowest MSE
sqrt(m1$mse[which.min(m1$mse)])


#--- 2. ANALYSIS: RF model - Advanced and tuning ----
# tuning to find optima number of variable for each split
set.seed(601)
m2 <- tuneRF(
  x          = attr.train[attr.list],
  y          = attr.train$mean_ki,
  ntreeTry   = 500,
  mtryStart  = 4,
  stepFactor = 1.5,
  improve    = 0.01,
  trace      = FALSE      # to not show real-time 
)

#--- grid search for best RF parameter search ---
hyper.grid <- expand.grid(
  mtry = seq(2, length(attr.list), by=1),
  node_size = seq(3, 9, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)
# total number of combinations
nrow(hyper.grid)

for(i in 1:nrow(hyper.grid)){
  
  # train model
  model <- ranger(
    formula = mean_ki ~ ., 
    data = attr.train,
    num.trees = 500,
    mtry = hyper.grid$mtry[i],
    min.node.size = hyper.grid$node_size[i],
    sample.fraction = hyper.grid$sample_size[i],
    seed = 123
  )
  # add OOB error to grid
  hyper.grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
}

set <- hyper.grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

# let repeat RF with best parameter settings
OOB_RMSE <- vector(mode = "numeric", length = 100)

for(i in seq_along(OOB_RMSE)){
  opt.ranger <- ranger(
    formula = mean_ki ~ ., 
    data = attr.train,
    num.trees = 500,
    mtry = set$mtry[1],
    min.node.size = set$node_size[1],
    sample.fraction = set$sample_size[1],
    importance = "impurity"
  )
  OOB_RMSE[i] <- sqrt(opt.ranger$prediction.error)
}
# histogram of OOB error for random boostrapping
hist(OOB_RMSE, breaks = 20)

imp.df <- data.frame(opt.ranger$variable.importance)
imp.df.row <- row.names(imp.df)
ord <- order(imp.df[,1], decreasing = T)

imp.df <- imp.df[ord,]
imp.df.row <- imp.df.row[ord]

imp.df <- as.data.frame(cbind(imp.df.row,imp.df))
colnames(imp.df) <- c("variables", "x")

imp.df[(1:length(attr.list)), ] %>%
  ggplot(aes(reorder(variables, as.numeric(x)), as.numeric(x)))+
  geom_col() +
  xlab("Variables") +
  ylab("x") +
  coord_flip() +
  ggtitle("Importance of catchment attributes")
ggsave(filename = "./Results/Regression_analysis/attributes-importance.pdf")
