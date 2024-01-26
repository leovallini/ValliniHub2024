library(dismo)
library(raster)

set.seed(31)

# directories ####

targetdir <- './work/maxent_model_classification'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# load models ####

maxent_model <- readRDS('./work/maxent_model_selection/maxent_eval_select.RDS')
sites <- read.csv('./work/maxent_model_prediction/sites_with_prediction.csv')
rasterlist <- list.files('./work/maxent_model_prediction',full.names=T,pattern='pred_mean')

# Reclassify ####

# 5% percentile of predicted training sites is used as threshold
thrs <- quantile(sites$prediction,0.05)
for(r in rasterlist){
  basename(r)
  rfile <- raster(r) > thrs
  writeRaster(rfile,paste(targetdir,'/range_',basename(r),sep=''),overwrite=TRUE)
}

