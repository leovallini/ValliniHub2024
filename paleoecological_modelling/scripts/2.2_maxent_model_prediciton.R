library(tidyverse)
library(dismo)
library(raster)

set.seed(31)

# directories ####

targetdir <- './work/maxent_model_prediction'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# run prediction ####

# load model object
maxent_model <- readRDS('./work/maxent_model_selection/maxent_eval_select.RDS')

# Load predictors
rasterlist <- list.files('./work/covariates', full.names=F) %>% as.data.frame() %>% rename(name=1) %>% 
  tidyr::extract(col=name, into=c('variable','interval'), "([[:alnum:]]+)_([[:alnum:]]+)", remove=F)

list_variables <- unique(rasterlist$variable)
list_intervals <- unique(rasterlist$interval)

# predict all intervals and produce raster files for mean and sd value
for(interval_subset in list_intervals){
  print(interval_subset)
  predictors <- rasterlist %>% filter(interval==interval_subset) %>% 
    pull(name) %>% paste('./work/covariates',., sep='/') %>% stack()
  names(predictors) <- rasterlist %>% filter(interval==interval_subset) %>% 
    pull(name) %>% as.data.frame %>% tidyr::extract(1, "([[:alnum:]]+)_[[:alnum:]]+") %>% pull()
  predictors
  prediction <- predict(maxent_model, predictors)
  writeRaster(mean(prediction),paste(targetdir,'/pred_mean_',interval_subset,'.tif',sep=''),overwrite=TRUE)
  writeRaster(calc(prediction,sd),paste(targetdir,'/pred_std_',interval_subset,'.tif',sep=''),overwrite=TRUE)
}

# predict on observation sites
sites <- readRDS('./work/sample/sample_points.rds')
sites$prediction <- predict(maxent_model,sites)
sites_pred <- sites@data %>% group_by(ID) %>% summarise(prediction=mean(prediction))
sites <- read.table('./data/data4model.txt',header = T) %>% left_join(sites_pred)

write.csv(sites,paste(targetdir,'sites_with_prediction.csv',sep = '/'),row.names = F)
