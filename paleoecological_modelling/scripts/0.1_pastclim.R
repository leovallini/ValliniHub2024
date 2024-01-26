library(raster)
library(dplyr)
devtools::install_github("EvolEcolGroup/pastclim")
library(pastclim)
library(sf)
library(sp)
library(rgdal)

# variables and directories ####

studyarea <- extent(-10,140,0,75) # study area from Vallini L.
if(!dir.exists('./work')){dir.create('./work')}
sourcedir <- './work/pastclim'
if(!dir.exists(sourcedir)){dir.create(sourcedir)}
targetdir <- './work/covariates'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# download data ####

# careful! it is a 1.5GB file
download_dataset(dataset = 'Beyer2020',path_to_nc = sourcedir)
download_dataset(dataset = 'Beyer2020',
                 bio_variables = c('altitude'),
                 path_to_nc = sourcedir)

# check downloaded data
get_downloaded_datasets()
get_vars_for_dataset('Beyer2020')
# variable names is in block letters

# load just one variable and check crs, extension, etc.
dummy <- brick('work/pastclim/Beyer2020_all_vars_v1.0.0.nc',var='BIO1')
dummy
plot(dummy)
# time scale is:
# - every 2k years from 120kya to 22kya
# - every 1k years from 22kya to present

# load, crop and re-project data ####

# set projection
ESRI53031 <- st_crs('ESRI:53031')
# this projection is an equidistant projection good for Asia and Europe

# define target time frames you want; if you want all just use:
# timeslices <- names(dummy)
timeslices <- paste("X.",seq(70000,30000,-2000),sep="")

## load, crop and re-project bio data and then write rasters ####
covariates <- c('BIO1',paste('BIO',4:19,sep = ''))
covariates.new <- tolower(covariates)
nc_vars <- list.files(path = sourcedir,pattern = 'vars',full.names = T)
for(i in 1:length(covariates)){
  iras <- brick(nc_vars,var=covariates[i]) # load raster
  iras <- crop(iras,studyarea) # crop to study area
  iras <- projectRaster(iras,crs = ESRI53031$wkt) # re-project raster
  for(j in timeslices){
    tras <- iras[[j]] # select layer (i.e., time)
    tname <- paste(substr(j,3,(nchar(j)-3)),'kya',sep = '') # get time name
    rasname <- paste(paste(covariates.new[i],tname,sep = '_'),'.tif',sep = '') # make a name
    writeRaster(tras,filename=paste(targetdir,rasname,sep = '/'),
                format="GTiff",overwrite=TRUE) # save new raster
  }
}
rm(i,j,iras,tras,tname,rasname)

## load, crop and re-project npp data and then write rasters ####
iras <- brick(nc_vars,var='npp') # load raster
iras <- crop(iras,studyarea) # crop to study area
iras <- projectRaster(iras,crs = ESRI53031$wkt) # re-project raster
for(j in timeslices){
  tras <- iras[[j]] # select layer (i.e., time)
  tname <- paste(substr(j,3,(nchar(j)-3)),'kya',sep = '') # get time name
  rasname <- paste(paste('npp',tname,sep = '_'),'.tif',sep = '') # make a name
  writeRaster(tras,filename=paste(targetdir,rasname,sep = '/'),
              format="GTiff",overwrite=TRUE) # save new raster
}
rm(j,iras,tras,tname,rasname)

## load, crop and re-project elevation data and then write rasters ####
nc_vars <- list.files(path = sourcedir,pattern = 'top',full.names = T)
iras <- brick(nc_vars,var='altitude') # load raster
iras # time has a different name here
names(iras) <- names(dummy)
iras <- crop(iras,studyarea) # crop to study area
iras <- projectRaster(iras,crs = ESRI53031$wkt) # re-project raster
for(j in timeslices){
  tras <- iras[[j]] # select layer (i.e., time)
  tname <- paste(substr(j,3,(nchar(j)-3)),'kya',sep = '') # get time name
  rasname <- paste(paste('elev',tname,sep = '_'),'.tif',sep = '') # make a name
  writeRaster(tras,filename=paste(targetdir,rasname,sep = '/'),
              format="GTiff",overwrite=TRUE) # save new raster
}
rm(j,iras,tras,tname,rasname)
