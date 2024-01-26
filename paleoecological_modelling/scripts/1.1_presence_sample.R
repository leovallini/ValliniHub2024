library(tidyverse)
library(readxl)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(spdplyr)

# directories ####

sourcedir <- './work/covariates'
targetdir <- './work/sample'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# Coordinate Systems ####

ESRI53031 <- st_crs('ESRI:53031') # New CRS
EPSG4326 <- st_crs(4326)  # Original CRS

# Load arch. sites ####

sites <- read.table('./data/data4model.txt',header = T)
sites <- sites %>%
  mutate(max2=ceiling(max.age/2000)*2,min2=floor(min.age/2000)*2,
         avg.age=round((max.age+min.age)/2),
         periods=map2(map2(max2,min2,seq,by=-2),'kya',paste,sep = ''),
         nperiod=unlist(map(periods,length))) %>%
  unnest(cols = c(periods))
head(sites)
coordinates(sites) <- ~long+lat
crs(sites) <- CRS(EPSG4326$wkt)
sites <- sites %>% spTransform(CRS(ESRI53031$wkt))
plot(sites)

writeOGR(sites, dsn=targetdir, layer='arch_sites', driver = "ESRI Shapefile")

# Create Sample ####

# Create an empty raster with the template layout
dummy <- raster() %>%
  projectRaster(to=raster('./work/covariates/bio1_30kya.tif')) %>% 
  rasterToPoints(spatial=T)
# plot(dummy)

# Construct sample points buffering each coordinates by 50km
buffer <- buffer(sites, width=50000, dissolve=F)
plot(buffer)
sample <- raster::intersect(dummy, buffer)
plot(sample)
table(sample$ID)

sample@data <- sample@data %>% mutate(x=as.numeric(sample@coords[,'x']), y=as.numeric(sample@coords[,'y']), id=rownames(.)) %>% 
  tidyr::extract(col=id, into=c('id_sub'), "[[:alnum:]]+[[:punct:]]([[:alnum:]]+)", remove=F) %>%
  mutate(id_sub=ifelse(is.na(id_sub), 0, id_sub))

# Get environmental values and rename layers in the stack as files
rasterstack <- list.files(sourcedir, full.names=T) %>% stack()
names(rasterstack) <- sub('[.](.*)$','',list.files(sourcedir,full.names=F))
rasterstack

# Extract environmental values for each point at each time step
sample2 <- raster::extract(rasterstack,sample,sp=T)@data %>% 
  pivot_longer(cols=-c(ID:id_sub),names_to='feat',values_to='val') %>% 
  tidyr::extract(col=feat,into=c('variable','ts'),"([[:alnum:]]+)_([[:alnum:]]+)") %>% 
  pivot_wider(id_cols=c(id, ts),names_from = variable, values_from = val)

# Join environmental data to sample data.frame keeping only info for the target time steps
sample3 <- sample %>% left_join(sample2,by=c('id'='id','periods'='ts')) %>% 
  data.frame() %>% drop_na(bio1:npp)

# Average data aged between two time steps
sample3.1 <- sample3 %>% filter(nperiod==2) %>%
  mutate(age.dev = abs(avg.age-1000*as.numeric(substr(periods,1,2)))) %>%
  group_by(ID,id_sub,species,x,y,max.age,min.age,avg.age,n.dates) %>%
  summarize(periods=paste(periods,collapse = ''),
            across(bio1:npp,weighted.mean,w=age.dev))

# Merge data
sample4 <- sample3 %>% filter(nperiod != 2) %>%
  dplyr::select(-optional,-id,-nperiod,-min2,-max2,-x.1,-y.1,-age.range) %>%
  full_join(sample3.1)

# transform to sp object and save
sample4 <- cbind(sample4$x,sample4$y) %>% SpatialPointsDataFrame(sample4)
crs(sample4) <- CRS(ESRI53031$wkt)

saveRDS(sample4,paste(targetdir,'sample_points.rds',sep = '/'))
writeOGR(sample4,dsn=targetdir,layer='sample_points',driver = "ESRI Shapefile")
