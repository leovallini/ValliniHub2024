library(tidyverse)
library(cowplot)
library(sf)
library(sp)
library(raster)
library(spatialEco)
library(ks)
set.seed(31)

# directories ####

if(!dir.exists('./results')){dir.create('./results')}
targetdir <- './work/sample/background_sample_validation'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# load full record file ####

dataset <- readRDS('./work/sample/full_record.RDS')
head(dataset)
summary(dataset[-c(1:4)])

# produce probability surface for background sampling based on effort ####

# load sampling effort
effort <- read.table('./data/TimmermanNature2022_ST1.csv',sep = ';',header = T)
head(effort)
# keep only homo sapiens data with average age equal or lower than 70ky
unique(effort$species) # there is both "Homo sapiens" and "Homo_sapiens"
effort <- effort %>% filter(species %in% c("Homo sapiens","Homo_sapiens")) %>%
  mutate(avg.age = (max.age+min.age)/2) %>% filter(avg.age <= 70000)
effort.xy <- effort[3:4]
effort.xy$count <- 1
# make a spatial object
ESRI53031 <- st_crs('ESRI:53031') # New CRS
EPSG4326 <- st_crs(4326)  # Original CRS
coordinates(effort.xy) <- ~long+lat
crs(effort.xy) <- CRS(EPSG4326$wkt)
effort.xy <- effort.xy %>% spTransform(CRS(ESRI53031$wkt))
plot(effort.xy)
# load model raster
dummy <- raster('./work/covariates/bio1_30kya.tif')
effort.r <- raster::rasterize(effort.xy,dummy,field='count',fun=sum,background=0)
plot(effort.r)
effort.xy <- rasterToPoints(effort.r) %>% data.frame() %>% filter(layer != 0)
# scale weights so that sum of all sampling weights = total sample size
eff.scale <- length(effort.xy[,3])/sum(effort.xy[,3])
effort.xy <- effort.xy %>% rename(nVisits = layer) %>%
  mutate(Weights = nVisits*eff.scale)
head(effort.xy)
# Do a 2d kernel density estimate (using scaled sampling effort as weighting) 
# to create a probability surface of sampling effort.
eff.density <- kde(effort.xy[1:2],w = effort.xy$Weights)
# now rasterize it
effort.r <- raster(eff.density)
crs(effort.r) <- CRS(ESRI53031$wkt)
plot(effort.r)
# make mask from background points
bg.mask <- unique(dataset[1:2])
coordinates(bg.mask) <- ~x+y
crs(bg.mask) <- CRS(ESRI53031$wkt)
bg.mask <- rasterize(bg.mask,dummy)
plot(bg.mask)
# Clip sample effort to same extent and resolution as background
effort.r <- resample(effort.r,bg.mask,method='bilinear')
# mask sample effort to remove outside area
effort.r <- mask(effort.r,bg.mask)
plot(effort.r)
# Normalize bias file
effort.r <- terra::rast(effort.r) %>% raster.transformation(trans="norm")
plot(effort.r)
# save bias file
writeRaster(effort.r,'./work/sample/effort_bias.tif',overwrite=T)

# clean up env
rm(effort,effort.xy,ESRI53031,EPSG4326,dummy,eff.scale,eff.density,bg.mask)

# Test the effect of sample size ####

d1 <- dataset %>% dplyr::select(c(bio5,bio7,bio8,bio15,bio17,bio18,bio19,elev))
N_mean <- d1 %>% colMeans()
sample_error <- setNames(rep(NA, times=9), c('bio5','bio7','bio8','bio15','bio17','bio18','bio19','elev','n'))
wgts <- raster::extract(effort.r$layer,dataset[,1:2])[,2]

for(samplesize in c(1,10,100,1000,10000,100000,500000)){
  for(i in seq(1,1000)){
    n_mean <- d1 %>% slice_sample(n = samplesize,weight_by = wgts) %>% colMeans()
    nN_Diff <- abs(N_mean - n_mean)
    nN_Diff['n'] <- samplesize
    sample_error <- bind_rows(sample_error, nN_Diff)
  }
}
sample_error <- sample_error %>% drop_na() %>%  mutate(n=as.character(n))
saveRDS(sample_error,paste(targetdir,'background_sample_size_error.rds',sep = '/'))

# plot sample error distribution by sample size
plt_sample_error <- ggplot(sample_error %>% drop_na() %>% pivot_longer(-n) %>% filter(n != '5e+05'), aes(x=n, y=value))+
  facet_grid(name~., scale='free')+
  geom_boxplot()+
  labs(
    title='Effect of background sample size',
    x='Sample size (n)',
    y='Absolute sample error')+
  scale_x_discrete(labels = c('1','10','100','1,000','10,000','100,000'))+
  theme_bw()
plot(plt_sample_error)

tiff(paste(targetdir,'background_sample_size_error.tiff',sep = '/'),
     units="cm",width=9,height=18,res=600)
plt_sample_error
dev.off()

# Test random vs. stratified vs. effort-based background sampling ####

d2 <- dataset %>%
  dplyr::select(c(interval,bio5,bio7,bio8,bio15,bio17,bio18,bio19,elev))

## random ####
# This is a random sample over all time slices
sample_2.random <- d2 %>% sample_n(10000)

## stratified ####
sample <- read.table('./data/data4model.txt',header = T)
sample <- sample %>%
  mutate(max2=ceiling(max.age/2000)*2,min2=floor(min.age/2000)*2,
         avg.age=round((max.age+min.age)/2),
         periods=map2(map2(max2,min2,seq,by=-2),'kya',paste,sep = ''),
         nperiod=unlist(map(periods,length)),fract=1/nperiod) %>%
  unnest(cols = c(periods))
sample %>% group_by(periods) %>% summarise(nsites=sum(fract)) %>%
  as.data.frame()

# Here sampling is stratified by the number of arch. sites per time slice
sample_2.stratified <- bind_rows('30kya'= d2 %>% filter(interval=="30kya") %>% sample_n(1.42/140*10000),
                                 '32kya'= d2 %>% filter(interval=="32kya") %>% sample_n(1.92/140*10000),
                                 '34kya'= d2 %>% filter(interval=="34kya") %>% sample_n(11.17/140*10000),
                                 '36kya'= d2 %>% filter(interval=="36kya") %>% sample_n(23.75/140*10000),
                                 '38kya'= d2 %>% filter(interval=="38kya") %>% sample_n(32/140*10000),
                                 '40kya'= d2 %>% filter(interval=="40kya") %>% sample_n(34/140*10000),
                                 '42kya'= d2 %>% filter(interval=="42kya") %>% sample_n(21.25/140*10000),
                                 '44kya'= d2 %>% filter(interval=="44kya") %>% sample_n(9.83/140*10000),
                                 '46kya'= d2 %>% filter(interval=="46kya") %>% sample_n(3.83/140*10000),
                                 '48kya'= d2 %>% filter(interval=="48kya") %>% sample_n(0.83/140*10000)
)

## effort-based ####
sample_2.effort <- d2 %>% slice_sample(n = 10000,weight_by = wgts)

## plot comparison ####
sample_2 <- bind_rows('uniform'=sample_2.random,
                      'stratified'=sample_2.stratified,
                      'effort-weighted'=sample_2.effort,.id='strategy')
saveRDS(sample_2,paste(targetdir,'background_sample_strategy.rds',sep = '/'))
x <- sample_2 %>% pivot_longer(-c(strategy, interval)) %>%
  mutate(strategy = factor(strategy,ordered = T,levels = c('uniform','stratified','effort-weighted')))

plt_sample2 <- ggplot(x, aes(x=strategy, y=value))+
  facet_grid(name~., scale='free')+
  geom_violin()+
  labs(
    title='Effect of sampling strategy',
    x='Sampling type',
    y='Value')+
  theme_bw()
plot(plt_sample2)

tiff(paste(targetdir,'background_sample_strategy.tiff',sep = '/'),
     units="cm", width=9, height=18, res=600)
plt_sample2
dev.off()

## produce figures for manuscript ####
tiff("./results/fig_background_sample_validation.tiff",
     units="cm", width=18.3, height=18, res=600)
plot_grid(plt_sample_error,plt_sample2, labels = c('a','b'), rel_widths = c(1,0.65), align='h')
dev.off()

ggsave("./results/fig_background_sample_validation.eps",
       units="cm", width=18.3, height=18, dpi=600)

# Draw and save the background sample ####

set.seed(31)

background_sample <- dataset %>% slice_sample(n = 10000,weight_by = wgts)
saveRDS(background_sample, './work/sample/background_sample_points.rds')

save.image(paste(targetdir,'bg.RData',sep = '/'))
