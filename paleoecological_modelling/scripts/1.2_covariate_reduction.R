library(tidyverse)
library(raster)

# directories ####

covdir <- './work/covariates'
sampledir <- './work/sample'
targetdir <- './work/correlation'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# create full record object ####
# do not run it if you already have the full_record.csv file

rasterlist <- list.files(covdir, full.names=T)
rasterstack <- stack(rasterlist)
names(rasterstack) <- sub('[.](.*)$','',list.files(covdir,full.names=F))
variable_levels = c('bio1','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19','npp','elev')
interval_levels = c()
points <- rasterToPoints(rasterstack) %>% data.frame() %>%
  filter(complete.cases(.)) %>%  # Remove NAs
  mutate(id=as.numeric(rownames(.))) %>% # create a point ID
  pivot_longer(cols=-c(id,x,y),names_to='feat',values_to='val') %>%
  tidyr::extract(col=feat,into=c('variable','interval'),"([[:alnum:]]+)_([[:alnum:]]+)") %>% # separate variables and timeslices
  mutate(variable=factor(variable,levels = variable_levels)) %>%
  pivot_wider(names_from=variable,values_from = val)

write.csv(points,paste(sampledir,'full_record.csv',sep = '/'))
saveRDS(points,paste(sampledir,'full_record.rds',sep = '/'))

# load full record object (only if already present in ./work/sample)
# points <- read.csv('./work/sample/full_record.csv')
# points$X <- NULL

# Correlation ####

library(Hmisc)
library(corrplot)

## Pearson's correlation coefficients and p-values ####
cor1 <- rcorr(points %>% dplyr::select(-c(x,y,id,interval,npp)) %>% as.matrix(), type=c("pearson"))
write.csv(cor1$r,paste(targetdir,'pearson_coeff.csv',sep = '/'))
write.csv(cor1$p,paste(targetdir,'pearson_p-value.csv',sep = '/'))

## Pearson's correlation plot ####
tiff(paste(targetdir,'pearson_coeff.tiff',sep = '/'),units="cm",width=15,height=15,res=300)
corrplot(cor1$r, method = "color",type = "upper",order = "hclust",
         addCoef.col = "black",number.cex = .5, # add labels
         tl.col = "black",tl.srt = 90)
dev.off()

# Collinearity ####

library(fuzzySim)

collinearity_all <- multicol(vars=points %>% dplyr::select(-c(x,y,id,interval,npp)) %>% as.data.frame()) %>% mutate(feat=rownames(.))
collinearity_r90 <- multicol(vars=points %>% dplyr::select(-c(x,y,id,interval,npp,bio4,bio6,bio9,bio10,bio11,bio13,bio14,bio16)) %>% as.data.frame()) %>% mutate(feat=rownames(.))
collinearity_VIF5 <- multicol(vars=points %>% dplyr::select(bio5,bio7,bio8,bio12,bio15,bio17,bio18,bio19,elev) %>% as.data.frame()) %>% mutate(feat=rownames(.))
collinearity_VIF5 <- multicol(vars=points %>% dplyr::select(bio5,bio7,bio8,bio15,bio17,bio18,bio19,elev) %>% as.data.frame()) %>% mutate(feat=rownames(.))
# selected variables bio5 bio7 bio8 bio15 bio17 bio18 bio19 elev

collinearity_table <- collinearity_all %>% 
  left_join(collinearity_r90, by='feat') %>%
  left_join(collinearity_VIF5, by='feat') %>% 
  dplyr::rename(Rsquared.filter0=Rsquared.x, Tolerance.filter0=Tolerance.x,VIF.filter0=VIF.x,
                Rsquared.filter1=Rsquared.y, Tolerance.filter1=Tolerance.y,VIF.filter1=VIF.y,
                Rsquared.filter2=Rsquared, Tolerance.filter2=Tolerance,VIF.filter2=VIF)

write.csv(collinearity_table,paste(targetdir,'collinearity_table.csv',sep = '/'))
