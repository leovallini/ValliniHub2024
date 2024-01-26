library(tidyverse)
library(dismo)
options(java.parameters = "-Xmx20g" )
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
library(ENMeval) #v2.0.3

set.seed(31)

# variables and directories ####

targetdir <- './work/maxent_model_selection'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# Covariates selected in script 1.2_covariates_reduction.R
covariates <- c('x','y','bio5','bio7','bio8','bio15','bio17','bio18','bio19','elev')

# Evaluate combinations of FCs and RMs ####

## load data and set type, where 0=background, 1=sample, as indicated in dismo::maxent ####
# presence sample
sample <- readRDS('./work/sample/sample_points.rds')@data %>% 
  as.data.frame() %>% 
  drop_na(bio1:npp) %>% 
  mutate(type=1)
occurrence <- sample %>% dplyr::select(all_of(covariates))
# background sample
background <- readRDS('./work/sample/background_sample_points.rds') %>% 
  mutate(type=0) %>% 
  dplyr::select(all_of(covariates))
# convert sites IDs to avoid errors in the maxent run
sample$ID <- sample %>% pull(ID) %>% as.factor() %>% as.numeric()
# define grouped presence samples
user.grp <- list(occs.grp = sample %>% pull(ID), # Make sure that NO number is missing
                 bg.grp = rep(unique(sample$ID), length.out=nrow(background)))

## Evaluate models ####
maxent_eval <- ENMevaluate(occs=occurrence,
                           bg=background,
                           algorithm='maxent.jar',
                           partitions='user',
                           user.grp=user.grp,
                           tune.args = list(
                             fc=c('L','LQ','LQH', 'LQP', 'LQHP', 'LQHPT'),
                             rm=c(0.2,0.4,0.6,0.8,1,1.5,2,3,4)),
                           parallel=T)

View(maxent_eval@results)
saveRDS(maxent_eval,paste(targetdir,'maxent_eval.RDS',sep = '/'))

## Plot results ####
# plot model evaluation
library(RColorBrewer)

evaldata <- maxent_eval@results %>% pivot_longer(-c(fc:tune.args), names_to='criterion', values_to='val') %>%
  filter(criterion %in% c('delta.AICc','or.10p.avg','auc.val.avg')) %>% 
  mutate(fc=factor(fc, levels=c('L','LQ','LQH', 'LQP', 'LQHP', 'LQHPT')),
         rm=as.numeric(as.character(rm)),
         criterion=factor(criterion, levels=c('delta.AICc','or.10p.avg','auc.val.avg'), labels=c('Delta AICc','OR .10','AUC test')))

evalplot <- ggplot(data=evaldata, aes(x=rm, y=val, group=fc, color=fc))+
  geom_vline(xintercept=1, linetype='dashed')+
  geom_point()+
  geom_line()+
  facet_wrap(~criterion,
             scales='free_y',
             strip.position ='top')+
  labs(x='Regularization multiplier', y='', color='Feature class')+
  scale_x_continuous(limits=c(0,4))+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(
    strip.text = element_text(size = rel(1)),
    strip.background = element_rect(fill = NA, colour = 'black', linewidth = 0.4)
  )
plot(evalplot)

tiff("./results/fig_model_evaluation.tiff",units="cm",width=19,height=6,res=600)
evalplot
dev.off()

## Export best results (after investigating the figure above) ####

selection <- maxent_eval@results %>%
  dplyr::select(tune.args,delta.AICc,or.10p.avg,auc.val.avg) %>%
  arrange(delta.AICc) %>% mutate(rank1=1:54) %>%
  arrange(or.10p.avg) %>% mutate(rank2=1:54) %>%
  arrange(desc(auc.val.avg)) %>% mutate(rank3=1:54) %>%
  mutate(totrank=rank1+rank2+rank3) %>% arrange(totrank)

maxent_eval_select <- maxent_eval@models$fc.LQP_rm.0.2

saveRDS(maxent_eval@models$fc.LQP_rm.0.2,
        paste(targetdir,'maxent_eval_select.RDS',sep = '/'))
