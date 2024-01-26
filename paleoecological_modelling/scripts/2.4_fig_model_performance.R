library(tidyverse)
library(dismo)
library(cowplot)

# directories ####

targetdir <- './results/fig_model_performance'
if(!dir.exists(targetdir)){dir.create(targetdir)}

# load selected model ####

maxent_model <- readRDS('./work/maxent_model_selection/maxent_eval_select.RDS')

# Plot part A+C ####

## Plot feature importance (C) ####
View(maxent_model@results)

feature_importance <- maxent_model@results %>% as.data.frame() %>% rownames_to_column() %>%
  slice(7:22) %>% tidyr::extract(col=rowname,into=c('feature','indicator'),"([[:alnum:]]+).([[:alnum:]]+)",remove=T) %>%
  rename(value = V1) %>% mutate(feature = toupper(feature))
feature_importance$feature[feature_importance$feature=='ELEV'] <- 'Elevation'
env_lvls <- feature_importance %>% filter(indicator == 'permutation') %>%
  arrange(desc(value)) %>% pull(feature)
feature_importance <- feature_importance %>%
  mutate(feature=factor(feature,levels=env_lvls,ordered = T))
feature_importance$indicator[feature_importance$indicator == 'permutation'] <- 'permutation importance'

plt_feature_importance <- ggplot(data=feature_importance, aes(x=feature, y=value, fill=indicator))+
  geom_bar(stat='identity', position=position_dodge())+
  geom_text(aes(label=round(value)), vjust=0.5, hjust=1.1, color="black",
            position = position_dodge(0.9), size=2)+
  scale_fill_manual(limits=c('permutation importance','contribution'),
                    labels=c('Permutation importance','Contribution'),
                    values=c('#76b5c5','#eab676'))+
  labs(x='Feature', y="Percentage", fill='Measure')+
  scale_x_discrete(limits = rev(levels(feature_importance$feature)))+
  coord_flip()+
  theme_cowplot(font_size = 7, font_family = 'sans')+
  theme(text = element_text(size = 7, family='sans'))
plot(plt_feature_importance)

## Plot sample and background densities (A) ####
model_data <- bind_rows('Presence'=maxent_model@presence, 'Background'=maxent_model@absence, .id='type') %>% pivot_longer(-c(type), names_to = 'feature')
model_data$feature <- toupper(model_data$feature)
model_data$feature[model_data$feature == 'ELEV'] <- 'Elevation'
model_data <- model_data %>%
  mutate(feature=factor(feature,levels=env_lvls,ordered = T))

plt_dens <- ggplot(data=model_data, aes(x=value, fill=type))+
  geom_density(alpha=0.6)+
  facet_wrap(.~feature, scales='free', nrow=1)+
  scale_fill_manual(limits=c('Presence','Background'),
                    values=c('red','grey'))+
  labs(x='', y="Density", fill='')+
  theme_cowplot(font_size = 7, font_family = 'sans')+
  theme(text = element_text(size = 7, family='sans'),
        strip.text = element_text(size = rel(1)),
        strip.background = element_rect(fill = NA,colour = 'black',linewidth = 0.4),
        axis.text.x = element_text(size = 4))
plot(plt_dens)

## Combine C with free space for B ####
bottom_row <- plot_grid(
  NA, plt_feature_importance,
  labels = c('b', 'c'), rel_widths = c(1, 1.5), label_size = 7
)
plot(bottom_row)

## Export panel leaving free space for B ####
tiff(paste(targetdir,'fig_model_performance.tiff',sep = '/'),
     units="cm", width=18.3, height=10, res=600)
plot_grid(plt_dens, bottom_row, labels = c('a', ''), label_size = 7, ncol = 1, rel_heights = c(0.75,1))
dev.off()

# 3D scatterplot (B) to be added with a graphics program ####

library(plot3D)
presence_points <- maxent_model@presence
absence_points <- maxent_model@absence

# run with tiff and comment 1st lab row, then with svg an comment the 2nd
tiff(paste(targetdir,'fig_3dscatter.tiff',sep = '/'),units="cm",width=7.32,height=5.71,res=600,pointsize=7)
# svg(paste(targetdir,'fig_3dscatter.svg',sep = '/'),width=7.32,height=5.71,pointsize=7)
par(mar = c(0, 0, 0, 0))
plt_scatter <- scatter3D(x=absence_points$bio7, y=absence_points$bio8, z=absence_points$elev,
                         col = 'black',
                         pch = 20,
                         cex = 0.15,
                         colkey = FALSE,
                         bty = "b2",
                         ticktype = "detailed",
                         # xlab = 'BIO7', ylab='BIO8',zlab = "Elevation",
                         xlab = '', ylab='',zlab = "",
                         cex.lab=1, cex.axis=1, line=2,
                         alpha=.2)
plt_scatter <- scatter3D(x=presence_points$bio7, y=presence_points$bio8, z=presence_points$elev,
                         col = 'red',
                         pch = 20,0,
                         cex = 0.15,
                         add = TRUE)
dev.off()
