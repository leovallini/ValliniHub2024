library(raster)

# read rasters with range and npp ####

rasterlist_range <- list.files('./work/maxent_model_classification',full.names=T,pattern='range')
rasterlist_npp <- list.files('./work/covariates',full.names=T,pattern='npp')
dates = seq(30,70,2)
range <- stack()
j=1
pdf('range.pdf', height = 10, width = 10)
  for(i in rasterlist_range) { 
    range <- addLayer(range, raster(i))
    terra::plot(raster(i), axes = F, col = c("grey80", "forestgreen"), legend = F, main = paste(dates[j], "kya - Range")) 
    j=j+1
  }
dev.off()
npp <- stack()
for(i in rasterlist_npp){ npp <- addLayer(npp, raster(i)) }

# plot range all dates ####

# identify cells that are always sea (always NA)
range_sea <- stack()
for(i in rasterlist_range){
  rasterlayer = raster(i)
  rasterlayer[is.na(rasterlayer)] <- 100
  range_sea <- addLayer(range_sea, rasterlayer) 
}
always_sea = raster::overlay(range_sea, fun = sum)
always_sea[always_sea != 100*nlayers(range_sea)] = NA
always_sea[always_sea == 100*nlayers(range_sea)] = 1

# get number of temporal intervals for which each cell is suitable for human occupation
range_overlay <- stack()
for(i in rasterlist_range){
  rasterlayer = raster(i)
  rasterlayer[is.na(rasterlayer)] <- 0
  rasterlayer[always_sea == 1] <- NA
  range_overlay <- addLayer(range_overlay, rasterlayer) 
}
range_combined = raster::overlay(range_overlay, fun = sum)
remove(range_sea, range_overlay, always_sea, i)

range_combined_reproj <- projectRaster(range_combined,crs = '+proj=longlat +datum=WGS84 +no_defs') # re-project raster
writeRaster(range_combined_reproj,'range_combined_reproj.tif')

green.colors = colorRampPalette(c("#E4EEE4", "darkgreen"), interpolate = "linear") # function for grey to green color palette
pdf('range_combined.pdf', height = 10, width = 10) # Figure 3
  terra::plot(range_combined, col = c("#F1F1F1", green.colors(20)), axes = F, 
              sub = 'Number of periods through which an area is predicted as suitable for human occupation')
dev.off()

# single figure with all predictions ####

rasterlist_range <- list.files('work/maxent_model_prediction/',full.names=T,pattern='pred_mean')
dates = seq(30,70,2)
j=1
pdf('pred.pdf', height = 21, width = 15) # Supplementary Figure 11
par(mfrow = c(7,3))
par(oma=c(0.25, 0.25, 0.25, 0.25)) # b, l, t, r
par(mar=c(0.25, 0, 0.25, 3.5)) # b, l, t, r
for(i in rasterlist_range) { 
  raster = raster(i)
  terra::plot(raster, axes = F, legend = T)
  legend("bottom", legend = paste(dates[j], "kya    "), bg = "white")
  j=j+1
}
dev.off()

rasterlist_range <- list.files('work/maxent_model_classification/',full.names=T,pattern='range')
dates = seq(30,70,2)
j=1
pdf('range_all.pdf', height = 21, width = 15) # Supplementary Figure 12
par(mfrow = c(7,3))
par(oma=c(0.25, 0.25, 0.25, 0.25)) # b, l, t, r
par(mar=c(0.25, 0, 0.25, 3.5)) # b, l, t, r
for(i in rasterlist_range) { 
  raster = raster(i)
  terra::plot(raster, axes = F, legend = F, col = c("grey80", "forestgreen"))
  legend("bottom", legend = paste(dates[j], "kya    "), bg = "white")
  j=j+1
}
dev.off()
