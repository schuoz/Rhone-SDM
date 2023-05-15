library(biomod2)
library(dplyr)
library(randomForest)
library(caret)
library(raster)

wd <- ('/cluster/project/ele/shuo')

setwd(wd)
#--------------------
# prepare a raster stack
#--------------------
# Get raster file from path
chl_a_med_layer <- paste0("rhone_RW_chl_a_med_v3.tif")
chl_a_std_layer <- paste0("rhone_RW_chl_a_std_v4.tif")

hum_max_layer <- paste0("rhone_R_hum_modi_max_v3.tif")
hum_med_layer <- paste0("rhone_R_hum_modi_v3.tif")
hum_std_layer <- paste0("rhone_R_hum_modi_std_v4.tif")

ep_std_layer <- paste0("rhone_R_ep_std_v3.tif")
ep_med_layer <- paste0("rhone_R_ep_v3.tif")

evi_med_layer <- paste0("rhone_R_EVI_v3.tif")
evi_std_layer <- paste0("rhone_R_EVI_std_v3.tif")

gpp_med_layer <- paste0("rhone_R_GPP_v3.tif")
gpp_std_layer <- paste0("rhone_R_GPP_std_v4.tif")

ndvi_med_layer <- paste0("rhone_R_NDVI_v3.tif")
ndvi_std_layer <- paste0("rhone_R_NDVI_std_v3.tif")

lst_med_layer <- paste0("rhone_RW_lst_med_v3.tif")
lst_std_layer <- paste0("rhone_RW_lst_std_v3.tif")

sd_med_layer <- paste0("rhone_RW_sd_med_v3.tif")
sd_std_layer <- paste0("rhone_RW_sd_s2_std_v4.tif")

dem_layer <- paste0("rhone_R_dem_v3.tif")
dem_std_layer <- paste0("rhone_R_dem_std_v4.tif")

slope_layer <- paste0("rhone_R_slope_v3.tif")
slope_layer <- paste0("rhone_R_slope_v3.tif")

tsi_med_layer <- paste0("rhone_RW_tsi_s2_med_v4.tif")
tsi_std_layer <- paste0("rhone_RW_tsi_s2_std_v4.tif")

water_sum_layer <- paste0("rhone_RW_water_sum_v3.tif")

# read file as raster
chl_a_med <- raster(chl_a_med_layer)
chl_a_std <- raster(chl_a_std_layer)
hum_max <- raster(hum_max_layer)
hum_med <- raster(hum_med_layer)
hum_std <- raster(hum_std_layer)
evi_med <- raster(evi_med_layer)
evi_std <- raster(evi_std_layer)
ndvi_med <- raster(ndvi_med_layer)
ndvi_std <- raster(ndvi_std_layer)
ep_med <- raster(ep_med_layer)
ep_std <- raster(ep_std_layer)
gpp_med <- raster(gpp_med_layer)
gpp_std <- raster(gpp_std_layer)
lst_med <- raster(lst_med_layer)
lst_std <- raster(lst_std_layer)
sd_med <- raster(sd_med_layer)
sd_std <- raster(sd_std_layer)
tsi_med <- raster(tsi_med_layer)
tsi_std <- raster(tsi_std_layer)
dem <- raster(dem_layer)
dem_std <- raster(dem_std_layer)
slope <- raster(slope_layer)
water_sum <- raster(water_sum_layer)

# stack them together for prediction (predictors)
dem_re <- projectRaster(dem,sd_med,method = 'bilinear')
#dem_re <- projectRaster(dem,sd_med,method = 'bilinear')
slope_re <- projectRaster(slope,sd_med,method = 'bilinear')



# Read "excluded" and "sp.regular" from rds file
excluded <- readRDS("output/excluded.rds")

regular <- readRDS("output/regular.rds")

sp.list <- readRDS("output/full_sp_df.rds")

load(file = "output/biomod_esm.RData")
load(file = "output/biomod_out.RData")
# Define the grouping variable



species_list_regular <- regular #list of unique species
regular <- readRDS("output/regular.rds")
var_balannce_list <- readRDS("output/var_balannce_list.rds")

#------------------

#---------------------
# Extract the name of the variables
#---------------------
var <- c('tsi_l8_med','tsi_s2_med','rst_med','sd_s2_med',
         'sd_l8_med','chl_a_med','evi_med','ndvi_med','ep_med','gpp_med',
         'hum_modi_med','slope_med',
         'tsi_l8_std','tsi_s2_std','rst_std','sd_s2_std',
         'sd_l8_std','chl_a_std','evi_std','ndvi_std','ep_std','gpp_std',
         'hum_modi_std','ele_std','slope_std','water_sum')
var.rmv <- c("ndvi_med","tsi_s2_std","tsi_l8_std","sd_l8_std","sd_l8_med","tsi_l8_med"
  ,"ele_std", "chl_a_med","slope_std",'ele_med')
var <- setdiff(var, var.rmv)
#---------------------

#-----------------------
# stack them together for prediction (predictors)
#-----------------------
preds <- stack(tsi_med,lst_med,sd_med,evi_med,ep_med,gpp_med,hum_med,slope_re,lst_std,sd_std,chl_a_std,evi_std,
               ndvi_std,ep_std,gpp_std,hum_std,water_sum)
names(preds) <- var
#-------------


# Define empty list to store biomod data
current.projections=list()

for(i in species_list_regular){

  # take the variable names accordingly
  var <- var_balannce_list[[i]]

  current.projections[[i]]=BIOMOD_EnsembleForecasting(EM.output=biomod.ensembles[[i]],
                                                      new.env=preds[[var]],
                                                      proj.name="current",
                                                      binary.meth="TSS",
                                                      output.format=".img",
                                                      do.stack=FALSE,
                                                      build.clamping.mask = F)
  # save map
  pdf(paste0(i,"/",i,"distribution_map.pdf"))
  pred_current=get_predictions(current.projections[[i]])
  plot(pred_current)
  dev.off()
}

# Save list of current projections
save(current.projections,file="current_projections.RData")
#--------------

#-----------------------
#plot
#-----------------------
for(i in species_list_regular){

  # save map
  pdf(paste0(i,"/",i,"distribution_map_with_train_data_PA.pdf"))
  pred_current <- raster(current.projections[[i]]@proj@link[1])
  plot(pred_current,main = paste0(i,"distribution_map_with_train_data.pdf"))
  xy_species <- cbind(biomod.data[[i]]@coord,species = biomod.data[[i]]@data.species)
  points(x = xy_species$lon[xy_species$species == 1],
         y = xy_species$lat[xy_species$species == 1],
         pch = 6, cex = 0.5)
  points(x = xy_species$lon[xy_species$species == 0],
         y = xy_species$lat[xy_species$species == 0],
         pch = 0, cex = 0.5)
  dev.off()
}
#--------------

print("done")