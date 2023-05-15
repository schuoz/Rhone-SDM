#Set library
library(biomod2)
library(dplyr)
library(randomForest)
library(caret)
library(reshape2)
library(ecospat)
library(splitTools)
library(FNN)
library(doParallel)

convert_to_myBiomodCV <- function(nom, k = 5) {
  cv_folds <- groupKFold(nom, k)

  # Convert folds to true/false data frame
  tf_df <- data.frame(matrix(NA, nrow = length(nom), ncol = length(cv_folds)))
  colnames(tf_df) <- paste0("fold", seq_along(cv_folds))

  for (i in 1:length(cv_folds)) {
    tf_df[cv_folds[[i]], i] <- TRUE
  }

  tf_df[is.na(tf_df)] <- FALSE

  # Return the resulting matrix
  return(as.matrix(tf_df))
}

wd <- ('/cluster/home/shzong/Revision_sdm/wd_test')
setwd(wd)

var_unbalannce_list <- readRDS("output/var_unbalannce_list.rds")
sp.unbalance <- readRDS("output/full_sp_df.rds")
excluded <- readRDS("output/excluded.rds")
sp.list <- readRDS("output/full_sp_df.rds")
BiomodOptions_unbalance_list <- readRDS("output/BiomodOptions_unbalance_list.rds")
ecospat_data_list <- list()
myCCV_GLM_Models_list <- list()
myCCV_GBM_Models_list <- list()
myCCV_RF_Models_list <- list()
species_list_unbalance <- excluded


# read BIOMOD data from rds data
myBiomodData <- readRDS("output/unbalance_biomod_data_eval.rds")


for (species in species_list_unbalance) {

  cl <- makeCluster(4)
  doParallel::registerDoParallel(cl)

  glm.bm.tuning <- BIOMOD_Tuning(bm.format = myBiomodData[[species]], models = "GLM")
  gbm.bm.tuning <- BIOMOD_Tuning(bm.format = myBiomodData[[species]], models = "GBM")
  rf.bm.tuning <- BIOMOD_Tuning(bm.format = myBiomodData[[species]], models = "RF")

  myBiomodCV <- convert_to_myBiomodCV(sp.list$nom, k=5)
  # Get corresponding presence/absence data
  myResp <- sp.unbalance[, species, drop = FALSE]

  myResp$copy_sp <- sp.unbalance[, species]
  # take the variable names accordingly
  # Get the variable names for the species
  var.unbalance <- var_unbalannce_list[[species]]

  # Get the environmental data
  env_data <- sp.unbalance[, var.unbalance]

  # Get the coordinates
  coord <- sp.unbalance[, c("X_WGS84", "Y_WGS84")]

  #Running all the models for all species
  tryCatch({
    myCCV_GLM_Models <- ecospat.CCV.modeling(sp.data = myResp,
                                             env.data = env_data,
                                             xy = coord,
                                             DataSplitTable = myBiomodCV,
                                             models.esm = "GLM",
                                             minNbPredictors = 2,
                                             VarImport = 6,
                                             ESM = "ALL",
                                             modeling.options.esm = glm.bm.tuning$models.options,
                                             parallel = TRUE)
    myCCV_GBM_Models <- ecospat.CCV.modeling(sp.data = myResp,
                                             env.data = env_data,
                                             xy = coord,
                                             DataSplitTable = myBiomodCV,
                                             models.esm = "GBM",
                                             minNbPredictors = 2,
                                             VarImport = 6,
                                             ESM = "ALL",
                                             modeling.options.esm = gbm.bm.tuning$models.options,
                                             parallel = TRUE)
    myCCV_RF_Models <- ecospat.CCV.modeling(sp.data = myResp,
                                            env.data = env_data,
                                            xy = coord,
                                            DataSplitTable = myBiomodCV,
                                            models.esm = "RF",
                                            minNbPredictors = 2,
                                            VarImport = 6,
                                            ESM = "ALL",
                                            modeling.options.esm = rf.bm.tuning$models.options,
                                            parallel = TRUE)
    myCCV_GLM_Models_list[[species]] <- myCCV_GLM_Models
    myCCV_GBM_Models_list[[species]] <- myCCV_GBM_Models
    myCCV_RF_Models_list[[species]] <- myCCV_RF_Models
  }, error = function(e) {
    message("An error occurred: ", e)
  })

  # Add the results to the ecospat data list
  ecospat_data_list[[species]] <- cbind(myResp, env_data, coord)
}
saveRDS(ecospat_data_list, file = "output/ecospat_data_list.rds")
saveRDS(myCCV_GLM_Models_list, file = "output/myCCV_GLM_Models_list.rds")
saveRDS(myCCV_GBM_Models_list, file = "output/myCCV_GBM_Models_list.rds")
saveRDS(myCCV_RF_Models_list, file = "output/myCCV_RF_Models_list.rds")
