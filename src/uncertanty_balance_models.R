
# install libraries
library(biomod2)
library(dplyr)
library(randomForest)
library(caret)
library(reshape2)
library(ecospat)
library(splitTools)
library(FNN)
library(foreach)
library(doParallel)
library(parallel)

average_biomod_prediction <- function( sp.balance, exclude, species, predictors, i) {
  # Function to perform normal prediction using biomod
  # sp.balance: data frame containing the full data for build sdm
  # exclude: vector containing the rows are not using in the traning dataset
  # species: species to build sdm
  # predictors: vector containing predictors of SDMs
  # i: intger indicating the number of row
  myResp <- as.numeric(sp.balance[-exclude, species])

  # Format the data for biomod
  biomod.data <- BIOMOD_FormatingData(resp.var = myResp,
                                      expl.var = sp.balance[-exclude, predictors],
                                      resp.xy = sp.balance[-exclude,c("X_WGS84","Y_WGS84")],
                                      resp.name = species)

  bm.tuning <- BIOMOD_Tuning(bm.format = biomod.data, models = c("GLM", "GBM", "RF"))
  myBiomodOptions <- bm.tuning$models.options

  # Model biomod
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = biomod.data,
                                      models = c('RF', 'GBM','GLM'),
                                      bm.options = myBiomodOptions,
                                      nb.rep = 2,
                                      metric.eval = c('TSS','ROC','ACCURACY','KAPPA'),
                                      var.import = 3,
                                      do.full.models = FALSE,
                                      seed.val = 42)

  # Prepare data for prediction
  myExpl <- sp.balance[i,predictors]

  # Project using biomod
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    proj.name = 'Current',
                                    new.env = myExpl,
                                    models.chosen = 'all')

  # Get the predictions
  prediction <- get_predictions(myBiomodProj)$pred

  return(prediction)
}



binomial_cv_group <- function(var.balance, data, group_var, species) {
  # Function to perform leave-one-group-out and jackknife cross-validation
  # var.balance: predictors of SDMs
  # data: data frame containing the full data for build sdm
  # group_var: column name of the group variable
  # species: species to build sdm

  # Initialize a data frame to store the predictions and observations
  predictions <- data.frame(observation = data[[species]], prediction = rep(NA, nrow(data)))

  # Loop through each row in the data and fit the model on the rest of the data
  for (i in 1:nrow(data)) {
    # Exclude the i-th row and its group from the training set
    group_i <- data[[group_var]][i]
    exclude <- which(data[[group_var]] == group_i)
    predictions$prediction[i] <- average_biomod_prediction (sp.balance = data, exclude = exclude, species = species, predictors = var.balance, i = i)
  }

  # Return the predictions and observations data frame
  return(predictions)
}





binomial_cv_group_slurmarray <- function(var.balance, data, group_var, species) {
  # Function to perform leave-one-group-out and jackknife cross-validation using slurm arrary
  # var.balance: predictors of SDMs
  # data: data frame containing the full data for build sdm
  # group_var: column name of the group variable
  # species: species to build sdm

  # Initialize a data frame to store the predictions and observations
  i <- as.numeric(args[1])
  # Exclude the i-th row and its group from the training set

  temp_dir <- tempdir()

  setwd(temp_dir)
  # Save the temporary files in the directory
  saveRDS(i, paste0("text", i, ".rds"))

  predictions <- data.frame(observation = data[[species]], prediction = rep(NA, nrow(data)))

  # Loop through each row in the data and fit the model on the rest of the data

  # Exclude the i-th row and its group from the training set
  group_i <- data[[group_var]][i]
  exclude <- which(data[[group_var]] == group_i)
  tryCatch({
    predictions$ID <- args[1]
    predictions$sp <- species
    preds <- average_biomod_prediction(sp.balance = data, exclude = exclude, species = species, predictors = var.balance, i = i)
    # Append the predicted value to the predictions data frame
    predictions <- cbind(predictions[i,], preds)
  }, error = function(e) {
    message("An error occurred: ", e)
  })


  wd <- ('/cluster/project/ele/shuo')

  setwd(wd)
  saveRDS(predictions, paste0("predicton_sdm_", species, "_" ,i, ".rds"))

  # Return the predictions and observations data frame
  return(predictions)
}




binomial_cv_group_foreach <- function(var.balance, data, group_var, species) {
  # Function to perform leave-one-group-out and jackknife cross-validation using paralell
  # var.balance: predictors of SDMs
  # data: data frame containing the full data for build sdm
  # group_var: column name of the group variable
  # species: species to build sdm


  # Initialize a data frame to store the predictions and observations
  predictions <- data.frame(observation = data[[species]], prediction = rep(NA, nrow(data)))

  # Loop through each row in the data and fit the model on the rest of the data
  predictions <- foreach(i=1:nrow(data), .combine=rbind) %dopar% {
    # Exclude the i-th row and its group from the training set
    group_i <- data[[group_var]][i]
    exclude <- which(data[[group_var]] == group_i)
    predictions$prediction[i] <- average_biomod_prediction (sp.balance = data, exclude = exclude, species = species, predictors = var.balance, i = i)
    # Return the predictions data frame for the i-th row
    return(predictions)
  }



  # Return the predictions data frame
  return(predictions)
}






args = commandArgs(trailingOnly = TRUE);

# set wd
wd <- ('/cluster/home/shzong/Revision_sdm/wd_test')
setwd(wd)



# read input from wd
sp_df <- readRDS("output/full_sp_df.rds")
regular <- readRDS("output/regular.rds")
var_balannce_list <- readRDS("output/var_balannce_list.rds")
biomod_option <- readRDS("output/BiomodOptions_list.rds")
species_list_regular <- regular

# initialize a list for save predictions
predicton_list <- list()

for (species in species_list_regular) {
  # Define the cluster for parallel processing

  # take the variable names accordingly
  var.balance <- var_balannce_list[[species]]

  predicton_list[[species]] <- binomial_cv_group_slurmarray(var.balance = var.balance, data = sp_df, group_var = "nom", species = species)
}


