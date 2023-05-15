
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

average_esm_prediction <- function( sp.unbalance, exclude, species, predictors, i) {
  # Function to perform normal prediction using ecospat model
  # sp.unbalance: data frame containing the full data for build sdm
  # exclude: vector containing the rows are not using in the traning dataset
  # species: species to build sdm
  # predictors: vector containing predictors of SDMs
  # i: intger indicating the number of row
  myResp <- as.numeric(sp.balance[-exclude, species])

  # Format the data for biomod
  biomod.data <- BIOMOD_FormatingData(resp.var = myResp,
                                      expl.var = sp.unbalance[-exclude, predictors],
                                      resp.xy = sp.unbalance[-exclude,c("X_WGS84","Y_WGS84")],
                                      resp.name = species)


  # Prepare data for prediction
  myExpl <- sp.balance[i,predictors]


  # Model ecospat
  my_ESM <- ecospat.ESM.Modeling(data = biomod.data, models = c("GLM", "GBM", "RF"), parallel = TRUE)

  my_ESM_EF <- ecospat.ESM.EnsembleModeling(my_ESM, weighting.score = c("TSS"), threshold = 0.4)

  my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my_ESM,
                                              new.env=myExpl)

  myESMProj <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
                                                 ESM.EnsembleModeling.output=my_ESM_EF)

  # Get the predictions
  prediction <- get_predictions(myESMProj)$pred

  return(prediction)
}




binomial_cv_group <- function(var.unbalance, data, group_var, species) {
  # Function to perform leave-one-group-out and jackknife cross-validation
  # var.unbalance: predictors of SDMs
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
    predictions$prediction[i] <- average_esm_prediction (sp.unbalance = data, exclude = exclude, species = species, predictors = var.unbalance, i = i)
  }

  # Return the predictions and observations data frame
  return(predictions)
}




binomial_cv_group_slurmarray <- function(var.unbalance, data, group_var, species) {
  # Function to perform leave-one-group-out and jackknife cross-validation using slurm arrary
  # var.unbalance: predictors of SDMs
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
    preds <- average_esm_prediction(sp.unbalance = data, exclude = exclude, species = species, predictors = var.unbalance, i = i)
    # Append the predicted value to the predictions data frame
    predictions <- cbind(predictions[i,], preds)
  }, error = function(e) {
    message("An error occurred: ", e)
  })


  wd <- ('/cluster/project/ele/shuo/output')

  setwd(wd)
  saveRDS(predictions, paste0("predicton_sdm_un_", species, "_" ,i, ".rds"))

  # Return the predictions and observations data frame
  return(predictions)
}




binomial_cv_group_foreach <- function(var.unbalance, data, group_var, species) {
  # Function to perform leave-one-group-out and jackknife cross-validation using paralell
  # var.unbalance: predictors of SDMs
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
    predictions$prediction[i] <- average_esm_prediction (sp.unbalance = data, exclude = exclude, species = species, predictors = var.unbalance, i = i)
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
excluded <- readRDS("output/excluded.rds")
var_unbalannce_list <- readRDS( "output/var_unbalannce_list_divide10.rds")
species_list_unbalance <- excluded

# initialize a list for save predictions
predicton_list <- list()

for (species in species_list_unbalance) {
  # Define the cluster for parallel processing

  # take the variable names accordingly
  var.unbalance <- var_unbalannce_list[[species]]

  predicton_list[[species]] <- binomial_cv_group_slurmarray(var.balance = var.unbalance, data = sp_df, group_var = "nom", species = species)

}


