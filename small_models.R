

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


calibrate_model <- function(data, model_type, data_split_table, weighting_score, threshold) {
  my_ESM <- ecospat.ESM.Modeling(data = data, models = model_type, DataSplitTable = data_split_table, weighting.score = weighting_score, parallel = TRUE)
  my_ESM_EF <- ecospat.ESM.EnsembleModeling(my_ESM, weighting.score = weighting_score, threshold = threshold)
  my_ESM_thresholds <- ecospat.ESM.threshold(my_ESM_EF)
  my_Var_contribution <- ecospat.ESM.VarContrib(my_ESM, my_ESM_EF)

  return(list(model = my_ESM, ensemble_model = my_ESM_EF, thresholds = my_ESM_thresholds, var_contribution = my_Var_contribution))
}


run_biomod_model <- function(myBiomodData, nom, species_list, model_types, weighting_score = c("TSS"), threshold = 0, k=5) {
  # Create empty lists to store evaluation, variable contribution, and threshold results
  eval_list <- list()
  var_contribt_list <- list()
  thresholds_list <- list()

  # Loop through the species list and model types
  for (species in species_list) {
    myBiomodCV <- convert_to_myBiomodCV(nom, k)
      # Call the calibrate_model function and store the evaluation, variable contribution, and threshold results in their respective lists
      calibrate_model_output <- calibrate_model(data = myBiomodData[[species]], model_type = model_types, data_split_table = myBiomodCV, weighting_score = weighting_score, threshold = threshold)
      eval_list[[paste(species, model_types, sep="_")]] <- calibrate_model_output$ensemble_model$ESM.evaluations
      var_contribt_list[[paste(species, model_types, sep="_")]] <- setNames(data.frame(t(calibrate_model_output$var_contribution)), rownames(calibrate_model_output$var_contribution))
      thresholds_list[[paste(species, model_types, sep="_")]] <- calibrate_model_output$thresholds
    esm_em_list [[i]] <- calibrate_model_output$ensemble_model
    esm_out_list[[i]] <- calibrate_model_output$my_ESM
  }

  # Combine the evaluation, variable contribution, and threshold results into data frames
  eval_df <- do.call(rbind, eval_list)
  var_contribt_df <- do.call(rbind, var_contribt_list)
  thresholds_df <- do.call(rbind, thresholds_list)

  # Return the data frames as a list
  return(list(eval_df = eval_df, var_contribt_df = var_contribt_df, thresholds_df = thresholds_df, esm_em_list = esm_em_list, esm_out_list = esm_out_list))
}

wd <- ('/cluster/home/shzong/Revision_sdm/wd_test')
setwd(wd)
# Read "excluded" and "sp.regular" from rds file

sp.list <- readRDS("output/full_sp_df.rds")
species_list <- readRDS("output/excluded.rds")
# read BIOMOD data from rds data
myBiomodData <- readRDS("output/unbalance_biomod_data.rds")
myBiomodData_null <- readRDS("output/unbalance_null_biomod_data.rds")
# read data split table from RDS file
cv_list_balance <- readRDS(file = "output/cv_list_unbalance.rds")

esm_out_list <- list()
esm_em_list <- list()
esm_nullout_list <- list()


results <- run_biomod_model(myBiomodData = myBiomodData, nom = sp.list$nom , species_list = species_list, model_types = c("GLM", "GBM", "RF"), weighting_score = c("TSS"), threshold = 0, k=5)

# Write the evaluation, variable contribution, and threshold results to separate CSV files
write.csv(results$eval_df, file = "output/evaluation_results.csv", row.names = FALSE)
write.csv(results$var_contribt_df, file = "output/var_contribt_results.csv", row.names = FALSE)
write.csv(results$thresholds_df, file = "output/thresholds_results.csv", row.names = FALSE)
save(results$esm_em_list, file="esm_em_out.RData")
save(results$esm_out_list, file="esm_out.RData")

results_null <- run_biomod_model(myBiomodData = myBiomodData_null, nom = sp.list$nom , species_list = species_list, model_types = c("GLM", "GBM", "RF"), weighting_score = c("TSS"), threshold = 0, k=5)

# Write the evaluation, variable contribution, and threshold results to separate CSV files
write.csv(results_null$eval_df, file = "output/evaluation_results_null.csv", row.names = FALSE)
write.csv(results_null$var_contribt_df, file = "output/var_contribt_results_null.csv", row.names = FALSE)
write.csv(results_null$thresholds_df, file = "output/thresholds_results_null.csv", row.names = FALSE)
save(results_null$esm_em_list, file="esm_em_null_out.RData")
save(results_null$esm_out_list, file="esm_null_out.RData")

