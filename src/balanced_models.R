# install libraries
library(biomod2)
library(dplyr)
library(randomForest)
library(caret)
library(reshape2)
library(ecospat)
library(splitTools)
library(FNN)


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

ver <- "all"


# Read "excluded" and "sp.regular" from rds file
excluded <- readRDS("output/excluded.rds")

regular <- readRDS("output/regular.rds")

sp.list <- readRDS("output/full_sp_df.rds")


# Define the grouping variable



species_list_regular <- regular #list of unique species

# Create empty lists to store results
var_importance_list <- list()
var_importance_em_list <- list()
eval_list <- list()
eval_em_list <- list()
eval_null_list <- list()

biomod_out_list <- list()
biomod_em_list <- list()
biomod_nullout_list <- list()

# read BIOMOD data from rds data
myBiomodData <- readRDS("output/balance_biomod_data_eval.rds")

# read BIOMOD null data for model comparison
biomod.null.data <- readRDS(file = "output/balance_null_biomod_data_eval.rds")

# read data split table from RDS file
cv_list_balance <- readRDS(file = "output/cv_list_balance.rds")


for(species in species_list_regular){

  myBiomodCV <- cv_list_balance[[species]]

  bm.tuning <- BIOMOD_Tuning(bm.format = myBiomodData[[species]], models = c("GLM", "GBM", "RF"))
  myBiomodOptions <- bm.tuning$models.options


  # Model single models
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData[[species]],
                                      modeling.id = 'mod.CV',
                                      models = c('RF', 'GBM','GLM'),
                                      bm.options = myBiomodOptions,
                                      nb.rep = 5,
                                      data.split.table = myBiomodCV,
                                      metric.eval = c('TSS','ROC','ACCURACY','KAPPA'),
                                      var.import = 3,
                                      do.full.models = FALSE,
                                      seed.val = 42)

  # Get evaluation scores & variables importance
  myEval <- get_evaluations(myBiomodModelOut)
  myEval$sp <- species

  myBiomodEM=BIOMOD_EnsembleModeling(bm.mod= myBiomodModelOut,
                                     models.chosen="all",
                                     metric.select = "TSS",
                                     metric.select.thresh = 0.3,
                                     em.by="all",
                                     em.algo = c('EMmean', 'EMmedian', 'EMcv', 'EMci', 'EMca', 'EMwmean'),
                                     metric.eval = c('TSS','ROC','ACCURACY','KAPPA'),
                                     EMci.alpha = 0.05,
                                     var.import= 3)

  # Get evaluation scores & variables importance
  myEval_em <- get_evaluations(myBiomodEM)
  myEval_em$sp <- species
  var_imp_em <- get_variables_importance(myBiomodEM)

  null.myBiomodModelOut <- BIOMOD_Modeling(bm.format = biomod.null.data[[species]],
                                           modeling.id = 'mod.CV.null',
                                           models = c( 'RF','GBM','GLM'),
                                           bm.options = myBiomodOptions,
                                           nb.rep = 5,
                                           data.split.table = myBiomodCV,
                                           metric.eval = c('TSS','ROC','ACCURACY','KAPPA'),
                                           var.import = 3,
                                           do.full.models = FALSE,
                                           seed.val = 42)


  # Get evaluation scores & variables importance
  myEval.null <- get_evaluations(null.myBiomodModelOut)
  myEval.null$sp <- species


  myEval$CV.strategy <- "grouped CV"
  var_imp <- get_variables_importance(myBiomodModelOut)

  var_importance_list[[species]] <- var_imp
  var_importance_em_list[[species]] <- var_imp_em
  eval_list[[species]] <- myEval
  eval_em_list[[species]] <- myEval_em
  eval_null_list[[species]] <- myEval.null
  biomod_out_list [[species]] <- myBiomodModelOut
  biomod_nullout_list [[species]] <- null.myBiomodModelOut
  biomod_em_list[[species]] <- myBiomodEM
}

# Combine all variable importance results into a data frame
var_importance_df <- do.call(rbind, var_importance_list)
var_importance_em_df <- do.call(rbind, var_importance_em_list)

# Summarize variable importance across all species
var_importance_summary <- aggregate(var.imp ~ expl.var, data = var_importance_df, FUN = mean)
var_importance_summary <- var_importance_summary[order(var_importance_summary$var.imp, decreasing = TRUE), ]

var_importance_em_summary <- aggregate(var.imp ~ expl.var, data = var_importance_em_df, FUN = mean)
var_importance_em_summary <- var_importance_em_summary[order(var_importance_em_summary$var.imp, decreasing = TRUE), ]

# Combine all evaluation and evaluation difference results into a data frame
eval_df <- do.call(rbind, eval_list)
eval_null_df <- do.call(rbind, eval_null_list)
eval_em_df <- do.call(rbind, eval_em_list)

wd <- ('/cluster/project/ele/shuo')
setwd(wd)
#
save(biomod_em_list, file="biomod_esm.RData")
save(biomod_out_list, file="biomod_out.RData")
write.csv(eval_df, paste0("output/validation_", ver, ".csv"), row.names = TRUE)# change output
write.csv(var_importance_df, paste0("output/var_importance_df", ver, ".csv"), row.names = TRUE)# change output
write.csv(eval_em_df, paste0("output/eval_em_df", ver, ".csv"), row.names = TRUE)# change output
write.csv(var_importance_summary, paste0("output/var_importance_summary", ver, ".csv"), row.names = TRUE)# change output
write.csv(var_importance_em_df, paste0("output/var_importance_em_df", ver, ".csv"), row.names = TRUE)# change output

write.csv(var_importance_em_summary, paste0("output/var_importance_em_summary", ver, ".csv"), row.names = TRUE)# change output
write.csv(eval_null_df, paste0("output/eval_null_df", ver, ".csv"), row.names = TRUE)
