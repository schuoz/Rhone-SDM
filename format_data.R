sp.list  <- readRDS(file = "output/full_sp_df.rds")
excluded <- readRDS("output/excluded.rds")
species_list_unbalance <- excluded

valid_list <- list()
for (species in species_list_unbalance) {

  # Define the grouping variable
  group <- sp.list[,species]

  # Create the folds
  set.seed(123) # For reproducibility
  folds <- partition(group, p = c(0.8, 0.2), type = "stratified")

  # get the length of each vector in the list
  lengths <- sapply(folds, length)

  # find the maximum length
  max_length <- max(lengths)

  # get the longest vector
  splitIndices <- folds[[which.max(lengths)]]

  valid_list[[species]] <- splitIndices
}

saveRDS(valid_list, "output/valida_list_unbalance.rds")


excluded <- readRDS("output/excluded.rds")
sp.unbalance <- readRDS(file = "output/full_sp_df.rds")
valid_list <- readRDS(file = "output/valida_list_unbalance.rds")
var_unbalannce_list <- readRDS( "output/var_unbalannce_list.rds")
species_list_unbalance <- excluded #list of unique species


# creat a empty list
biomod.data <- list()
for(species in species_list_unbalance){

  # set an example species
  # species <- "Sander.lucioperca"

  #validtion index
  val_ind <- valid_list[[species]]
  # Get corresponding presence/absence data
  myResp <- as.numeric(sp.unbalance[val_ind, species])



  # take the variable names accordingly
  var.unbalance <- var_unbalannce_list[[species]]
  vali_sp <- sp.unbalance[val_ind,]
  test_sp <- sp.unbalance[-val_ind,]
  # Put the data in biomod format
  biomod.data[[species]] <- BIOMOD_FormatingData(resp.var=myResp,
                                                 expl.var=vali_sp[,var.unbalance],
                                                 resp.xy=vali_sp[,c("X_WGS84","Y_WGS84")],
                                                 resp.name=species,
                                                 eval.resp.var = as.numeric(test_sp[, species]),
                                                 eval.expl.var = test_sp[, var.unbalance],
                                                 eval.resp.xy = test_sp[,c("X_WGS84","Y_WGS84")])
}
saveRDS(biomod.data, file = "output/unbalance_biomod_data_eval.rds")

# creat a empty list
biomod.null.data <- list()
for(species in species_list_unbalance){

  val_ind <- valid_list[[species]]
  # Get corresponding presence/absence data
  myResp <- as.numeric(sp.unbalance[val_ind, species])

  # take the variable names accordingly
  var.unbalance <- var_unbalannce_list[[species]]
  null.expl <- sp.unbalance[,var.unbalance]
  # change all the values in the data frame to 1
  null.expl[] <- runif(length(null.expl))
  # Put the data in biomod format
  biomod.null.data [[species]] <- BIOMOD_FormatingData(resp.var=myResp,
                                                       expl.var= null.expl[val_ind,],
                                                       resp.xy=sp.unbalance[val_ind,c("X_WGS84","Y_WGS84")],
                                                       resp.name=species,
                                                       eval.resp.var = as.numeric(sp.unbalance[-val_ind, species]),
                                                       eval.expl.var = null.expl[-val_ind,],
                                                       eval.resp.xy = sp.unbalance[-val_ind,c("X_WGS84","Y_WGS84")])
}
saveRDS(biomod.null.data , file = "output/unbalance_null_biomod_data_eval.rds")