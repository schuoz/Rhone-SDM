# Rhone-SDM
This repository has species distribution models and other essential codes for running the species distribution models (SDM) in Rhone. The scripts were written by Shuo Zong (Shzong AT ethz dot ch), please contact him if you have any questions.
## Format
The code was written in R mainly using species distribution model packages ([BIOMOD](https://biomodhub.github.io/biomod2/index.html) and [Ecospat](https://www.unil.ch/ecospat/en/home.html))

### Step 1 : Data preparation - [format_data.R]
The data need to be formated in a biomod data format to input to BIOMOD and Ecospat models. This step includes split the data into traning and testing dataset, choose the varibles lists for SDM for each species.

### Step 2 : Model building - [balanced_models.R, small_models.R and small_models_importance.R]
Choose the model according to the occurence of the species along the river (balanced spread species and unbalanced spread species), Evaluate the SDM and Run variable importance for each SDM.

### Step 3 : Uncertainty analysis - [uncertainty_small_models.R, and uncertanty_balance_models.R]
Use jack-knife to evaluate the uncertainties, predict the occurrences for each location for all the species using the three modelling techniques (GLMs, GBMs, and RFs), without considering the sampling location.

### Step 4 : Prediction and Mapping -[balance_prediction.R and small_models_prediction.R]
Map SDM for both balanced spread and unblanced spread species using qualified model


