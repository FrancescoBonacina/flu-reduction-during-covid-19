################################################################################
###
###   Selection of the optimal parameters (cp, minbucket) explored through
###   cross-validation and implemetation of the optimal regression tree to
###   predict flu reduction.
###   
###   1) Selection of the optimal parameters:
###   Following the Breiman's rule we'll select the simplest tree with 
###   prediction error smaller than the minimum prediction error augmented by
###   its standard deviation.
###   As simplest tree we mean the tree which has the smallest average number 
###   of splits, the largest minbucket and the largest cp, respectively.
###
###   2) Implementation of the optimal tree:
###   Model: regression tree implemented using the rpart package.
###   Data: observations (countries, seasons) considered in the variable selection.
###   Covariates: the important predictors selected in the variable selection.
###   Hyperparameters: 'cp' and  'minbucket', selected through cross-validation.
###
################################################################################

rm(list=ls())

### libraries
library(rpart)
library(rpart.plot)
library(dplyr)

# useful functions -------------------------------------------------------------
tree_scoreR2 <- function(tree, yobs){
  # compute the coefficient of determination R² given the tree object of rpart
  # and the vector of the observatios to predict
  yobs_mean <- mean(yobs)
  v <- sum((yobs-yobs_mean)*(yobs-yobs_mean))
  y_residuals <- residuals(tree)
  u <- sum(y_residuals*y_residuals)
  R2 <- 1 - u/v
  return(R2)
}

tree_nsplits <- function(tree){
  # compute the number of splits in the tree
  df_cptable <- as.data.frame(tree$cptable)
  nsplits <- df_cptable$nsplit[nrow(df_cptable)]
  return(nsplits)
}
# ------------------------------------------------------------------------------

### set directory for output files
wd <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir <- sprintf('outputs/')
dir.create(file.path(wd, outdir), showWarnings = FALSE)
outdir_path <- paste0(wd, '/', outdir)

### 1. Load input data
# dataset of observations
source_obs <- paste0(wd, sprintf('/inputs/df_observations.txt'))
df_obs_ <- read.csv2(source_obs, header = T, sep = ',', dec='.')
df_obs_$period <- as.factor(df_obs_$period)
df_obs_$country <- as.factor(df_obs_$country)

# load important predictors
source_selected_predictors <- paste0(wd, sprintf('/outputs/1_selected_variables.txt'))
df_selected_predictors <- read.csv2(source_selected_predictors, header = F)
selected_predictors <- df_selected_predictors$V1

# load results of cross-validation for hyperparameter tuning
source_cv <- paste0(wd, sprintf('/outputs/2_results_CV.csv'))
df_cv <- read.csv2(source_cv, header=T, sep = ',', dec='.')
str(df_cv)  

### 2. Re-arrange data
# select only important predictors
variables_to_keep <- c(selected_predictors, c('country', 'period', 'log_rel_flu_level'))
df_obs <- df_obs_[, names(df_obs_) %in% variables_to_keep]

### 3. Select optimal hyperparameters
# sort trees according to their prediction errors
df_CVsorted <- df_cv[order(df_cv$mean_err),]

# compute the Breiman's threshold: minimum error + its standard deviation
threshold <- df_CVsorted[1,'mean_err'] + df_CVsorted[1,'sd_err']

# select all models that are 'equivalently' optimal
df_optimal_trees <- df_CVsorted[(df_CVsorted$mean_err <= threshold),]   

# sort those trees from the simplest one to the most complex.
df_optimal_trees <- df_optimal_trees[with(df_optimal_trees, order(nsplits, -minbucket, -cp)),]

# select the hyperparameters of the optimal tree (the Breiman's tree)
breiman_minb <- df_optimal_trees[1,'minbucket']
breiman_cp <- df_optimal_trees[1,'cp']
mean_pred_err <- df_optimal_trees[1,'mean_err']

### 4. Implement the optimal tree
# run the model
df_tree <- df_obs[, !(names(df_obs) %in% c('country', 'period'))]
str(df_tree)
tree_breiman = rpart(data = df_tree,
                     formula = log_rel_flu_level ~. ,
                     control = rpart.control(cp=breiman_cp, 
                                             minbucket=breiman_minb),
                     xval=0)

# compute the performance score and the number of splits and plot the tree
scoreR2_breiman <- tree_scoreR2(tree_breiman, df_tree$log_rel_flu_level)
nb_splits <- tree_nsplits(tree_breiman)
title <- sprintf('Tree: cp=%4.3f, minbucket=%d, R^2=%4.3f, nb_splits=%d', breiman_cp, breiman_minb, scoreR2_breiman, nb_splits)
rpart.plot(tree_breiman, main = title, extra = 1, digits = 3)

### 5. Save results
df_pred1 <- df_obs[rownames(df_tree), c('country', 'period')]
df_pred2 <- rpart.predict(tree_breiman, df_tree, nn=TRUE)  # prediction for the observations
df_predictions <- cbind(df_pred1, df_pred2)
colnames(df_predictions) <- c('country', 'period', 'leaf_log_rel_flu_level', 'leaf_label')

out_file <- sprintf('%s/3_tree_predictions.csv', outdir_path)
write.table(df_predictions, out_file, sep = ',', dec = '.', row.names = F, col.names = T)
