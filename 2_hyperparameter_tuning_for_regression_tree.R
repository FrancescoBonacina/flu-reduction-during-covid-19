################################################################################
###
###   Tuning of hyperparameters of the regression tree to predict the flu
###   reduction. 
###   
###   Model: regression tree implemented using the rpart package.
###   Data: observations (countries, seasons) considered in the variable selection.
###   Covariables: the important predictors selected in the variable selection.
###   Hyperparameters to be tuned, in order to avoid overfitting: 
###      - 'cp', the complexity parameter which penalizes the addition of leaves;
###      - 'minbucket', the minimum number of observations which must be in each
###         leaves.
###
###   The hyparameter tuning is implemented through a cross-validation procedure.
###   This scrip provides as result a dataframe with performances of all the
###   models explored through the cross-validation. The optimal parameters are 
###   selected in the next script.
###
################################################################################

rm(list=ls())

### libraries
library(rpart)
library(dplyr)
library(splitstackshape)

### useful functions -------------------------------------------------------------
tree_scoreR2 <- function(tree, yobs){
  # compute the coefficient of determination R² given the tree object of rpart
  # and the vector of the observatios to predict
  yobs_mean <- mean(yobs)
  v <- sum((yobs-yobs_mean)*(yobs-yobs_mean))
  y_residuals <- residuals(tree)
  u <- sum(y_residuals*y_residuals)
  R2 <- 1 - u/v
  return(R2)
} # --> serve ala fine dela parte 3 del codice

scoreR2 <- function(yobs, ypred){
  # compute the coefficient of determination R² from y_observed and y_predicted
  yobs_mean <- mean(yobs)
  v <- sum((yobs-yobs_mean)*(yobs-yobs_mean))
  y_residuals <- yobs-ypred
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
# dataset of observations and covariates
input_source <- paste0(wd, '/inputs/df_observations.txt')
df_obs_ <- read.csv2(input_source, header = T, sep = ',', dec='.')
df_obs_$period <- as.factor(df_obs_$period)   # fix columns 'period' and 'country'
df_obs_$country <- as.factor(df_obs_$country)

# load important predictors
source_selected_predictors <- paste0(wd, sprintf('/outputs/1_selected_variables.txt'))
df_selected_predictors <- read.csv2(source_selected_predictors, header = F)
selected_predictors <- df_selected_predictors$V1

# load list of countries and fluregions.
# Fluregions are macroregions defined by W.H.O which group countries with similar 
# patterns of influenza epidemics. These regions are considered in the following
# in order to perform a stratified cross-validation
source_fluregions <- paste0(wd, '/inputs/FluNet_countries_fluregions.csv')
df_fluregions <- read.csv2(source_fluregions, header = T, sep = ',')

### 2. Re-arrange data
# select only important predictors
variables_to_keep <- c(selected_predictors, c('country', 'period', 'log_rel_flu_level'))
df_obs <- df_obs_[, names(df_obs_) %in% variables_to_keep]

# define groups of countries (based on fluregions) to perform stratified cross-validation
group1 <- c("Temperate South America", "Tropical South America", "Central America and Caribbean")
group2 <- c("North America", "Northern Europe", "South West Europe", "Eastern Europe", "Northern Africa")
group3 <- c("Southern Africa", "Middle Africa", "Eastern Africa", "Western Africa")
group4 <- c("Western Asia", "Central Asia", "Southern Asia", "South-East Asia", "Eastern Asia", "Oceania Melanesia Polynesia")

# add column 'group' to the dataframe of observations
df_obs$group <- 0
for (i in 1:nrow(df_obs)){
  country <- df_obs[i, 'country']
  flureg <- as.character(df_fluregions[df_fluregions$country == country,]$fluregion)
  if (flureg %in% group1){ 
    group <- 1
  } else if (flureg %in% group2){
    group <- 2
  } else if (flureg %in% group3){
    group <- 3
  } else if (flureg %in% group4){
    group <- 4
  }
  df_obs[i,'group'] <- group
}
str(df_obs)

### 3. CROSS-VALIDATION
# Phase space
cpmin <- 0
cpmax <- 0.025
cp_space <- seq(cpmin, cpmax, 0.001)
minbmin <- 4
minbmax <- 18
minbucket_space <- seq(minbmin, minbmax,1)

# parameters for cross-validation
niterations <- 20 #4000
train_proportion <- 0.7
nresearches <- 1  # This can be incremented in order to repeat the procedure several times

# set SEED for reproducibility
set.seed(291121, kind = "L'Ecuyer-CMRG")

# run research
for (nr in 1:nresearches){

  research_results <- list()
  for (m in 1:length(minbucket_space)){
    minb <- minbucket_space[m]
    
    for (c in 1:length(cp_space)){
      cp <- cp_space[c]

      # for this configuration of the parameters fit several trees and compute the 
      # prediction errors
      vec_errors <- numeric(niterations)
      vec_nsplits <- numeric(niterations)
      for (n in 1:niterations){
        # select training and test sets
        data <- tibble::rownames_to_column(df_obs, var = "id") %>% mutate_at(vars(id), as.integer)  # add column with number of the row
        training <- as.data.frame(data %>% stratified(., group = "group", size = train_proportion))  # select training set
        testing <- data[-training$id,]                                                               # select test set
        train_set <- training[, !(names(training) %in% c('id', 'group', 'country', 'period'))] 
        test_set <- testing[, !(names(testing) %in% c('id', 'group', 'country', 'period'))]
      
        # fit
        tree = rpart(data = train_set,
                     formula = log_rel_flu_level ~. ,
                     control = rpart.control(minbucket=minb, cp=cp),
                     xval=0)
        nsplits <- tree_nsplits(tree) 

        # prediction error: err=1-R²=1-ceoff. determination. Estimated on the tests set
        ypred_ <- predict(object = tree, test_set)
        ypred <- as.data.frame(ypred_)$ypred_
        pred_err <- 1 - scoreR2(yobs = test_set$log_rel_flu_level, ypred = ypred)

        # save info
        vec_errors[[n]] <- pred_err
        vec_nsplits[[n]] <- nsplits
      }
    
      # compute mean and sd of the prediction errors
      mean_err <- mean(vec_errors)
      sd_err <- sd(vec_errors)
      
      # compute the average number of splits
      nsplit <- round(mean(vec_nsplits))
      
      # save info
      idx <- length(minbucket_space)*(m-1)+c
      research_results[[idx]] <- data.frame('minbucket'=minb,
                                            'cp'=cp,
                                            'mean_err'=mean_err,
                                            'sd_err'=sd_err,
                                            'nsplits'=nsplit)
    }
    message_minb <- sprintf('*** minb=%d',minb)
    print(message_minb)
  }
  # create dataframe of the results
  df_CV <- do.call(rbind, research_results)
  
  # save dataframe of results
  if (nresearches==1){
    file_name <- sprintf('%s/2_results_CV.csv', outdir_path)
  } else if (nresearches>1){
    file_name <- sprintf('%s/2_results_CV_%d.csv', outdir_path, nr)
  }
  write.table(df_CV, file_name, sep = ',', dec = '.', row.names = F, col.names = T)
}
