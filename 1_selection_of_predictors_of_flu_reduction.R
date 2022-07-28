################################################################################
###
###   Perform the selection of predictors of flu reduction through 
###   models based on Random Forests
###
###   Ref: Robin Genuer, Jean-Michel Poggi, Christine Tuleau-Malot. VSURF: An R 
###        Package for Variable Selection Using Random Forests. The R Journal, 
###        R Foundation for Statistical Computing, 2015, 7(2), pp.19-33. 
###        hal-01251924
###
################################################################################

rm(list=ls())

### libraries
library(VSURF)

### set directory for output files
wd <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
outdir <- sprintf('outputs/')
dir.create(file.path(wd, outdir), showWarnings = FALSE)
outdir_path <- paste0(wd, '/', outdir)

### 1. load input files: data and set of variables to consider   ###############
input_source <- paste0(wd, '/inputs/df_observations.txt')
df_obs <- read.csv2(input_source, header = T, sep = ',', dec='.')
str(df_obs)
colnames(df_obs)
### 2. Variable selection
# set SEED for reproducibility
set.seed(121221, kind = "L'Ecuyer-CMRG")
  
# data for regression analysis: keep only columns of covariates
df_data <- df_obs[, !names(df_obs) %in% c('country', 'period', 'iso_code')]
str(df_data)

# nb of observations and nb of variables
nobs <- nrow(df_data)
nvariables <- ncol(df_data)   # all variables, including the target variable
npredictors <- nvariables-1   # excluding the target variable
message1 <- sprintf('***   Number of predictors = %d, nb of observations = %d', npredictors, nobs)
print (message1)
  
# perform future selection
flu.VSURF <- VSURF(formula = log_rel_flu_level ~ .,
                   data = df_data,
                   ntree = 8000,       # default 2000
                   #mtry = 7,          # default nb of predictors/3,
                   nfor.thres = 100,   # default 50
                   nfor.interp = 100,  # default 25
                   nfor.pred = 100,    # default 25
                   na.action = na.omit,  # to deal with missing values: na.omit = drop the rows with missing values
                   parallel = TRUE,      #                              na.roughfix = impute missing values starting from median values and then the proximity matrix from the randomForest is used to update the imputation
                   clusterType='FORK',
                   ncores = 8)

# plot variable importance 
plot(flu.VSURF, step = "thres", imp.sd = FALSE, var.names = TRUE)

# keep important variables and export them 
interpretation_var <- flu.VSURF$varselect.interp
selected_variables <- c()
for (i in interpretation_var){
  selected_variables <- c(selected_variables, colnames(df_data)[i+1])
}
selected_variables
df_selected <- data.frame(selected_variables)
dfname <- sprintf('%s/1_selected_variables.txt', outdir_path)  # sensitivity about nb processed per week 
write.table(df_selected, dfname, sep = ',', dec = '.', row.names = F, col.names = F)
