# README

This repository contains a demo version of the code used in the analyses presented in the paper *Global patterns and drivers of influenza decline during the COVID-19 pandemic*, by Bonacina, F et al.

### Overview of the analysis

This study analyzes the global influenza reduction during the first year and a half of the Covid-19 pandemic (March 2020 to September 2021). Influenza samples by country are aggregated by trimesters, and after removing country-trimesters not satisfying the inclusion critera, 93 countries were considered, up to 6 trimesters each, for a total of 330 data points (country-trimesters). 

Potential predictors of the influenza reduction were identified by analyzing 20 covariates, including country factors (geographic, meteorological, demographic, and health preparedness factors) and variables associated with COVID-19 pandemic. 

In a first step, the most important predictors of influenza reduction are identified. Then, influenza reduction is predicted for each country-trimester by using a Classification and Regression Tree that includes covariates identified as significant. 

### Content

- The `inputs` folder contains two input files:
    - `df_observations.txt`: log relative influenza level and 19 covariates for 330 countries and trimesters.  
      Specifically the following covariates are included:
        - `log_rel_flu_level`
        - `cov_daily_cases`
        - `workplace_presence_reduction`
        - `idvi`
        - `age`
        - `T`
        - `RH`
        - `latitude`
        - `longitude`
        - `nb_days_school_closure`
        - `nb_days_workplace_closure`
        - `nb_days_public_event_restrictions`
        - `nb_days_gathering_restrictions`
        - `nb_days_public_transport_restrictions`
        - `nb_days_stay_at_home_requirements`
        - `nb_days_facial_covering_requirements`
        - `nb_days_elderly_shielding`
        - `nb_days_international_travel_restrictions`
        - `nb_days_testing_implementation`
        - `nb_days_contact_tracing_implementation`

      Air-travel data is not included since they are commercially available and restrictions apply to redistribution of this data. Source of data are listed in the document [data_sources.md](https://github.com/FrancescoBonacina/flu-reduction-during-covid-19/blob/main/data_sources.md).
    
    - `FluNet_countries_fluregions.csv`: list of countries included in the FluNet repository and influenza region they belong to. Influenza regions are macro-regions defined by W.H.O which group countries with similar influenza epidemics. This information is used to implement stratified cross-validation in script `2_hyperparameters_tuning_for_regression_tree.R`.

- The `outputs` folder contains the outputs returned by the R scripts.

- The main steps of the analyses are implemented in three R scripts:
  1) In `1_selection_of_predictors_of_flu_reduction.R` the *VSURF* algorithm is used to estimate the importance of each predictors and select the most significant ones.

      *Output*: `1_selected_variables.txt` with the list of the important predictors.     
      *Computation time*: ~10 min.

  2) `2_hyperparameters_tuning_for_regression_tree.R`. The regression tree is implemented by using the *rpart* package. To control overfitting we play with tune two hyperparameters, *minbucket* and *cp*, which regularize the structure of the tree. A cross-validation procedure is implemented in order to explore numerous combinations of the two parameters. 

      *Output*: `2_results_CV.csv`, values tested for *minbucket* and *cp* are reported here. In addition, estimates of the mean and standard deviation of the prediction error and the mean number of tree subdivisions are given for each combination of hyperparameters.    
      *Computation time*: ~ 6-8 hours.

  3) `3_optimal_tree.R`. First, the optimal combination of the previously explored parameters is identified following the Breiman’s rule. The optimal tree is the simplest tree with prediction error less than the minimum prediction error plus its standard deviation. The simplest tree is determined by looking at the smallest number of splits on average, the largest *minbucket*, and the largest *cp*, in order. Second, the *rpart* package is used to implement the optimal tree.

      *Output*: `3_tree_predictions.csv`, list of the predicted values of flu reduction and the leaf classification of each country-trimester.   
      *Computation time*: < 1 min.

### Software requirements

This code is supported for *Linux*, *macOS* and *Windows*. The scripts has been tested on the following systems:
- Linux: Ubuntu 18.04
- macOS: Monterey 12.2

The code is implemented in the `R` language (version 3.6.1) and uses the following packages:
- VSURF (1.1.0)
- rpart (4.1.16)
- rpart.plot (3.1.0)
- dplyr (1.0.9)
- splitstackshape (1.4.8)

The packages are publicly available for installation from the [CRAN](https://cran.r-project.org/) web servers.

### License

This project is covered under the **MIT License**.
