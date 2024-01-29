##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################


## This code is example code only and will not replicate results as published
## because the observer data we used is confidential. But this code provides
## the model structures and the fundamentals to repeat similar analyses.


###### Load packages ######
library(mgcv)  #GAMs
library(caret) #Performance functions
library(pROC)  #AUC function
library(randomForest)  #RFs
library(smotefamily)  #Smote for RFs
library(Hmsc)  #Joint SDMs
library(abind)  #use with Hmsc


###### Load functions ######
source("Performance functions.R")
source("Cross validation function.R")
source("GAMM hurdle function.R")
source("GAMM tweedie function.R")
source("Random forest hurdle function.R")


###### Load example data ######
# these are mostly dummy data to illustrate formats 
mydata <- readRDS("Example_data.rds")
mytraits <- readRDS("Example_traits.rds")
myphylo <- readRDS("Example_phylo.rds")
#^this is the actual taxonomic tree from study; so won't work with the dummy data


###### Prepare data ######
# Species list
colx <- which(names(mydata)=="Spp1")  #first spp
spp_list <- names(mydata[,colx:ncol(mydata)])

# Create binary and > 0 data for hurdle models
Y <- mydata[,spp_list]
Ypres <- Y
Ypres[Ypres>0] <- 1
mydata_pres <- cbind(mydata[,1:(colx-1)],Ypres)
Yabund <- Y
Yabund[Yabund==0] <- NA
mydata_abund <- cbind(mydata[,1:(colx-1)],Yabund)


###### Fit full models (all data) ######
# GAMMs
gam_hurdle <- fit_gam_stack_hurdle(train_data_pres = mydata_pres,
                                   train_data_abund = mydata_abund,
                                   test_data_pres = NA,
                                   test_data_abund = NA,
                                   species_list = spp_list,
                                   k = 5,
                                   pres_threshold = 0,
                                   pres_prevalence = T,
                                   save_model = T,
                                   save_model_suffix = "runNumber_date")
gam_hurdle$dev_explained  #goodness of fit per taxa
gam_hurdle$predicted_totals  #fitted total biomasses

gam_tweedie <- fit_gam_stack_tweedie(train_data = mydata,
                                     test_data = NA,
                                     species_list = spp_list,
                                     k = 5,
                                     pres_threshold = 0.1,
                                     biomass_threshold = T,
                                     save_model = T,
                                     save_model_suffix = "runNumber_date")
gam_tweedie$dev_explained  #goodness of fit per taxa
gam_tweedie$predicted  #fitted total biomasses

# Random forests
rf_hurdle <- fit_rf_stack_hurdle(train_data_pres = mydata_pres,
                                 train_data_abund = mydata_abund,
                                 test_data_pres = NA,
                                 test_data_abund = NA,
                                 species_list = spp_list,
                                 mtry = 2,  #2 or 3
                                 ntrees = 1200,
                                 pres_threshold = 0,
                                 pres_prevalence = T,
                                 down_smote = F,  #T = class balancing is used; F = not
                                 smote_thresh = 140,  #only works if above is T
                                 save_model = T,
                                 save_model_suffix = "runNumber_date")
rf_hurdle$oob_performance  #out-of-bag (OOB) goodness of fit per taxa

rf_single <- fit_rf_stack(train_data = mydata,
                          test_data = NA,
                          species_list = spp_list,
                          mtry = 2,  #2 or 3
                          ntrees = 1200,
                          pres_threshold = 0.1,  #only used to estimate AUC etc
                          biomass_threshold = T,  #truncates predictions
                          save_model = T,
                          save_model_suffix = "runNumber_date")
rf_single$oob_performance

# Hmsc JSDMs


