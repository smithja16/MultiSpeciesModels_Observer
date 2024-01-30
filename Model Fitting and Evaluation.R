##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################

## This code is example code only and will not replicate results as published
## because the observer data we used is confidential. But this code provides
## the model structures and the fundamentals to repeat similar analyses.
## This is the 'top level' script that calls other functions to fit models and
## extract/save results.


###### Load packages ######
library(mgcv)  #GAMs
library(caret) #Performance functions
library(pROC)  #AUC function
library(randomForest)  #RFs
library(smotefamily)  #Smote for RFs
library(Hmsc)  #Joint SDMs
library(abind)  #used with Hmsc
library(ape)  #might be needed for phylogeny


###### Load functions ######
source("Performance functions.R")
source("Cross validation function.R")
source("GAMM hurdle function.R")
source("GAMM tweedie function.R")
source("Random forest hurdle function.R")
source("Random forest function.R")
source("HMSC hurdle function.R")
source("HMSC spatial hurdle function.R")


###### Load example data ######
## These are mostly dummy data to illustrate formats 
mydata <- readRDS("Example_data.rds")
mytraits <- readRDS("Example_traits.rds")
myphylo <- readRDS("Example_phylo.rds")
#^this is the actual taxonomic tree from study; it won't work with these dummy data


###### Prepare data ######
## Species list
colx <- which(names(mydata)=="Spp1")  #first spp
spp_list <- names(mydata[,colx:ncol(mydata)])

## Create binary and > 0 data for hurdle models
Y <- mydata[,spp_list]
Ypres <- Y
Ypres[Ypres>0] <- 1
mydata_pres <- cbind(mydata[,1:(colx-1)],Ypres)
Yabund <- Y
Yabund[Yabund==0] <- NA
mydata_abund <- cbind(mydata[,1:(colx-1)],Yabund)


###### Fit full models (all data) ######
## GAMMs
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

## Random forests
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

## Hmsc JSDMs (can take DAYS to fit, so start small and build up)
hmsc_vanilla <- fit_hmsc_hurdle(train_data_pres = mydata_pres,
                                train_data_abund = mydata_abund,
                                species_list = spp_list,
                                first_species_name = "Spp1",
                                nChains = 3,  #final = 3
                                thin = 20,  #final = 20
                                samples = 1000,  #final = 1000
                                transient = 10000,  #final >= 10,000 (evaluate trace plots)
                                nParallel = 3,  #same as nChains
                                YScale = T,
                                calc_fit = T,
                                save_model = T,
                                save_model_suffix = "runNumber_date") 
hmsc_vanilla$model_fit

hmsc_spatial <- fit_hmsc_spatial_hurdle(train_data_pres = mydata_pres,
                                        train_data_abund = mydata_abund,
                                        species_list = spp_list,
                                        first_species_name = "Spp1",
                                        phylo_tree = myphylo,
                                        traits = mytraits,
                                        gpp_nngp = "gpp",
                                        nChains = 3,  #final = 3
                                        thin = 20,  #final = 20
                                        samples = 1000,  #final = 1000
                                        transient = 10000,  #final >= 10,000
                                        nParallel = 3,  #same as nChains
                                        YScale = T,
                                        calc_fit = T,
                                        save_model = T,
                                        save_model_suffix = "runNumber_date")
hmsc_spatial$model_fit


###### Run 5-fold cross-validation on all models ######

## Set up folds
myfolds <- CreateFolds(data=mydata, n_folds=5, n_repeats=3, seed=117)
folds <- myfolds$folds
fold_rows <- myfolds$fold_rows

## CV on GAMM hurdle
save_auc <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_auc), ncol=length(spp_list)))
colnames(dfspp) <- spp_list
save_auc <- cbind(save_auc, dfspp)
save_rmae <- save_auc

for (nn in 1:(n_folds*n_repeats)) {
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  test_dat_presx <- mydata_pres[fold_rows[[nn]],]
  train_dat_presx <- mydata_pres[-fold_rows[[nn]],]
  test_dat_abundx <- mydata_abund[fold_rows[[nn]],]
  train_dat_abundx <- mydata_abund[-fold_rows[[nn]],]
  
  gam_results <- fit_gam_stack_hurdle(train_data_pres = train_dat_presx,
                                      train_data_abund = train_dat_abundx,
                                      test_data_pres = test_dat_presx,
                                      test_data_abund = test_dat_abundx,
                                      species_list = spp_list,
                                      k = 5,
                                      pres_threshold = 0,
                                      pres_prevalence = T,
                                      save_model = F,  #don't save model objects during CV
                                      save_model_suffix = "NA")
  
  save_auc[nn, 3:ncol(save_auc)] <- gam_results$test_r2_rmse$auc_pres
  save_rmae[nn, 3:ncol(save_auc)] <- gam_results$test_r2_rmse$rmae_total
}
boxplot(save_auc[,3:ncol(save_auc)]); abline(h=0.70, col="red")
boxplot(save_rmae[,3:ncol(save_rmae)], ylim=c(0,15)); abline(h=1, col="red")

## CV on GAMM Tweedie
save_auc <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_auc), ncol=length(spp_list)))
colnames(dfspp) <- spp_list
save_auc <- cbind(save_auc, dfspp)
save_rmae <- save_auc

for (nn in 1:(n_folds*n_repeats)) {
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  test_datx <- mydata[fold_rows[[nn]],]
  train_datx <- mydata[-fold_rows[[nn]],]
  
  gam_results <- fit_gam_stack_tweedie(train_data = train_datx,
                                       test_data = test_datx,
                                       species_list = spp_list,
                                       k = 5,
                                       pres_threshold = 0.1,
                                       biomass_threshold = T,
                                       save_model = F,  #don't save model objects during CV
                                       save_model_suffix = "NA")
  
  save_auc[nn, 3:ncol(save_auc)] <- gam_results$test_r2_rmse$auc
  save_rmae[nn, 3:ncol(save_auc)] <- gam_results$test_r2_rmse$rmae
}
boxplot(save_auc[,3:ncol(save_auc)]); abline(h=0.70, col="red")
boxplot(save_rmae[,3:ncol(save_rmae)], ylim=c(0,15)); abline(h=1, col="red")

## CV on Random Forest Hurdle (include class balancing or not)
save_auc <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_auc), ncol=length(spp_list)))
colnames(dfspp) <- spp_list
save_auc <- cbind(save_auc, dfspp)
save_rmae <- save_auc

for (nn in 1:(n_folds*n_repeats)) {
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  test_dat_presx <- mydata_pres[fold_rows[[nn]],]
  train_dat_presx <- mydata_pres[-fold_rows[[nn]],]
  test_dat_abundx <- mydata_abund[fold_rows[[nn]],]
  train_dat_abundx <- mydata_abund[-fold_rows[[nn]],]
  
  rf_results <- fit_rf_stack_hurdle(train_data_pres = train_dat_presx,
                                    train_data_abund = train_dat_abundx,
                                    test_data_pres = test_dat_presx,
                                    test_data_abund = test_dat_abundx,
                                    species_list=spp_list,
                                    mtry = 2,
                                    ntrees = 1200,
                                    pres_threshold = 0,
                                    pres_prevalence = T,
                                    down_smote = T,  #T or F
                                    smote_thresh = 140*0.8,  #only works if above is T; 0.8 bc this is 5-fold CV
                                    save_model = F,  #don't save model objects during CV
                                    save_model_suffix = "NA")
  
  save_auc[nn, 3:ncol(save_auc)] <- rf_results$test_r2_rmse$auc_pres
  save_rmae[nn, 3:ncol(save_auc)] <- rf_results$test_r2_rmse$rmae_total
}
boxplot(save_auc[,3:ncol(save_auc)]); abline(h=0.70, col="red")
boxplot(save_rmae[,3:ncol(save_rmae)], ylim=c(0,15)); abline(h=1, col="red")

## CV on Random Forest
save_auc <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_auc), ncol=length(spp_list)))
colnames(dfspp) <- spp_list
save_auc <- cbind(save_auc, dfspp)
save_rmae <- save_auc

for (nn in 1:(n_folds*n_repeats)) {
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  test_datx <- mydata[fold_rows[[nn]],]
  train_datx <- mydata[-fold_rows[[nn]],]
  
  rf_results <- fit_rf_stack(train_data = train_datx,
                             test_data = test_datx,
                             species_list = spp_list,
                             mtry = 2,
                             ntrees = 1200,
                             pres_threshold = 0.1,
                             biomass_threshold = T,
                             save_model = F,  #don't save model objects during CV
                             save_model_suffix = "NA")
  
  save_auc[nn, 3:ncol(save_auc)] <- rf_results$test_r2_rmse$auc
  save_rmae[nn, 3:ncol(save_auc)] <- rf_results$test_r2_rmse$rmae
}

## CV on Hmsc JSDM vanilla
# This uses the Hmsc appropach to CV due to speed and simplicity
# Load the fitted models saved when using 'fit_hmsc_hurdle'
Mhm_pres <- readRDS()
Mhm_abund <- readRDS()
# For saving data
test_r2_rmse <- data.frame(species=spp_list, auc_pres=NA, r2_abund=NA, r2_total=NA, rmae_total=NA)
# Create partition (* must use same partition for both components)
partitionP = createPartition(hM = Mhm_pres, nfolds = 5, column = "sample")

# Presence component
t1 <- Sys.time()
predsP = computePredictedValues(hM = Mhm_pres,
                                partition = partitionP,
                                nParallel = 3,
                                expected = FALSE)
t2 <- Sys.time(); t2-t1
cv_res <- evaluateModelFit(hM = Mhm_pres, predY = predsP)
test_r2_rmse$auc_pres <- cv_res$AUC 
barplot(cv_res$AUC); abline(h=0.7, col="red")

# Abundance component
t1 <- Sys.time()
predsA = computePredictedValues(hM = Mhm_abund,
                                partition = partitionP,
                                nParallel = 3,
                                expected = TRUE)
t2 <- Sys.time(); t2-t1
cvA_res <- evaluateModelFit(hM = Mhm_abund, predY = predsA)
test_r2_rmse$r2_abund <- cvA_res$R2

# Multiply components for total biomass
preds_t <- predsA
for (ff in 1:dim(predsA)[3]) {
  preds_t[,,ff] <- predsP[,,ff] * exp(predsA[,,ff])  #undo the log
}
preds_t_mean <- apply(preds_t, c(1,2), mean)

# Calculate species-level performance for total biomass
for (ss in 1:length(spp_list)) {
  sppx <- spp_list[ss]
  obs_sppx <- mydata[,sppx]
  pred_sppx <- preds_t_mean[,ss]
  r2_total <- R2(pred_sppx, obs_sppx)
  test_r2_rmse$r2_total[ss] <- r2_total
  rmae_total <- RMAE(pred_sppx, obs_sppx)
  test_r2_rmse$rmae_total[ss] <- rmae_total
}
mean(test_r2_rmse$r2_total)
mean(test_r2_rmse$rmae_total)

## CV on Hmsc JSDM Spatial
# This uses the Hmsc appropach to CV due to speed and simplicity
# Load the fitted models saved when using 'fit_hmsc_spatial_hurdle'
MhmS_pres <- readRDS()
MhmS_abund <- readRDS()
# For saving data
test_r2_rmse <- data.frame(species=spp_list, auc_pres=NA, r2_abund=NA, r2_total=NA, rmae_total=NA)

# Presence component
t1 <- Sys.time()
predsPS = computePredictedValues(hM = MhmS_pres,
                                partition = partitionP,
                                nParallel = 3,
                                expected = FALSE)
t2 <- Sys.time(); t2-t1
cv_res <- evaluateModelFit(hM = MhmS_pres, predY = predsPS)
test_r2_rmse$auc_pres <- cv_res$AUC
barplot(cv_res$AUC); abline(h=0.7, col="red")

# Abundance component
t1 <- Sys.time()
predsAS = computePredictedValues(hM = MhmS_abund,
                                partition = partitionP,
                                nParallel = 3,
                                expected = TRUE)
t2 <- Sys.time(); t2-t1
cvA_res <- evaluateModelFit(hM = MhmS_abund, predY = predsAS)
test_r2_rmse$r2_abund <- cvA_res$R2

# Multiply components for total biomass
preds_t <- predsAS
for (ff in 1:dim(predsAS)[3]) {
  preds_t[,,ff] <- predsPS[,,ff] * exp(predsAS[,,ff])  #undo the log
}
preds_t_mean <- apply(preds_t, c(1,2), mean)

# Calculate species-level performance for total biomass
for (ss in 1:length(spp_list)) {
  sppx <- spp_list[ss]
  obs_sppx <- mydata[,sppx]
  pred_sppx <- preds_t_mean[,ss]
  r2_total <- R2(pred_sppx, obs_sppx)
  test_r2_rmse$r2_total[ss] <- r2_total
  rmae_total <- RMAE(pred_sppx, obs_sppx)
  test_r2_rmse$rmae_total[ss] <- rmae_total
}
mean(test_r2_rmse$r2_total)
mean(test_r2_rmse$rmae_total)

