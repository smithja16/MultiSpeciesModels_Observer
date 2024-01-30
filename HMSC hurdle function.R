##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################

## This script contains a function to fit a 'vanilla' Hurdle Hmsc Joint species model
## This is vanilla because it is non-spatial and does not include traits or phylogeny
## This function is not designed around cross-validation, because HMSC does that differently
## For Hmsc we use logNormal for abundance component (Y=log(abundance), family=gaussian)


fit_hmsc_hurdle <- function(train_data_pres,
                            train_data_abund, #use NAs in place of zeros in abundance part
                            species_list,
                            first_species_name, #name of first species, e.g. "Spp1"
                            nChains,  #number of MCMC chains
                            thin,  #MCMC thinning rate
                            samples,  #desired MCMC samples
                            transient,  #MCMC burn in samples
                            nParallel,  #Cores; usually make cores = chains
                            YScale,  #T or F; Scale the response of the abundance component (T recommended) 
                            calc_fit, #T or F (F saves lots of time, but we don't know goodness of fit)
                            save_model,  #T or F
                            save_model_suffix) {
  # total number of iterations = thin * samples + transient
  
  # Save data
  model_fit <- data.frame(species=species_list, AUC_pres=0, R2_abund=0)
  
  # Set up models
  X_pres = train_data_pres[,c("Latitude","Depth_ftm","glorys_sst_C","glorys_mld_m","lunar_illum","Area_swept_km2")]
  X_abund = train_data_abund[,c("Latitude","Depth_ftm","glorys_sst_C","glorys_mld_m","lunar_illum","Area_swept_km2")]
  
  XFormula = ~Latitude + Depth_ftm + poly(glorys_sst_C, degree=2, raw=T) +
    poly(glorys_mld_m, degree=2, raw=T) + poly(lunar_illum, degree=2, raw=T) + Area_swept_km2
  
  XFormula_abund = ~Latitude + Depth_ftm + poly(glorys_sst_C, degree=2, raw=T) +
    poly(lunar_illum, degree=2, raw=T) + Area_swept_km2  #MLD removed from abundance part to avoid overfitting
  
  colx <- which(names(train_data_pres)==first_species_name)  #ID first spp location
  Y_pres <- as.matrix(train_data_pres[,colx:ncol(train_data_pres)])
  Y_abund <- as.matrix(train_data_abund[,colx:ncol(train_data_abund)])
  Y_abund_log <- log(Y_abund)
  
  studyDesign = data.frame(sample=as.factor(1:nrow(X_pres)),
                           vessel=droplevels(train_data_pres$Boat_ID))  #droplevels essential for data subsets
  rL_s <- HmscRandomLevel(units=studyDesign$sample)
  rL_v <- HmscRandomLevel(units=studyDesign$vessel)
  
  Mhm_pres <- Hmsc(Y=Y_pres, XData=X_pres,
                   XFormula=XFormula,
                   studyDesign=studyDesign,
                   ranLevels=list(sample=rL_s, vessel=rL_v),
                   distr="probit")
  
  Mhm_abund <- Hmsc(Y=Y_abund_log, XData=X_abund,
                    XFormula=XFormula_abund,
                    studyDesign=studyDesign,
                    ranLevels=list(sample=rL_s, vessel=rL_v),
                    distr="normal",
                    YScale=YScale)  #YScale only used for normal distribution
  # *^ note the logged response variable
  
  # Fit models
  print("Fitting Mhm_pres model")
  t1 <- Sys.time(); t1
  Mhm_pres <- sampleMcmc(Mhm_pres,
                         thin=thin,
                         samples=samples,
                         transient=transient,
                         nChains=nChains,
                         nParallel=nParallel)
  t2 <- Sys.time()
  print(t2-t1)
  
  print("Fitting Mhm_abund model")
  t1 <- Sys.time(); t1
  Mhm_abund <- sampleMcmc(Mhm_abund,
                          thin=thin,
                          samples=samples,
                          transient=transient,
                          nChains=nChains,
                          nParallel=nParallel)
  t2 <- Sys.time()
  print(t2-t1)
  
  if (save_model == T) {
    saveRDS(Mhm_pres, paste0("Mhm_pres_",save_model_suffix,".rds"))
    saveRDS(Mhm_abund, paste0("Mhm_abund_",save_model_suffix,".rds"))
  }
  
  if (calc_fit == T) {
    # Model fit
    preds_p <- computePredictedValues(Mhm_pres, expected=F)  #expected=F predicts 0/1
    preds_a <- computePredictedValues(Mhm_abund, expected=T)
    
    Mhm_pres_fit <- evaluateModelFit(hM=Mhm_pres, predY=preds_p)
    Mhm_abund_fit <- evaluateModelFit(hM=Mhm_abund, predY=preds_a)

    model_fit$AUC_pres <- Mhm_pres_fit$AUC
    model_fit$R2_abund <- Mhm_abund_fit$R2
    
    # Calculate total biomass
    preds_t <- preds_a
    for (ff in 1:dim(preds_a)[3]) {
      preds_t[,,ff] <- preds_p[,,ff] * exp(preds_a[,,ff])  #undo the log
    }
    preds_t_mean <- apply(preds_t, c(1,2), mean)

    saveRDS(Mhm_pres_fit, paste0("Mhm_pres_fit_",save_model_suffix,".rds"))
    saveRDS(Mhm_abund_fit, paste0("Mhm_abund_fit_",save_model_suffix,".rds"))
  }

  return(list(model_fit=model_fit,
              fitted_totals=preds_t_mean))
  
}

