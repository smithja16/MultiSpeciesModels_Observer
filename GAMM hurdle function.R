##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################

## This script contains a function to fit a Hurdle (delta) GAMM
## It has function arguments to allow for a full model fit (test data is "NA") or
## allow for cross validation fits (test data is specified)


fit_gam_stack_hurdle <- function(train_data_pres,
                                 train_data_abund, #both sets of train data must have equal rows (use NAs in abund part)
                                 test_data_pres,  #can be "NA" to ignore testing
                                 test_data_abund,  #can be "NA" to ignore testing
                                 species_list,
                                 k,               #wiggliness
                                 pres_threshold,  #keep = zero; used only if below = F
                                 pres_prevalence,  #T or F (use T if 'species richness' is important)
                                 save_model,  #T or F; don't save models during a CV process
                                 save_model_suffix) {  #only if above is T
  # * note that deviance explained can include random effects
  # * use k=5 to avoid more coefficients than data for rarest species
  # * for abundance (>0) component, some rows will now have no observations and give warnings
  
  ## Create files to save data
  save_fit_pres <- as.data.frame(matrix(data=NA, nrow=nrow(train_data_pres), ncol=length(species_list)))
  save_fit_abund <- save_fit_pres
  save_fit_total <- save_fit_pres
  colnames(save_fit_pres) <- species_list
  colnames(save_fit_abund) <- species_list
  colnames(save_fit_total) <- species_list
  dev_expl <- data.frame(species=species_list, dev_expl_pres=0, dev_expl_abund=0)
  test_r2_rmse <- data.frame(species=species_list, auc_pres=0, accuracy_pres=0, specificity_pres=0, f1_score_pres=0, 
                             r2_abund=0, r2_total=0, rmse_abund=0, rmse_total=0, rmae_total=0)
  if (is.na(test_data_pres)[1] == T) { 
    save_test_pres <- 0
    save_test_abund <- 0
    save_test_total <- 0
  } else {
    save_test_pres <- save_fit_pres[1:nrow(test_data_pres),]
    save_test_abund <- save_test_pres
    save_test_total <- save_test_pres
  }
  
  pb <- txtProgressBar(min=0, max=length(species_list))
  
  for (ss in 1:length(species_list)) {
    
    setTxtProgressBar(pb,ss)
    sppx <- species_list[ss]
    
    Mg_pres <- gam(get(sppx) ~ s(Latitude, k=k) +
                     s(glorys_sst_C, k=k) +
                     s(glorys_mld_m, k=k) +
                     s(Depth_ftm, k=k) +
                     s(lunar_illum, k=3) +
                     s(Area_swept_km2, k=k) + 
                     s(Boat_ID, bs="re"),
                   data=train_data_pres,
                   family=binomial)

    Mg_abund <- gam(get(sppx) ~ s(Latitude, k=k) +
                      s(glorys_sst_C, k=k) +
                      s(lunar_illum, k=3) +
                      s(Depth_ftm, k=k) +
                      s(Area_swept_km2, k=3),
                    data=train_data_abund,  
                    family=Gamma(link="log"))  #or default "inverse" works similarly
    #covariates and some 'k' reduced here to avoid too many coefficients

    if (save_model == T) {  #Save the models for prediction (maps)
      saveRDS(Mg_pres, paste0("Spp_",ss,"Mg_pres_",save_model_suffix,".rds"))
      saveRDS(Mg_abund, paste0("Spp_",ss,"Mg_abund_",save_model_suffix,".rds"))
    }
    
    if (pres_prevalence == T) {  #calculate presence threshold for this species based on its prevalence
      colxx <- which(names(train_data_pres) == sppx)
      pres_thresholdx <- train_data_pres[,colxx]
      pres_threshold <- length(pres_thresholdx[pres_thresholdx==1])/nrow(train_data_pres)
      if (pres_threshold > 0.95) { pres_threshold <- 0.95 }  #max is 0.95 - added to avoid issues with very prev. species
    }
    
    #Calculate fitted values
    #predict can give warning about missing factor levels (for REs), but they are ignored in the prediction
    P_pres_l <- predict(Mg_pres, train_data_pres, type="link", exclude="s(Boat_ID)")
    P_abund_l <- predict(Mg_abund, train_data_abund, type="link")
    P_pres <- Mg_pres$family$linkinv(P_pres_l)
    P_pres[P_pres < pres_threshold] <- 0
    P_abund <- Mg_abund$family$linkinv(P_abund_l)
    P_total <- P_pres*P_abund
    
    #Calculate predictions on test data
    if (is.na(test_data_pres)[1] == T) { 
    } else {
      
      P_pres_test_l <- predict(Mg_pres,  test_data_pres, type="link", exclude="s(Boat_ID)")
      P_abund_test_l <- predict(Mg_abund, test_data_abund, type="link")
      P_pres_test <- Mg_pres$family$linkinv(P_pres_test_l)
      P_abund_test <- Mg_abund$family$linkinv(P_abund_test_l)
      P_pres_test[P_pres_test < pres_threshold] <- 0
      P_total_test <- P_pres_test*P_abund_test
      
      save_test_pres[,ss] <- P_pres_test
      save_test_abund[,ss] <- P_abund_test
      save_test_total[,ss] <- P_total_test
      
      test_data_abund[,sppx][is.na(test_data_abund[,sppx])] <- 0  #put zeros back in where NAs currently are
      obs_total_test <- test_data_pres[,sppx]*test_data_abund[,sppx]
      
      r2_abund <- R2(P_abund_test, test_data_abund[,sppx])
      r2_total <- R2(P_total_test, obs_total_test)
      
      P_pres_test_pa <- P_pres_test
      P_pres_test_pa[P_pres_test_pa > 0] <- 1
      
      if (sum(test_data_pres[,sppx]) == 0 | mean(test_data_pres[,sppx]) == 1) {
        auc_pres <- NA  #AUC can't be calculated if test data are all zeros or all ones
        accuracy <- NA
        specificity <- NA
        f1_score <- NA
      } else {
        auc_pres <- suppressMessages(auc(test_data_pres[,sppx],P_pres_test))  #1 is best; measures performance of ranking TNs and TPs by their predicted probabilities (TNs lower prob than TPs)
        conf_mat <- caret::confusionMatrix(data=as.factor(P_pres_test_pa),
                                           reference=as.factor(test_data_pres[,sppx]),
                                           positive="1",
                                           mode="everything") #all the confusion stuff is subject to the 'pres_threshold' so probs less reliable than AUC
        accuracy <- conf_mat$overall["Accuracy"] #1 is best, biased by lots of zeros
        specificity <- conf_mat$byClass["Specificity"]  #how good are we at occurrences (1s)
        f1_score <- conf_mat$byClass["F1"]  #1 is best, balances precision and recall - how good are we at absences (0s)
      }
      rmse_abund <- RMSE(P_abund_test, test_data_abund[,sppx])  #predicted then observed
      rmae_total <- RMAE(P_total_test, obs_total_test)
      rmse_total <- RMSE(P_total_test, obs_total_test)
      
      #Save out-of-sample performance
      test_r2_rmse$r2_abund[ss] <- r2_abund
      test_r2_rmse$r2_total[ss] <- r2_total
      test_r2_rmse$auc_pres[ss] <- auc_pres
      test_r2_rmse$accuracy_pres[ss] <- accuracy
      test_r2_rmse$specificity_pres[ss] <- specificity
      test_r2_rmse$f1_score_pres[ss] <- f1_score
      test_r2_rmse$rmse_abund[ss] <- rmse_abund
      test_r2_rmse$rmse_total[ss] <- rmse_total
      test_r2_rmse$rmae_total[ss] <- rmae_total
    }
    
    #Save predictions and in-sample performance
    save_fit_pres[,ss] <- P_pres
    save_fit_abund[,ss] <- P_abund
    save_fit_total[,ss] <- P_total
    dev_expl$dev_expl_pres[ss] <- round(summary(Mg_pres)$dev.expl,3)
    dev_expl$dev_expl_abund[ss] <- round(summary(Mg_abund)$dev.expl,3)
    
  }
  
  close(pb)
  
  return(list(predicted_presences=save_fit_pres,
              predicted_abundances=save_fit_abund,
              predicted_totals=save_fit_total,
              dev_explained=dev_expl,
              test_predicted_presences=save_test_pres,
              test_predicted_abundances=save_test_abund,
              test_predicted_totals=save_test_total,
              test_r2_rmse=test_r2_rmse))
}

