##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################

## This script contains a function to fit a Tweedie GAMM
## It has function arguments to allow for a full model fit (test data is "NA") or
## allow for cross validation fits (test data is specified)


fit_gam_stack_tweedie <- function(train_data,
                                  test_data,
                                  species_list,
                                  k,              #wiggliness
                                  pres_threshold,  #only used for AUC etc approximations; 0.1 works well
                                  biomass_threshold,  #use a spp-specific minimum biomass threshold (below that is a zero)
                                  save_model,
                                  save_model_suffix) {
  #pres_threshold is only to generate zeros for AUC etc estimates, which aren't very reliable here
  
  ## Create files to save data
  save_fit <- as.data.frame(matrix(data=NA, nrow=nrow(train_data), ncol=length(species_list)))
  colnames(save_fit) <- species_list
  dev_expl <- data.frame(species=species_list, dev_expl=0)
  test_r2_rmse <- data.frame(species=species_list, auc=0, accuracy=0, specificity=0, f1_score=0,
                             r2=0, rmse=0, rmae=0)
  if (is.na(test_data)[1] == T) { 
    save_test <- 0
  } else {
    save_test <- save_fit[1:nrow(test_data),]
  }
  
  pb <- txtProgressBar(min=0, max=length(species_list))
  
  for (ss in 1:length(species_list)) {
    
    setTxtProgressBar(pb,ss)
    sppx <- species_list[ss]
    
    Mg_tw <- gam(get(sppx) ~ s(Latitude, k=k) +
                   s(glorys_sst_C, k=k) +
                   s(glorys_mld_m, k=k) +
                   s(Depth_ftm, k=k) +
                   s(lunar_illum, k=3) +
                   s(Area_swept_km2, k=k) +
                   s(Boat_ID, bs="re"),
                 data=train_data,
                 family=tw())
    #k for lunar-Illum reduced to allow only basic domed relationships and avoid overfitting

    if (save_model == T) {
      saveRDS(Mg_tw, paste0("Spp_",ss,"Mg_tw_",save_model_suffix,".rds"))
    }
    
    #Calculate fitted values
    #predict can give warning about missing factor levels (for REs), but they are ignored in any case
    P <- predict(Mg_tw, train_data, type="response", exclude="s(Boat_ID)")
    
    if (biomass_threshold == T) {  #calculate threshold for this species based on its prevalence
      colxx <- which(names(train_data) == sppx)
      biom_thresholdx <- train_data[,colxx]
      biom_thresholdx <- biom_thresholdx[biom_thresholdx>0]
      biom_threshold <- min(biom_thresholdx)*0.5  #threshold is half the minimum, to avoid overly penalizing small values
      
      P[P <= biom_threshold] <- 0  #truncate fitted values
    }
    
    #Calculate predictions on test data
    if (is.na(test_data)[1] == T) { 
    } else {
      
      P_test <- predict(Mg_tw,  test_data, type="response", exclude="s(Boat_ID)")
      if (biomass_threshold == T) {
        P_test[P_test <= biom_threshold] <- 0  #truncate test values
      }
      save_test[,ss] <- P_test
      
      P_test_pa <- P_test
      P_test_pa[P_test_pa < pres_threshold] <- 0  #this is only used to calc. AUC etc
      P_test_pa[P_test_pa > 0] <- 1  #convert to pres/abs
      test_data_pa <- test_data[,sppx]
      test_data_pa[test_data_pa > 0] <- 1  #convert to pres/abs
      
      if (sum(test_data_pa) == 0 | mean(test_data_pa) == 1) {
        Auc <- NA  #AUC can't be calculated if test data are all zeros or all ones
        accuracy <- NA
        specificity <- NA
        f1_score <- NA
      } else {
        Auc <- suppressMessages(auc(test_data_pa, P_test_pa))
        conf_mat <- confusionMatrix(data=as.factor(P_test_pa),
                                    reference=as.factor(test_data_pa),
                                    positive="1",
                                    mode="everything")
        accuracy <- conf_mat$overall["Accuracy"] #1 is best, biased by lots of zeros
        specificity <- conf_mat$byClass["Specificity"]  #how good are we at zeros (1 is best)
        f1_score <- conf_mat$byClass["F1"]  #1 is best, balances precision and recall - how good are we at absences (0s)
      }
      
      r2 <- R2(P_test, test_data[,sppx])
      rmse <- RMSE(P_test, test_data[,sppx])
      rmae <- RMAE(P_test, test_data[,sppx])
      
      #Save out-of-sample performance
      test_r2_rmse$auc[ss] <- Auc
      test_r2_rmse$accuracy[ss] <- accuracy
      test_r2_rmse$specificity[ss] <- specificity
      test_r2_rmse$f1_score[ss] <- f1_score
      test_r2_rmse$r2[ss] <- r2
      test_r2_rmse$rmse[ss] <- rmse
      test_r2_rmse$rmae[ss] <- rmae
      
    }
    
    #Save predictions and in-sample performance
    save_fit[,ss] <- P
    dev_expl$dev_expl[ss] <- round(summary(Mg_tw)$dev.expl,3)
    
  }
  
  close(pb)
  
  return(list(predicted=save_fit,
              dev_explained=dev_expl,
              test_predicted=save_test,
              test_r2_rmse=test_r2_rmse))
  
}
