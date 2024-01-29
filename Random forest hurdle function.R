##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################

## This script contains a function to fit a Hurdle Random Forest
## It has function arguments to allow for a full model fit (test data is "NA") or
## allow for cross validation fits (test data is specified)


fit_rf_stack_hurdle <- function(train_data_pres,
                                train_data_abund,
                                test_data_pres,  #can be "NA" to ignore CV testing
                                test_data_abund,
                                species_list,
                                mtry,  #number of variables to try at each split
                                ntrees,  #number of trees
                                pres_threshold,  #leave = zero, only used when below = F
                                pres_prevalence,  #T or F
                                down_smote,  #if T uses downsampling and SMOTE
                                smote_thresh,  #number of occurrences below which SMOTE is used instead of downsampling
                                save_model,  #T or F
                                save_model_suffix) {  
  
  ## Create files to save data
  save_fit_pres <- as.data.frame(matrix(data=NA, nrow=nrow(train_data_pres), ncol=length(species_list)))
  save_fit_abund <- save_fit_pres
  save_fit_total <- save_fit_pres
  colnames(save_fit_pres) <- species_list
  colnames(save_fit_abund) <- species_list
  colnames(save_fit_total) <- species_list
  oob_perf <- data.frame(species=species_list, oob_error_rate_pres=NA, oob_r2_abund=NA)
  test_r2_rmse <- data.frame(species=species_list, auc_pres=NA, accuracy_pres=NA, specificity_pres=NA, 
                             f1_score_pres=NA, r2_abund=NA, r2_total=NA,
                             rmse_abund=NA, rmse_total=NA, rmae_total=NA)
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
    
    if (length(unique(train_data_pres[,species_list[ss]])) == 1) {  #if training data is all 1s or all 0s, can't do presence RF
      ss <- ss + 1 #skip this species for this data subset
    }
    
    setTxtProgressBar(pb,ss)
    
    sppx <- species_list[ss]
    train_data_pres[,sppx] <- as.factor(train_data_pres[,sppx])
    
    if (down_smote == F) {
      
      # Fit presence RF without class balancing
      Mrf_pres <- randomForest(as.formula(paste0(sppx," ~ Latitude + Depth_ftm + glorys_sst_C +
                                         glorys_mld_m + lunar_illum + Area_swept_km2")),
                               data=train_data_pres, importance=T,
                               keep.inbag=T, num.trees=ntrees, mtry=mtry,
                               keep.forest=T)
    }
    
    if (down_smote == T) {
      
      num_positive <- nrow(train_data_pres[train_data_pres[,sppx] == 1,])
      num_absent <- nrow(train_data_pres[train_data_pres[,sppx] == 0,])
      
      if (num_positive > smote_thresh & num_absent > smote_thresh) {  #if both classes (0s, 1s) are more than 'smote_thresh' use downsampling
        
        train_down_pres <- train_data_pres
        train_down_pres$Class_down <- train_down_pres[,sppx]
        train_down_pres <- caret::downSample(x=train_data_pres[,-ncol(train_down_pres)],  #remove 'Class'
                                             y=as.factor(train_down_pres$Class_down))
        
        # Fit presence RF with class balancing - downsampling
        Mrf_pres <- randomForest(as.formula(paste0(sppx," ~ Latitude + Depth_ftm + glorys_sst_C +
                                         glorys_mld_m + lunar_illum + Area_swept_km2")),
                                 data=train_down_pres, importance=T,
                                 keep.inbag=T, num.trees=ntrees, mtry=mtry,
                                 keep.forest=T)
        
      } else {  #if one class is < smote_thresh use SMOTE
        
        colsx <- which(names(train_data_pres) %in% c("Latitude","Depth_ftm","glorys_sst_C",
                                                     "glorys_mld_m","lunar_illum","Area_swept_km2"))
        
        if (num_positive > 5 & num_absent > 5) { 
          K <- 5 
        } else {
          K <- min(c(num_positive, num_absent)) - 1  #can't have more neighbours than data points
        }
        train_smote_pres <- smotefamily::SMOTE(X=train_data_pres[,colsx],
                                               target=as.factor(train_data_pres[,sppx]),
                                               K=K)  #number of neighbours to get data from
        train_smote_pres2 <- train_smote_pres$data
        names(train_smote_pres2)[names(train_smote_pres2) == 'class'] <- sppx
        train_smote_pres2[,sppx] <- as.factor(train_smote_pres2[,sppx])
        
        # Fit presence RF with class balancing - SMOTE
        Mrf_pres <- randomForest(as.formula(paste0(sppx," ~ Latitude + Depth_ftm + glorys_sst_C +
                                         glorys_mld_m + lunar_illum + Area_swept_km2")),
                                 data=train_smote_pres2, importance=T,
                                 keep.inbag=T, num.trees=ntrees, mtry=mtry,
                                 keep.forest=T)
      }
      
    }
    
    train_data_abund2 <- train_data_abund  #remove all rows that are non-zero for RESPONSE SPECIES ONLY (won't work with NAs)
    rows0 <- which(is.na(train_data_abund2[,sppx]))  #which rows are zero for this species
    if (length(rows0) > 0) {
      train_data_abund2 <- train_data_abund2[-rows0,]
    }
    
    # Fit abundance RF
    Mrf_abund <- randomForest(as.formula(paste0(sppx," ~ Latitude + Depth_ftm + glorys_sst_C +
                                         glorys_mld_m + lunar_illum + Area_swept_km2")),
                              data=train_data_abund2, importance=T,
                              keep.inbag=T, num.trees=ntrees, mtry=mtry,
                              keep.forest=T)

    if (save_model == T) {
      saveRDS(Mrf_pres, paste0("Spp_",ss,"Mrf_pres_",save_model_suffix,".rds"))
      saveRDS(Mrf_abund, paste0("Spp_",ss,"Mrf_abund_",save_model_suffix,".rds"))
    }
    
    if (pres_prevalence == T) {  #calculate threshold for this species based on its prevalence
      colxx <- which(names(train_data_pres) == sppx)
      pres_thresholdx <- train_data_pres[,colxx]
      pres_threshold <- length(pres_thresholdx[pres_thresholdx==1])/nrow(train_data_pres)
      if (pres_threshold > 0.95) { pres_threshold <- 0.95 }
      save_prev$prevalence[ss] <- pres_threshold
    }
    
    #calculate fitted values
    P_pres <- predict(Mrf_pres, newdata=train_data_pres, type="prob")[,2]  #probability of presence
    P_pres[P_pres < pres_threshold] <- 0
    P_abund <- predict(Mrf_abund, newdata=train_data_abund, type="response")
    P_total <- P_pres * P_abund
    
    #Calculate predictions on test data
    if (is.na(test_data_pres)[1] == T) { 
    } else {
      
      P_pres_test <- predict(Mrf_pres, test_data_pres, type="prob")[,2]
      P_pres_test[P_pres_test < pres_threshold] <- 0
      P_abund_test <- predict(Mrf_abund, test_data_abund, type="response")
      P_total_test <- P_pres_test*P_abund_test
      
      save_test_pres[,ss] <- P_pres_test
      save_test_abund[,ss] <- P_abund_test
      save_test_total[,ss] <- P_total_test
      
      test_data_abund[,sppx][is.na(test_data_abund[,sppx])] <- 0
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
        auc_pres <- suppressMessages(auc(test_data_pres[,sppx],P_pres_test))
        conf_mat <- confusionMatrix(data=as.factor(P_pres_test_pa),
                                    reference=as.factor(test_data_pres[,sppx]),
                                    positive="1",
                                    mode="everything")
        accuracy <- conf_mat$overall["Accuracy"] #1 is best, biased by lots of zeros
        specificity <- conf_mat$byClass["Specificity"]  #how good are we at occurrences (1s)
        f1_score <- conf_mat$byClass["F1"]
      }
      rmse_abund <- RMSE(P_abund_test, test_data_abund[,sppx])
      rmse_total <- RMSE(P_total_test, obs_total_test)
      rmae_total <- RMAE(P_total_test, obs_total_test)
      
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
    
    #Save predictions and model performance
    if (down_smote == F) {  #otherwise there are diff numbers of fitted values
      save_fit_pres[,ss] <- P_pres
      save_fit_abund[,ss] <- P_abund
      save_fit_total[,ss] <- P_total
    }
    oob_perf$oob_error_rate_pres[ss] <- Mrf_pres$err.rate[,1][length(Mrf_pres$err.rate[,1])]  #fraction of misclassified samples (smaller is better)
    oob_perf$oob_r2_abund[ss] <- round(Mrf_abund$rsq[length(Mrf_abund$rsq)],3)  #their OOB R2
    
  }
  
  close(pb)
  return(list(predicted_presences=save_fit_pres,
              predicted_abundances=save_fit_abund,
              predicted_totals=save_fit_total,
              oob_performance=oob_perf,
              test_predicted_presences=save_test_pres,
              test_predicted_abundances=save_test_abund,
              test_predicted_totals=save_test_total,
              test_r2_rmse=test_r2_rmse))
  
}

