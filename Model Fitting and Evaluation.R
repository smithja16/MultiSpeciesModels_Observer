##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################


## This code is example code only and will not replicate results as published
## because the observer data we used is confidential. But this code provides
## the model structures and the fundamentals to repeat similar analyses.


## Load packages
library(mgcv)  #GAMs
library(pROC)  #AUC function
library(randomForest)  #RFs
library(smotefamily)  #Smote for RFs
library(Hmsc)  #Joint SDMs
library(abind)  #use with Hmsc

## Load functions
source("Performance functions.R")
source("Cross validation function.R")
source("GAMM hurdle function.R")
source("GAMM tweedie function.R")

## Load example data (this is dummy data to illustrate format)
my_data <- readRDS("Example_data.rds")



