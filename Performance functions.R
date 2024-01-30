##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################

## This script contains some basic performance functions

## Performance functions
RMSE = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  sqrt(mean((p - o)^2))
}

RMAE = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  (sum(abs(p - o))/length(p)) / mean(o)
}

R2 = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  round(cor(p,o)^2,3)
}
