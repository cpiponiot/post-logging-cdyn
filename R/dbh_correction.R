# # # ## test
# dat = subset(data, i_arbre == "145494")[order(campagne)]
# X = dat$circonf/pi
# tm = dat$campagne
# limit = 10
# ladder = dat$ladder
# max_allowed_growth=4
# min_allowed_decrease = 0

dbh_correction <- function(X,
                           tm,
                           limit = 10,
                           ladder = NULL, 
                           max_allowed_growth = 5, 
                           min_allowed_decrease = -2) {
  # Xsav if for browser() use: save initial value of X
  Xsav <- X
  limit = min(limit)
  # first_census=min(c(first_census,tm))
  
  # cresc_abs: absolute annual diameter increment
  cresc <- rep(0, length(X) - 1)
  cresc_abs <- rep(0, length(X) - 1)
  if (sum(!is.na(X)) > 1) {
    cresc[which(!is.na(X))[-1] - 1] <-
      diff(X[!is.na(X)]) / diff(tm[!is.na(X)])
    cresc_abs[which(!is.na(X))[-1] - 1] <- diff(X[!is.na(X)])
  }
  
  if (length(cresc) > 0) {
    
    ### First, correct when there is a raised measurement
    if (sum(ladder,na.rm=TRUE) > 0) {
      raised = which(diff(c(NA, ladder)) == 1)
      if (length(raised) == 0 |
          isTRUE(X[raised] > X[raised - 1]) ) {
        raised = which.min(cresc) + 1
      }
      for (rs in raised)
        X[rs:length(X)] = X[rs:length(X)] + (X[rs - 1] - X[rs]) + mean(cresc[-(rs -1)])
      
      ## update cresc
      cresc <- rep(0, length(X) - 1)
      cresc_abs <- rep(0, length(X) - 1)
      if (sum(!is.na(X)) > 1) {
        cresc[which(!is.na(X))[-1] - 1] <-
          diff(X[!is.na(X)]) / diff(tm[!is.na(X)])
        cresc_abs[which(!is.na(X))[-1] - 1] <- diff(X[!is.na(X)])
      }
    }
    
    ####    if there is a DBH change > 5cm/year or < -2 cm   ####
    ### do as many corrections as there are abnormal DBH change values ###
    Ncresc_abn = sum(abs(cresc) >= max_allowed_growth | cresc_abs < min_allowed_decrease)
    
    if (Ncresc_abn > 0) {
      for (i in 1:Ncresc_abn) {
        
        # begin with the census with the highest DBH change
        cresc_abn = which(abs(cresc) >= max_allowed_growth | cresc_abs < min_allowed_decrease)
        ab <- cresc_abn[which.max(abs(cresc[cresc_abn]))]
        
        if (length(ab) == 1) {
          # values surrounding ab
          surround = c(ab - 2, ab - 1, ab + 1, ab + 2)
          # that have a meaning (no NAs or 0 values)
          surround = surround[surround > 0 &
                                surround <= length(cresc)]
          
          # mean DBH change around ab
          meancresc = max(mean(cresc[surround], na.rm = TRUE), 0)
          
          # moment of max and min DBH changes around ab (including ab, that should be one of the 2)
          sourround_ab = sort(c(surround, ab))
          up = sourround_ab[which.max(cresc[sourround_ab])]
          down = sourround_ab[which.min(cresc[sourround_ab])]
          
          if (length(surround) > 0) {
            # 1st case : excessive increase/decrease offset by a similar decrease in dbh, plus 5cm/yr
            # is there a value that could compensate the excessive DBH change?
            # check if removing those values would solve the problem (ie cresc < 5 & cresc_abs > -2 )
            if (isTRUE(down > up & cresc[up] * cresc[down] < 0 &
                       # first an increase and then a decrease in DBH
                       (X[down + 1] - X[up]) / (tm[down + 1] - tm[up])  <= max_allowed_growth &
                       X[down + 1] - X[up] >= min_allowed_decrease) |
                isTRUE(up > down & cresc[up] * cresc[down] < 0 &
                       # first an decrease and then a increase in DBH
                       (X[up + 1] - X[down]) / (tm[up + 1] - tm[down])  <= max_allowed_growth &
                       X[up + 1] - X[down] >= -min_allowed_decrease)) {
              # correction: abnormal values are deleted and will be replaced later on (see missing)
              first <- min(up, down) + 1
              last <- max(up, down)
              X[first:last] <- NA
            }
            
            
            # 2nd case: abnormal DBH change with no return to initial values
            # we trust the set of measurements with more values
            # if they are the same size, then we trust the last one
            # ladders?
            else {
              if ((sum(!is.na(X[1:ab])) > sum(!is.na(X))/2)) {
                X[(ab + 1):length(X)] <- X[(ab + 1):length(X)] - (X[ab + 1] - X[ab]) + mean(cresc[-ab]) * diff(tm)[ab]
              } else {
                X[1:ab] <-
                  X[1:ab] + (X[ab+1]-X[ab]) - meancresc * diff(tm)[ab]
              } 
            }
          }
          
          if (length(X[!is.na(X)]) == 2 & i==1){
            ## only two values, with abnormal difference
            # trust the second one
            X[!is.na(X)][1] <- X[!is.na(X)][2] 
          }
          
          # cresc_abs: absolute annual diameter increment
          cresc <- rep(0, length(X) - 1)
          cresc_abs <- rep(0, length(X) - 1)
          if (sum(!is.na(X)) > 1) {
            cresc[which(!is.na(X))[-1] - 1] <-
              diff(X[!is.na(X)]) / diff(tm[!is.na(X)])
            cresc_abs[which(!is.na(X))[-1] - 1] <- diff(X[!is.na(X)])
          }
        }
      }
    }
    
    ## replace missing values
    if (any(!is.na(X))) { X <- repl_missing(X, tm) } else { X = rep(0, length(X)) }
    
  }
  
  return(X)
}