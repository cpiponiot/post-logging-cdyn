cDyn_inference <- function(var, dG, dM, year, idplot = 1, log_year = NA, prelog_year = NA, 
                           BaG = c(0,Inf), BaM =  c(0,Inf), BbG = c(0,Inf), 
                           BbM = c(0,Inf), Btheta = c(0,Inf), 
                           Bti = c(0,Inf), Btlog = c(0,Inf), 
                           Bvmax =  c(0,Inf), Bmu_vmax = Bvmax, 
                           supSigma = c(Inf, Inf, Inf, Inf, Inf), 
                           n_iter = 2000, n_warmup = 1000, n_chains = 3) {
  
  if (length(var) != length(dG) | length(dG) != length(dM) | length(dM) != length(year))
    stop("Your data (var, dG, dM, year) don't have the same lenght")
  
  dataDyn <- data.table(var, dG, dM, year, idplot)
  
  if ( length(dim(log_year)) == 2 & all(c("idplot","log_year") %in% colnames(log_year))) {
    dataDyn <- merge(dataDyn, log_year, by="idplot", all.x=TRUE)
    
  } else if (length(log_year) %in% c(1,length(year)) ){
    dataDyn$log_year <- log_year
    
  } else stop("log_year is not in the right format")
  
  ### last prelogging year ###
  if ( length(dim(prelog_year)) == 2 & all(c("idplot","prelog_year") %in% colnames(prelog_year))) {
    dataDyn <- merge(dataDyn, prelog_year, by="idplot", all.x=TRUE)
    
  } else if (length(prelog_year) %in% c(1,length(year)) ){
    dataDyn$prelog_year <- prelog_year
    
  } else stop("prelog_year is not in the right format")
  
  ## for the inference, keep only censuses for which t > 0
  subdata <- subset(dataDyn, (!is.na(log_year) & year > log_year) | is.na(log_year))
  
  ## data needed for model inference
  ## np: plot index
  ## cumulative variable changes (for the data.table syntax, see: http://r-datatable.com)
  cumchanges <- subdata[,.(cumG = cumsum(dG),
                           cumM = cumsum(dM), 
                           first_year = min(year),
                           year = year), .(idplot)] 
  subdata <- merge(subdata, cumchanges, by = c("year","idplot"))
  subdata <- subset(subdata, !is.na(log_year) | year > first_year)
  
  setorder(subdata, idplot, year)
  
  ### Model inference using Stan ###
  
  ### data ###
  N <- nrow(subdata)
  np <- as.numeric(as.factor(subdata$idplot)) 
  P <- max(np)
  t <- subdata$year - subdata$log_year
  t[is.na(t)] <- subdata$year[is.na(t)] - subdata$first_year[is.na(t)]
  logged <- unique(np[!is.na(subdata$log_year)])
  L <- length(logged)
  cGrowth <- subdata$cumG
  cMort <- subdata$cumM
  Var <- subdata$var
  maxV <- tapply(dataDyn$var, dataDyn$idplot, max)
  ## deltaV: for all logged plots (treat > 0), the logging-induced loss is 
  ## the pre-logging variable (V[t==min(t)]) minus the agb on the first year after 
  ## logging / silvicultural treatments (V[t==1])
  df <- subset(dataDyn, !is.na(log_year))[order(idplot),.(pre_logging = var[year==min(c(prelog_year, year[year<log_year]), na.rm=TRUE)], 
                                                          post_logging = var[year==min(year[year>log_year])]),.(idplot)]
  deltaV <- df$pre_logging - df$post_logging
  
  Bounds <- data.frame(par = c("aG","aM","bG","bM","theta","ti","tlog","vmax","mu_vmax"),
                       rbind(BaG, BaM, BbG, BbM, Btheta, 
                             Bti, Btlog, Bvmax, Bmu_vmax))
  
  if (any(Bounds[,2:3] < 0) )
    stop("Parameters must be positive.")
  if (any(Bounds[,2] >= Bounds[,3]))
    stop("Lower bounds must be less than upper bounds.")
  
  lowBound <- Bounds[,2]; upBound <- Bounds[,3]
  
  ### inference ###
  stanDynV <- stan("stan models/vdyn_model.stan", iter = n_iter, chains = n_chains, warmup = n_warmup, verbose = FALSE)
  
  return(stanDynV)
  
}