########################################################
# Helpers for inference scripts
########################################################

infer_networks <- function(cooc_tab, lag, opt, thr, use_default_offset,
                           hourcuts, hour_offset ="none",
                           offset, merge_off_by, ncols = 2,
                           cov, merge_cov_by, 
                           formula){
  # This function infers a network given a cooc table and inference parameters.
  # cooc_tab: the coocurrence table frow which to infer the nk.
  # lag: the time lag to use (defaults to 1).
  # opt: character string, type of aggregation.
  # thr: filtering threshold.
  # use_offset: whether to use PLN automated offset (boolean).
  # offset: a custom offset to be used (eg for clusters). Named, names corresponding    #     to...
  # merge_off_by: string, culumn name in the df the offset corresponds to.
  # covariates: a custom covariates vector to be used (eg for habitat). Named, names    #     corresponding to...
  # formula is the formula to use in PLN (optional; if ommited, deduced
  if(hour_offset == "duration"){
    error("Cannot handle duration offset yet")
  }
  
  if(missing(lag)){
    lag <- ddays(1) # Anyways wont be used
  }
  processed <- process_data(cooc_tab, lag = lag, 
                            opt = opt, threshold = thr, hourcuts = hourcuts)
  
  # Prepare data using offset and covar
  if(use_default_offset | missing(offset) | hour_offset == "TSS"){ # No offset, or default
    if(missing(cov)){
      f.prep <- prepare_mydata(processed, ncols = ncols)
    }else{
      f.prep <- prepare_mydata(processed, ncols = ncols,
                               cov = cov, merge_cov_by = merge_cov_by)
    }
  }else{ # Use custom offset
    if(missing(cov)){
      f.prep <- prepare_mydata(processed, ncols = ncols, offset = offset,
                               merge_offset_by = merge_off_by)
    }else{
      f.prep <- prepare_mydata(processed, ncols = ncols, 
                               offset = offset, merge_offset_by = merge_off_by,
                               cov = cov, merge_cov_by = merge_cov_by)
    }
  }

  f.nks <- PLNnetwork(formula, f.prep$data)
  
  return(list("network" = f.nks,
              "prep" = f.prep))
}

pick_network <- function(nks, f.prep, crit){
  # Gets the best model according to selecion criterion.
  # nks is a PLN object containing inferred nks.
  # f.prep is the prepared df used to infer the nks.
  # crit: the inference criterion (StARS or BIC.)
  
  if(nrow(f.prep$data$Abundance) < 10 & crit == "StARS"){
    message("StARS cannot be used, using BIC instead.")
    f.nk <- getBestModel(nks, crit = "BIC")
  }
  else{
    f.nk <- getBestModel(nks, crit = crit)
  }
  return(f.nk)
}

active_duration <- function(df, span){
  # Compute active camera duration given a df with columns site ID and DateTimeOriginal.
  # span is a lubridate dd... object
  # Returns a df with total duration, active duration, their difference, 
  #     the min and max date and the proposed offset.
  
  # Group the observations by date (at 10 days)
  active_days <- df %>% dplyr::select(site_ID, DateTimeOriginal) %>%
    group_by(site_ID) %>% 
    mutate(dategroup = as.character(cut(DateTimeOriginal, 
                                        breaks = seq(min(DateTimeOriginal), max(DateTimeOriginal) + span, 
                                                     by = span))))
  # Cout active duration
  count_days <- active_days %>% group_by(site_ID) %>%
    summarise(active_duration = length(unique(dategroup))*span, # active duration defined as n timespans * timespan ( so +/- timespan)
              min_date = min(DateTimeOriginal),
              max_date = max(DateTimeOriginal)) %>% 
    mutate(total_duration = as.duration(max_date - min_date)) %>%
    mutate(difference = total_duration - active_duration) %>% # time difference between measured and 'real' timespan
    mutate(offset = ifelse(abs(difference) < span, as.numeric(total_duration, 'days'), 
                           as.numeric(active_duration, 'days'))) # offset is defined as the tital duration if difference is not too big, else approximated by active duration
 
  count_days <- count_days %>% dplyr::select(site_ID, offset, total_duration, active_duration, everything())
  
  return(count_days) 
}