########################################################
# Functions used to import and aggregate data
########################################################

library(tidyverse)
library(StreamMetabolism)
library(ggplot2)
library(lubridate)
library(sf)
library(chron)

# Define global variable to eliminate non species
NOT_SPP <- c("fire", "UNIDENTIFIED")
BYCATCH <- c("bat","birdofprey", "birdother", "bustardkori", "bustardludwigs",
             "rodent", "squirreltree", "reptilesamphibians", 
             "secretarybird","tortoise", "domesticanimal", "craneblue")


read_files <- function(patterns, path = "../indep_records_spp_30mins"){
  # path: the directory to look into (default: indep_records_spp_30mins)
  # patterns: a vector of patterns to match in the filename 
  # This function assumes that there are at least two files matching the given
  # (list of) pattern.
  regex <- paste0(path,"/*.csv")
  files = Sys.glob(regex)
  files = basename(files)
  r <- paste0(patterns,"|", collapse = "") %>% str_sub(1,-2)
  bool = str_detect(files, r) 
  names <- files[bool == TRUE]
  
  fullname = paste0(path,"/",names[1])
  d_old <- read.csv(fullname, header = TRUE, sep = ",",
                    colClasses = rep("character",29))
  
  for (n in names[2:length(names)]){
    fullname = paste0(path,"/",n)
    d_new <- read.csv(fullname, header = TRUE, sep = ",", 
                      colClasses = rep("character",29))
    d_old <- bind_rows(d_old,d_new)
  }
  return(d_old)
}

standardise_names <- function(df){
  # ---- INTENDED FOR PRIVATE USE ---
  # Function to standardise names used for the same spp across
  # Zooniverse datasets.
  # Returns a dataframe.
  
  res <- df
  res$question__species <- as.character(res$question__species)
  res$question__species[res$question__species == "birdsofprey"] <- "birdofprey"
  res$question__species[res$question__species == "duiker"] <- "duikercommon"
  res$question__species[res$question__species == "duikercommongrey"] <- "duikercommon"
  res$question__species[res$question__species == "harecape"] <- "hare"
  res$question__species[res$question__species == "catafricanwild"] <- "wildcat"
  res$question__species[res$question__species == "catblackfooted"] <- "blackfootedcat"
  res$question__species[res$question__species == "zebra"] <- "zebraburchells"
  
  res$question__species[res$question__species == "NOT identifiable"] <- "UNIDENTIFIED"
  res$question__species[res$question__species == "Unidentifiable"] <- "UNIDENTIFIED"
  
  return(res)
}

clean_table <- function(df, rm_bycatch = TRUE, pres_abs = TRUE){
  # Clean datatable.
  # rm_bycatch: should the spp defined as bycatch be included?
  # pres_abs: should the data be transformed to presence/absence before aggregating?
  # Returns a S3 custom class "cooc_table" with an attribute data with the cleaned data table.
  
  
  # Merge spp with the same names
  res <- standardise_names(df)
  
  res$capture_date_local <- as.Date(res$capture_date_local)
  res$capture_time_local <- times(res$capture_time_local)
  
  # Convert to numeric
  res$question__count_median <- as.numeric(res$question__count_median)
  
  # Dates and times
  res$capture_date_local <- as.Date(res$capture_date_local)
  res$DateTimeOriginal <- as.POSIXct(res$DateTimeOriginal)
  res$DateTimeOriginal <- force_tz(res$DateTimeOriginal, tz = "Etc/GMT-2")
  
  # Lion
  res$question__species <- as.character(res$question__species)
  res$question__species[res$question__species == "lionmale"] <- "lion"
  res$question__species[res$question__species == "lionfemale"] <- "lion"
  res$question__species[res$question__species == "lioncub"] <- "lion"
  
  # Rabbit
  res$question__species[res$question__species == "rabbitredrock"] <- "rabbit"
  res$question__species[res$question__species == "rabbitriverine"] <- "rabbit"
  
  
  # Fire and co
  res <- filter(res,
               !(question__species %in% NOT_SPP))
  # Remove non mammals 
  if(rm_bycatch){
    res <- filter(res,
                  !(question__species %in% BYCATCH))
  }
  # Convert to presence/absence data
  if(pres_abs){
    res <- mutate(res, 
                  question__count_median = ifelse(question__count_median > 1,
                                                  1,question__count_median))
  }
  
  # Define data table
  res$question__species <- as.factor(res$question__species)
  
  res_class <- cooc_table(data = res, 
                          rm_bycatch = ifelse(rm_bycatch, TRUE, FALSE),
                          pres_abs = pres_abs,
                          threshold = NA, 
                          round_hr = NA, 
                          timestep = NA, 
                          space_scale = NA)
  
  return(res_class)
}

filter_rare_spp <- function(cooc_tab, threshold = 10, filter = TRUE,
                            log = TRUE){
  # Function to filter species with less than (threshold)
  # individuals of one species in the table.
  # df is a cleaned table (S3 cooc-table object).
  # threshold is the lower bound of the counts to keep.
  # if filter = TRUE, returns a filtered table: else, returns
  # a plot of abundances.
  # if log = TRUE, plots abundances in log10.
  
  if(!is(cooc_tab, "cooc_table")){
    stop("df must be a S3 cooc_table.")
  }
  
  df <- cooc_tab$data
  
  d.summary <- df %>% group_by(question__species) %>% 
    summarise(count = sum(question__count_median)) %>% ungroup()
  
  if(filter){ # e want to return the whole object
    d.summary.filtered <- d.summary %>% filter(count > threshold)
    
    r <- filter(df, question__species %in% d.summary.filtered$question__species)
    
    res <- cooc_tab
    res$data <- r
    res$threshold <- threshold
  }
  else{ # We just want to plot
    res <- ggplot(d.summary) +
      geom_point(aes(y = count, x = reorder(question__species, count))) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
    if(log){
      res <- res + scale_y_log10()
    }
  }
  return(res)
}

set_round_hour <- function(datetime, h){
  # ---- INTENDED FOR PRIVATE USE ---
  # Set the time in the datetime object at hour h:00:00.
  # h is an integer between 0 and 23.
  
  hour(datetime) <- h
  minute(datetime) <- 0
  second(datetime) <- 0
  return(datetime)
}

find_closest_datetime <- function(t, datetime_df, datetime_unique){
  # ---- INTENDED FOR PRIVATE USE ---
  # Finds the closest datetime element from the given datetime, but beginnig
  # at a given time of day at one day interval max.
  # t: an integer between 0 and 23 (the hour)
  # datetime_df: a dataframe with columns mini (the datetimes to find a closest date from)
  #   and group_ID (the ID for each datetime)
  # datetime_unique: if provided, evaluates the closest start for a unique datetime.
  
  if(missing(datetime_unique)){
    # --- Get the day before and the day after at the same time.
    datetime_boundaries <- datetime_df %>% dplyr::select(group_ID, mini) %>%
      rename(Original = mini) %>% mutate(After = Original + ddays(1)) %>%
      mutate(Before = Original - ddays(1))
  }
  else if(missing(datetime_df)){
    # --- Get the day before and the day after at the same time.
    Original <- datetime_unique
    After <- Original + ddays(1)
    Before <- Original - ddays(1)
    datetime_boundaries <- data.frame(Original, Before, After)
    
  }
  # --- Round time of day (get starting candidates)
  candidates <- datetime_boundaries %>% 
    mutate(Before = set_round_hour(datetime = Before, h = t)) %>% 
    mutate(Day = set_round_hour(datetime = Original, h = t)) %>% 
    mutate(After = set_round_hour(datetime = After, h = t))
  if(missing(datetime_df)){
    candidates <- candidates %>% dplyr::select(Original, everything())
  }
  else if(missing(datetime_unique)){
    candidates <- candidates %>% dplyr::select(group_ID, Original, everything())
  }
  
  # Transform to long form
  candidates_long <- candidates %>% pivot_longer(cols = c(Before, Day, After), names_to = 'Pos',
                                                 values_to = "Date")
  
  
  # --- Compute lag between original and candidates
  candidates_long <- candidates_long %>% mutate(Duration = abs(Original - Date))
  
  # --- Get minimal duration
  closest <- candidates_long %>% group_by(Original) %>% dplyr::filter(Duration == min(Duration))
  
  # Rename column for clarity
  closest <- closest %>% dplyr::rename(Chosen = Date) %>% ungroup()
  
  return(closest)
}

get_best_compromise <- function(df, step){
  # ---- INTENDED FOR PRIVATE USE ---
  # Get the best compromise between closeness to the original start and exactitude
  # df: a dataframe with cols (group_ID), Original, Pos, Chosen, Duration (output of find_closest_datetime)
  # step: the timestep to use for aggregation
  # Returns a dataframe with added columns Chosen_compromise (the new chosen start date),
  #   Alternative (the date that should be chosen instead if the original is too far) 
  #   and Chosen_before (boolean, is the chosen date before the original start?).
  
  # --- Copy dataframe
  df2 <- df
  
  # --- Is the chosen date before the original start?
  df2 <- df2 %>% mutate(Chosen_before = ifelse(Chosen < Original, TRUE, FALSE))
  
  # --- Get alternative to chosen (always after original start)
  # If chosen is before: the alternative is day +1
  # Else: we don't czre coz it's after.
  df2 <- df2 %>% mutate(Alternative = if_else(Chosen_before, Chosen + ddays(1),
                                              Chosen)) 
  # Rq: if_else =/= ifelse keeps date type
  
  # --- Get the new chosen date
  # Change to alternative iff chosen is before and duration > step/2
  df2 <- df2 %>% mutate(Chosen_compromise = if_else((Chosen_before & (Duration > step/2)),
                                                    Alternative, Chosen))
  
  # --- Reorder columns
  if(df2 %>% has_name('group_ID')){
    df2 <- df2 %>% dplyr::select(group_ID, Original, Chosen_compromise, Pos, Duration, 
                                 everything())
  }
  else{
    df2 <- df2 %>% dplyr::select(Original, Chosen_compromise, Pos, Duration, 
                                 everything())
  }
  
  return(df2)
}

aggregate_table <- function(cooc_tab, timestep = ddays(1), option = 'time',
                            start_hour = 5, use_start = TRUE,
                            hourcuts = times(c('00:00:00', '5:00:00', '13:00:00', '19:00:00','23:59:59'))){
  # Aggregates the table df on the basis of the provided option.
  # df: the S3 cooc_table object which $data should be aggregated
  # timestep: a dxxx object (lubridate), defaults to one day.
  # option: if the data should be aggregated by site/camera too.
  # Possible values: 'site', 'cam', 'sitecam', 'sitetime', 'timecam',
  #     'clust', 'sitetimecam', 'siteclust', 'timeclust', 'sitetimeclust',
  #     ''sitehour', 'hourcam', 'timehour', 'sitehourcam', 'sitetimehour', 
  #     'timehourcam', 'hourclust', 'sitetimehourcam', 'sitehourclust',
  #     'timehourclust',  'sitetimehourclust', 'hour'.
  # start_hour: the hour to cut each day at (defaults to 5), 
  #     an integer between 0 and 23. Cannot be used with hour....
  # use_start: whether to use start hour.
  # hourcuts: the cuts to use for hours (must start at 00:00:00 and end at 23:59:59).
  #     Optional, only useful if aggregation is to be performed by hour.
  
  # If not specified, only time aggregation is performed.
  # Returns the S3 object with the aggregated table
  
  # By time + ... (site or cam): aggregation site or camera-wise
  # By time + site + cam: aggregation site and camera-wise (ie camera wise coz
  # not 2 cameras share the same ID)
  
  
  if(!is(cooc_tab, "cooc_table")){
    stop("df must be a S3 cooc_table.")
  }
  if(!(option %in% c('site', 'cam', 'time', 'sitecam', 'sitetime', 'timecam', 'clust', 
                   'sitetimecam', 'siteclust', 'timeclust', 'sitetimeclust',
                   'sitehour', 'hourcam', 'timehour', 'sitehourcam', 'sitetimehour', 
                   'timehourcam', 'hourclust', 'sitetimehourcam', 'sitehourclust',
                   'timehourclust',  'sitetimehourclust',
                   'hour'))){
    stop("option must be in: 
    c('site', 'cam', 'time, 'sitecam', 'sitetime', 'timecam', 'clust', 'sitetimecam', 
         'siteclust', 'timeclust', 'sitetimeclust', 'sitehour', 'hourcam', 'timehour', 'sitehourcam', 'sitetimehour', 
                   'timehourcam', 'hourclust', 'sitetimehourcam', 'sitehourclust',
                   'timehourclust',  'sitetimehourclust',
                   'hour')")
  }
  
  
  df <- cooc_tab$data
  
  if(grepl('hour', option)){ # First add hour colums in case we aggregate by hour
    # Then we don't want use_start
    use_start = FALSE
    
    # We need to add an hour_ID
    lab <- paste(substr(hourcuts[1:(length(hourcuts) - 1)],1,2),
                 substr(hourcuts[2:(length(hourcuts))],1,2), sep = "--")
    df <- df %>% mutate(hour_ID = cut(capture_time_local, 
                                      breaks = hourcuts, 
                                      labels = lab,
                                      include.lowest = TRUE
    ))
    
    # Collapse firts part of night and last part of night in one factor
    nightlab <- paste(substr(hourcuts[(length(hourcuts) - 1)],1,2),
                      substr(hourcuts[2],1,2), sep = "--")
    levels(df$hour_ID)[length(levels(df$hour_ID))] <- nightlab
    levels(df$hour_ID)[1] <- nightlab
  }
  
  if(grepl('site', option) | grepl('clust', option) | grepl('cam', option)){
    if(grepl('cam', option)){ # we want to agregate by camera (and time)
      df.g <- group_by(df, site_ID)
      group.name <- "site_ID"
    }
    else if(grepl('clust', option)){ # we want to agregate by cluster (and time)
      df.g <- group_by(df, cluster_ID)
      group.name <- "cluster_ID"
    }
    else if(grepl('site', option) & !(grepl('cam', option) | grepl('clust', option))){ # only if 2 above not present, we want to agregate by site (and time)
      df.g <- group_by(df, code_loc)
      group.name <- "code_loc"
    }
    
    # Define time boundaries for each site / camera / cluster
    boundaries_raw <- summarise(df.g, mini = min(DateTimeOriginal),
                                maxi = max(DateTimeOriginal))
    colnames(boundaries_raw) <- c("group_ID", "mini", "maxi")
    
    if(use_start){
      # Find closest start date for each unit
      closest_start <- find_closest_datetime(t = start_hour, datetime_df = boundaries_raw)
      
      # Get compromise for starting date
      compromise <- get_best_compromise(closest_start, timestep)
      compromise_subset <- as.data.frame(compromise) %>% dplyr::select(group_ID, Chosen_compromise) # subset for merge
      # Change boundaries to match desired hour
      
      boundaries <- inner_join(boundaries_raw, compromise_subset, by = 'group_ID') %>% 
        dplyr::select(-mini) %>% rename(mini = Chosen_compromise) %>% 
        dplyr::select(group_ID, mini, maxi)
    }
    else{
      boundaries <- boundaries_raw
    }
    
    boundaries <- column_to_rownames(boundaries, var = "group_ID")
    
    # Define sequence for each site / camera / cluster
    breaks <- by(data = boundaries, INDICES = rownames(boundaries), 
                 FUN = function(x, tstep) seq(from = as.POSIXct(x[[1]], tz = "Etc/GMT-2"), 
                                              to = as.POSIXct(x[[2]], tz = "Etc/GMT-2") + tstep, 
                                              by = tstep),
                 tstep = timestep, simplify = FALSE)

    # Get entire duration of cams (add after breaks calculation)
    boundaries$camdays <- as.duration(boundaries$maxi - boundaries$mini)
    
    # Delete last break if too short
    # -- Compute last intervals for each camera
    # beginning of last interval
    beginnings <- lapply(breaks, function(x) x[length(x) - 1]) # not last one which is max + timestep
    beginnings <- t(as.data.frame(beginnings))
    
    tab <- as.data.frame(beginnings)
    colnames(tab) <- "beginnings"
    
    # End of last interval
    tab <- tab %>% rownames_to_column(var = "group_ID")
    tab$ends <- boundaries[tab$group_ID,"maxi"]
    
    # Cast
    tab$ends <- as.POSIXct(tab$ends, tz = "Etc/GMT-2")
    tab$beginnings <- as.POSIXct(tab$beginnings, tz = "Etc/GMT-2")
    # Get last interval
    tab$intervals <- interval(tab$beginnings,tab$ends)
    
    # -- Compare to the last interval
    min_len <- timestep/2 # minimal length for the last interval
    tab$last_obs_duration <- as.duration(int_length(tab$intervals)) # Duration of last interval
    
    # Add camdays to tab to check duration
    camdays.join <- boundaries %>% rownames_to_column("group_ID") %>% dplyr::select(c("group_ID", "camdays"))
    tab <- full_join(tab, camdays.join,
                     by = "group_ID")
    
    # Check that all cameras will be used
    cams.discarded <- tab$group_ID[tab$camdays < timestep]
    
    if((length(cams.discarded) != 0)){
      w <- paste("Warning: camera(s) /site(s) / cluster(s)", paste(cams.discarded, collapse = ", "), 
                 "was / were discarded")
      print(w)
      tab <- filter(tab, !(group_ID %in% cams.discarded))
    }
    
    # Select sites for which there is a too short time window
    sites.too.short <- tab$group_ID[tab$last_obs_duration < min_len]
    
    breaks[sites.too.short] <- lapply(breaks[sites.too.short], 
                                      function(x) x[-(length(x) -1 )])
    
    sites <- tab$group_ID # get only not discarded sites
    
    for(i in 1:length(sites)){
      s <- sites[i]
      
      if(grepl('cam', option)){ # If groupped by camera
        rsite <- filter(df.g, site_ID == s) %>% dplyr::select(site_ID, DateTimeOriginal)
      }
      else if(grepl('clust', option)){
        rsite <- filter(df.g, cluster_ID == s) %>% dplyr::select(cluster_ID, DateTimeOriginal)
      }
      else if(grepl('site', option) & !(grepl('cam', option) | grepl('clust', option))){
        rsite <- filter(df.g, code_loc == s) %>% dplyr::select(code_loc, DateTimeOriginal)
      }
      
      if(use_start){ # If we use a custom start hour and did not aggregate by hour
        smalltimes <- rsite$DateTimeOriginal[rsite$DateTimeOriginal < breaks[[s]][1]]
        if(length(smalltimes != 0)){ # If value are excluded from the range because the start was rounded
          rsite$DateTimeOriginal[rsite$DateTimeOriginal %in% smalltimes] <- breaks[[s]][1] # Set them to first element of the break instead
        }
      }
      
      rsite$time_ID <- as.POSIXct(cut(rsite$DateTimeOriginal,
                                      breaks = breaks[[s]], include.lowest = TRUE), 
                                  tz = "Etc/GMT-2")
      if(i>1){
        res <- bind_rows(res,rsite)
      }
      else{
        res <- rsite
      }
    }
    res <- inner_join(df, res, by = c(group.name, "DateTimeOriginal"))
  }
  
  # Aggregate by time only: breaks on whole table
  else if(option == 'time' | option == 'timehour'){
    if(use_start){
      # Find closest start date for each unit
      closest_start <- find_closest_datetime(t = start_hour, 
                                             datetime_unique = min(df$DateTimeOriginal))
      
      # Get compromise for starting date
      compromise <- get_best_compromise(closest_start, timestep)
      
      mini <- compromise$Chosen_compromise
    }
    else{
      mini <- min(df$DateTimeOriginal)
    }
    
    breaks = seq(from = mini, to = max(df$DateTimeOriginal) + timestep, 
                 by = timestep) # Go to the timestep after the last day at orig time
    
    # Delete last break if too short
    # -- Compute intervals
    beginnings <- breaks[1:(length(breaks)-1)] # all but last one which is max + timestep
    ends <- c(breaks[2:(length(breaks)-1)], max(df$DateTimeOriginal)) # not first one and replace
    # last one by max of cameras
    
    intervals <- interval(beginnings,ends)
    
    # -- Compare to the last interval
    last_obs_len <- int_length(intervals[length(intervals)]) # Duration of last interval
    min_len <- timestep/2 # minimal length for the last interval
    
    last_obs_duration <- as.duration(last_obs_len) # cast for comparison
    
    if(!(last_obs_duration > min_len)){ # lasts less than timestep/2
      # merge last break with forelast
      breaks <- breaks[-(length(breaks)-1)]
    }
    
    if(use_start){
      smalltimes <- df$DateTimeOriginal[df$DateTimeOriginal < breaks[1]]
      if(length(smalltimes != 0)){ # If value are excluded from the range because the start was rounded
        df$DateTimeOriginal[df$DateTimeOriginal %in% smalltimes] <- breaks[1] # Set them to first element of the break instead
      }
    }
    res <- mutate(df,
                  time_ID = as.POSIXct(cut(DateTimeOriginal,
                                           breaks = breaks),
                                       tz = "Etc/GMT-2")) 
  }
  else{ # site or cam or sitecam or clust or siteclust
    res <- df
  }
  
  if(grepl('hour', option) & grepl('time', option)){ # If we groupped dy time and hour,
    # amongst others
    
    # We want only to keep date ID, time is not delevant
    res$time_ID <- date(res$time_ID)
  }
  
  res_total <- cooc_tab
  
  if(option == 'time' | option == 'hour'){ # Default: group only by time
    if(option == 'time'){
      res <- group_by(res, time_ID, question__species)
    }else if(option == 'hour'){
      res <- group_by(res, hour_ID, question__species)
    }
    ngroups <- 1
    res_total$space_scale <- list(ngroups = ngroups, type = NA)
  }
  else if(option == 'timehour'){ # Group by time and hour
      res <- group_by(res, time_ID, hour_ID, question__species)
    ngroups <- 1
    res_total$space_scale <- list(ngroups = ngroups, type = NA)
  }
  else if(option == "cam"){ # Group by cam
    res <- group_by(res, site_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "hourcam"){ # Group by cam and hour
    res <- group_by(res, site_ID, hour_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "clust"){ # Group by cluster
    res <- group_by(res, cluster_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  else if(option == "hourclust"){ # Group by cluster and hour
    res <- group_by(res, cluster_ID, hour_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  else if(option == "site"){ # Group by study site
    res <- group_by(res, code_loc, question__species)
    ngroups <- length(unique(res$code_loc)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "site")
  }
  else if(option == "sitehour"){ # Group by study site and hour
    res <- group_by(res, code_loc, hour_ID, question__species)
    ngroups <- length(unique(res$code_loc)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "site")
  }
  else if(option == "sitetime"){ # Group by study site and time
    res <- group_by(res, code_loc, time_ID, question__species)
    ngroups <- length(unique(res$code_loc)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "site")
  }
  else if(option == "sitetimehour"){ # Group by study site and time and hour
    res <- group_by(res, code_loc, time_ID, hour_ID, question__species)
    ngroups <- length(unique(res$code_loc)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "site")
  }
  else if(option == "sitecam"){ # Group by study site and camera
    res <- group_by(res, code_loc, site_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "sitehourcam"){ # Group by study site and camera and hour
    res <- group_by(res, code_loc, hour_ID, site_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "timecam"){ # Group by cam and time
    res <- group_by(res, time_ID, site_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "timehourcam"){ # Group by cam and time and hour
    res <- group_by(res, time_ID, hour_ID, site_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "timeclust"){ # Group by clust and time
    res <- group_by(res, time_ID, cluster_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  else if(option == "sitetimeclust"){ # Group by site, clust and time
    res <- group_by(res, code_loc, time_ID, cluster_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  else if(option == "timehourclust"){ # Group by cam and time and hour
    res <- group_by(res, time_ID, hour_ID, cluster_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  else if(option == "sitetimehourclust"){ # Group by cam and time and hour
    res <- group_by(res, code_loc, time_ID, hour_ID, cluster_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  else if(option == "sitetimecam"){ # Group by study site, camera and time
    res <- group_by(res, time_ID, code_loc, site_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "sitetimehourcam"){ # Group by study site, camera and time ad hour
    res <- group_by(res, code_loc, time_ID, hour_ID, site_ID, question__species)
    ngroups <- length(unique(res$site_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "cam")
  }
  else if(option == "siteclust"){ # Group by study site, camera and time
    res <- group_by(res, code_loc, cluster_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  else if(option == "sitehourclust"){ # Group by study site, camera and time and hour
    res <- group_by(res, code_loc, hour_ID, cluster_ID, question__species)
    ngroups <- length(unique(res$cluster_ID)) # Do it now not to count discarded sites
    res_total$space_scale <- list(ngroups = ngroups, type = "clust")
  }
  
  res <- summarise(res, count = sum(question__count_median))
  
  res_total$data <- res
  
  if(!grepl('time', option)){ # if data was not grouped by time
    res_total$timestep = NA
    res_total$round_hour = NA
  }
  else{  # if data was grouped by time
    res_total$timestep = timestep
    
    if(use_start){ # if a constant starting hour was used
      res_total$round_hr = start_hour
    }
    else{
      res_total$round_hr = NA
    }
  }
  
  if(grepl('hour', option)){ # Data was grouppde by hour
    res_total$round_hr = FALSE
    res_total$hour = TRUE
  }
  else{
    res_total$hour = FALSE
  }

  return(res_total)

}

expand_table <- function(cooc_tab, aggregated = TRUE){
  # Assumes the argument is a cooc_table with cleaned data
  # aggregated: whether the data is aggregated or not
  
  if(!is(cooc_tab, "cooc_table")){
    stop("cooc_tab must be a S3 cooc_table.")
  }
  
  df <- cooc_tab$data
  
  if(aggregated){
    res <- df %>% pivot_wider(names_from = question__species,
                              values_from = count,
                              values_fill = list(count = 0))
  }
  else{
    res <- df %>% pivot_wider(names_from = question__species,
                              values_from = question__count_median,
                              values_fill = list(question__count_median = 0))
  }
  
  res_final <- cooc_tab
  res_final$data <- res
  
  return(res_final)
}

sep_daynight <- function(cooc_tab, lat = -26.2, long = 28.0){
  # Compute sunrise/sunset time for the datetime series
  # cooc_tab = the cooc_table object with cleaned table
  
  if(!is(cooc_tab, "cooc_table")){
    stop("cooc_tab must be a S3 cooc_table.")
  }
  
  df <- cooc_tab$data
  
  rinf <- min(df$DateTimeOriginal)
  rsup <- max(df$DateTimeOriginal)
  
  n <- rsup-rinf
  n <- as.numeric(n)
  
  sunset.rise <- sunrise.set(lat, long, rinf, timezone = "Etc/GMT-2", num.days = n+1)
  
  # Coerce to SA timezone
  sunset.rise$sunrise <- force_tz(sunset.rise$sunrise, tz = "Etc/GMT-2")
  sunset.rise$sunset <- force_tz(sunset.rise$sunset, tz = "Etc/GMT-2")
  
  # Make intervals
  days <- interval(sunset.rise$sunrise,sunset.rise$sunset)
  
  res <- vector()
  for(i in 1:length(df$DateTimeOriginal)){
    d <- df$DateTimeOriginal[i]
    day <- days[date(int_start(days)) == date(d)]
    
    verdict <- d %within% day
    if(verdict){
      res <- c(res,"Day")
    }
    else{
      res <- c(res,"Night")
    }
  }
  
  rdf <- mutate(df,
                 Daytime = as.factor(res))
  
  res_final <- cooc_tab
  res_final$data <- rdf
  
  return(res_final)
}

process_data <- function(d, rm_bycatch = TRUE, filter = TRUE, threshold = 1, aggregate = TRUE, 
                         lag = ddays(1), opt, start_hour = 5, use_start = TRUE, 
                         hourcuts = times(c('00:00:00', '5:00:00', '13:00:00', '19:00:00','23:59:59')),
                         sep_dn = FALSE, pres_abs = TRUE, lat = -26.2, long = 28.0, expand = TRUE){
  # rm_bycatch : rm mammals before cleaning?
  # filter: should the data be filtered (delete spp with low abundances)
  # threshold: threshold of abundances (spp with less than threshold sightings will not
  #   be taken into account)
  # aggregate: whether the data should be aggregated
  # lag: the time lag for aggregation (defaults to 1 day)
  # opt: aggregate by ("cam" or "site"). If missing, aggregation only based on time.
  # sep_dn: should day and night be separated?
  # pres_abs: should the data be converted to presence absence first?
  # lat / long: Johannesburg coordinates
  # expand: return an expanded table?
  # ---------
  # Returns a named list "Day" and "Night" if separate dn is true, else returns a tibble.
  
  dc <- clean_table(d, rm_bycatch = rm_bycatch, pres_abs = pres_abs)
  
  if(filter){
    if(!sep_dn){
      dc <- filter_rare_spp(dc, threshold = threshold)
    }
  }
  
  if(sep_dn){ # Separate day/night
    dc <- sep_daynight(dc, lat, long)
    # Copy object in 2 separate objects
    dcd <- dc
    dcn <- dc
    
    dcd$data <- filter(dc$data, Daytime == "Day")
    dcd <- filter_rare_spp(dcd, threshold = threshold)
    
    dcn$data <- filter(dc$data, Daytime == "Night")
    dcn <- filter_rare_spp(dcn, threshold = threshold)
    
    if(aggregate){
      if(missing(opt)){
        dcd <- aggregate_table(dcd, lag, use_start = FALSE)
        dcn <- aggregate_table(dcn, lag, use_start = FALSE)
      }
      else{
        dcd <- aggregate_table(dcd,lag, option = opt, use_start = FALSE)
        dcn <- aggregate_table(dcn,lag, option = opt, use_start = FALSE)
      }
      if(expand){
        de <- expand_table(dcd)
        ne <- expand_table(dcn)
      }
    }
    else if(expand){
      de <- expand_table(dcd, aggregated = FALSE)
      ne <- expand_table(dcn, aggregated = FALSE)
    }
    res <- de
    res$data <- list(Day = de$data,
                     Night = ne$data)
    # res <- list("Day" = de,
    #             "Night" = ne)
  }
  
  else{  # Don't sep day/night
    if(aggregate){
      if(missing(opt)){
        res <- aggregate_table(dc,lag, start_hour = start_hour, use_start = use_start, hourcuts = hourcuts)
      }
      else{
        res <- aggregate_table(dc,lag, option = opt, start_hour = start_hour, use_start = use_start, hourcuts = hourcuts)
      }
      if(expand){
        res <- expand_table(res)
      }
      
    }
    else{ # if don't aggregate
      if(expand){
        res <- expand_table(dc, aggregated = FALSE)
      }
      else{
        res <- dc
      }
    }
    
  }
  return(res)
}

rebind_dn <- function(list){
  # Make a single table from a list with day/night partition
  d <- mutate(list$Day, Daytime = "Day")
  n <- mutate(list$Night, Daytime = "Night")
  
  res <- bind_rows(d, n)
  res[is.na(res)] <- 0 # consider that if a spp was not seen it was absent
  
  # Reorder
  res <- dplyr::select(res, Daytime, everything())
  return(res)
}

estimate_lag <- function(df, n = 100, thr){
  # df: a count df (ie only spp counts)
  # n = repetitions nb
  # thr: the richness threshold to reach
  # ------
  # Returns the vector of number of observations that were needed eacu repetition to reach thr
  # (as a proxy for time lag needed to see the required richness)
  R <- rowSums(df != 0)
  n_obs <- rep(0,n)
  for(rep in 1:n){
    cumsum <- 0
    while(cumsum < thr){
      i <- sample(1:length(R), 1, replace=TRUE)
      cumsum <- cumsum + R[i]
      n_obs[rep] <- n_obs[rep] +1
    }
  }
  return(n_obs)
}

join_clusters_to_df <- function(clusters, df){
  # Merges clusters as returned by cutree on the output of spatial_clustering_from_sfc
  #   with a snapshot dataframe.
  # Joins by site_ID.
  # clusters: a named list (names = site_ID and content = int for cluster ID)
  # Returns a merged dataframe.
  # Assumes that df has a site_ID column.
  
  # Build a dataframe with clusters
  site_ID <- names(clusters)
  cluster_ID <- clusters
  clusters.df <- data.frame(site_ID, cluster_ID)
  # Coerce factor to character
  clusters.df$site_ID <- as.character(clusters.df$site_ID)
  df$site_ID <- as.character(df$site_ID)
  
  warning <- unique(df$site_ID[!(df$site_ID %in% site_ID)])
  if(length(warning) != 0){ # Then it means at least one site that is in the df is not in the clustering
    print(paste("Site(s)", paste(warning, collapse = ", "), 
                "will be discarded because they are not in the clustering."))
    print("Hint: did you filter for a minimal captures nb? Are they cameras that have no coordinates?")
  }
  # Merge
  m <- inner_join(df, clusters.df, by = "site_ID")
  # Reorder columns
  res <- m %>% dplyr::select(code_loc, cluster_ID, site_ID, everything())
  return(res)
}

merge_spatial_info <- function(pts, vegetation){
  # Merges info in pts and vegetation
  # pts: a df with columns "long_X_deg","lat_Y_deg"to be used as points coordinates
  # vegetation: the spatial info to add to these points
  
  # Returns a spatial object  
  
  pt2 <- st_as_sf(pts, coords = c("long_X_deg","lat_Y_deg"))
  pt2 <- pt2 %>% st_set_crs("+proj=longlat +ellps=WGS84")
  
  p <- st_join(pt2, vegetation, join = st_intersects)
  return(p)
}

plot_bar <- function(df, x="time", ...){
  # Plot barchart of observed species by aggregate unit.
  # Data is plotted along time, grouped by site_ID / code_loc if relevant.
  # x can be space or time (indicates x axis to use). Useful when
  # comparing multiple sites for one time value (bc all bars are stacked anyway).
  # ... = graphical options to pass to ggplot
  
  # Define x axis
  if(x=="time"){
    gg <- ggplot(df, aes(x = time_ID, y = count, fill = question__species))
    if("site_ID" %in% colnames(df)){ # The table was aggregated on site
      gg <- gg+ facet_wrap(~site_ID)
    }
    else if("code_loc" %in% colnames(df)){
      gg <- gg+ facet_wrap(~code_loc)
    } 
  }
  else if(x == "site"){
    gg <- ggplot(df, aes(x = code_loc, y = count, fill = question__species)) +
      theme(axis.text.x = element_text(angle = 90))
  }
  else if(x == "cam"){
    gg <- ggplot(df, aes(x = reorder(site_ID, -count), y = count, fill = question__species)) +
      theme(axis.text.x = element_text(angle = 90))
  }
  
  gg <- gg + geom_bar(stat = "identity") + theme(legend.text = element_text(size = 10),
                                                 legend.title = element_text(size = 12),
                                                 legend.key.size = unit(.2,"cm")) +
    guides(fill=guide_legend(ncol=1))
  if(!missing(...)){
    gg <- gg + ...
  }
  return(gg)
}