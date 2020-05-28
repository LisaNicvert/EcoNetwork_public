########################################################
# Helpers for multivariate analyses
########################################################

library(ade4)
library(adegraphics)
library(dplyr)

pcoa <- function(df, aggregated = "time", meth = 1, 
                 scan = FALSE, do_pco = TRUE){
  # Performs PCoA on datatable.
  # df is an expanded datatable.
  # Aggreted is a keyword to specify the type of aggregation performed:
  # ("no" if none performed).
  # meth is an integer used to determine the 
  # distance method to use (default: 1 = Jaccard) => integers of ade4.
  # scan: whether to scan for eigenvalues (default: false, 2 eigenvalues kept).
  # do_pco: whether the pcoa should be perfrmed, or only the distance matrix returned
  
  if(aggregated == "no"){
    dat <- select(df, -c(question__count_min, question__count_max,
                         capture_id, site_ID, site, code_loc, season1, roll,
                         capture, DateTimeOriginal, capture_date_local,
                         capture_time_local))
  }
  else if(aggregated == "time"){
    df <- as.data.frame(df)
    dat <- select(df, -time_ID)
    rownames(dat) <- df$time_ID
  }
  else if(aggregated == "dncam"){
    df <- as.data.frame(df)
    dat <- select(df, -c(Daytime, site_ID))
    rownames(dat) <- paste(substr(df$Daytime,1,1), substr(df$site_ID,5,7))
  }
  else if(aggregated == "sitetime"){
    df <- as.data.frame(df)
    dat <- select(df, -c(time_ID, code_loc))
    rownames(dat) <- paste(df$code_loc, df$time_ID)
  }
  else if(aggregated == "timecam"){
    df <- as.data.frame(df)
    dat <- select(df, -c(time_ID, site_ID))
    rownames(dat) <- paste(substr(df$site_ID,5,7), df$time_ID)
  }
  else if(aggregated == "sitecamtime"){
    df <- as.data.frame(df)
    dat <- select(df, -c(time_ID, site_ID,
                         code_loc))
    rownames(dat) <- paste(df$site_ID, df$time_ID)
  }
  else if(aggregated == "sitecam"){
    df <- as.data.frame(df)
    dat <- select(df, -c(site_ID,
                         code_loc))
    rownames(dat) <- df$site_ID
  }
  else if(aggregated == "cam"){
    df <- as.data.frame(df)
    dat <- select(df, -site_ID)
    rownames(dat) <- substr(df$site_ID,5,7)
  }
  di <- dist.binary(dat, method = meth)
  if(do_pco){
    pco <- dudi.pco(d = di, scannf = scan)
    return(pco)
  }
  else{
    return(di)
  }
  
}

plot_pcoa <- function(pco, df, col = "cam", label = TRUE, 
                      class = FALSE, fac.c,
                      traject = FALSE, fac.t,
                      plot.p.c = TRUE,
                      lab = rownames(pco$tab), 
                      ...){
  # Takes in argument a pco object and plots it nicely.
  # pco: a pco object (ade4)
  # df: the complete dataframe associated with thze pco with covariates
  # col: the colors to use. Can be colored either by site, cam, date or month.
  # lab: the labels for the points (default: NULL)
  # label: plot s.label
  # class: plot s.class (requires fac.c)
  # traject: plot s.traject (requires fac.t)
  # plot.p.c : whether to plot points in the class
  # ... adegraphics graphical parameters.
  
  if(label){ # define colors for the s.label
    if(col == "cam"){ # For cameras
      pal <- rainbow(length(unique(df$site_ID)))
      cols <- pal[as.factor(df$site_ID)]
    }
    else if(col == "time"){ # For time
      pal <- rainbow(length(unique(df$time_ID)))
      cols <- pal[as.factor(df$time_ID)]
    }
    else if(col == "month"){ # For month
      mon <- month(date(df$time_ID))
      pal <- rainbow(length(unique(mon)))
      cols <- pal[as.factor(mon)]
    }
    else if(col == "site"){ # For sites
      pal <- rainbow(length(unique(df$code_loc)))
      cols <- pal[as.factor(df$code_loc)]
    }
    else if(col == "dn"){ # For day/night
      pal <- rainbow(length(unique(df$Daytime)))
      cols <- pal[as.factor(df$Daytime)]
    }
    
    if(col != "no"){
      la <- s.label(pco$li, labels = lab,
                    ppoints.col = cols, plabels.optim = TRUE,
                    plabels.col = cols, plot = FALSE,
                    ...)
    }
    else{
      la <- s.label(pco$li, labels = lab, 
                    plabels.optim = TRUE,
                    plot = FALSE, ...)
    }

    res <- la
  }
  
  if(class){
    if(plot.p.c){
      cla <- s.class(pco$li, fac.c,
                     ellipse = 0, starSize = 0, chullSize = 1,
                     col = TRUE, plabels.optim = TRUE,
                     plot = FALSE)
    }
    else{
      cla <- s.class(pco$li, fac.c,
                     ellipse = 0, starSize = 0, chullSize = 1,
                     ppoints.cex = 0,
                     col = TRUE, plot = FALSE)
    }
    
    if(exists("res")){
      res <- res+cla
    }
    else{
      res <- cla
    }
    
  }
  if(traject){
    traj <- s.traject(pco$li, fac.t,
                      col = TRUE, plot = FALSE)
    if(exists("res")){
      res <- res+traj
    }
    else{
      res <- traj
    }
  }
  
  res
}