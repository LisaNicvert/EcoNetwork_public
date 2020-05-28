# ----------------------------
# Function call to render Rmarkdown document
# ----------------------------

# ------------------------------------------ Parameters ----------------------------------------
# folder <- "../filter/filtered_spa/not_annoying/"
# fname <- "KGA_indep_records_30mins_2.csv"

folder <- "../filter/filtered_spa/not_annoying/seasons/"
fname <- "MTZ_dry.csv"

th <- 4
sel <- "BIC"
d <- c(lubridate::ddays(20),lubridate::ddays(30),lubridate::ddays(40),lubridate::ddays(60))

use_default_offset = FALSE
habcov = FALSE

use_hour = FALSE
hour_offset = "none"

hour_slices <- chron::times(c('00:00:00', '5:00:00', '9:30:00', 
                                '14:30:00', '19:00:00','23:59:59'))
spatial_scales <- c("cam",
                    "clust_fine")

sfolder <- "test/"

# ----------------------------------------------------------------------------------------------

# -------------------------------------- Function definition -----------------------------------
render_report = function(file, sfolder, 
                         thr, selcrit = "StARS", 
                         durations = c(lubridate::ddays(2), lubridate::ddays(5), 
                                       lubridate::ddays(10)), 
                         use_default_offset = FALSE, use_habcov = FALSE, 
                         use_hour = FALSE, 
                         hour_slices = chron::times(c('00:00:00', '5:00:00', '9:30:00', '14:30:00', '19:00:00','23:59:59')),
                         spatial_scales = c("site", "site_all", 
                                            "cam", "cam_all",
                                            "clust_coarse", "clust_coarse_all",
                                            "clust_fine", "clust_fine_all"),
                         hour_offset){
  if(length(file) > 1){
    codeloc <- substr(basename(file), 1, 3) # Extract codeloc
    ncodeloc = paste(codeloc, collapse = "_")
  }else{
    codeloc <- substr(basename(file), 1, 3) # Extract codeloc
    ncodeloc = codeloc
  }
  
  habcov <- ifelse(use_habcov, "_habcov", "")
  offset <- ifelse(use_default_offset, "_offset", "")
  hour <- ifelse(use_hour, "_byhour", "")

  saving_folder <- paste0(sfolder, ncodeloc, "_thr", thr, habcov, hour, "_", selcrit) # saving folder
  
  # Test if folder esists
  if (file.exists(saving_folder)){ # If so, create an other folder with suffix 
    s <- saving_folder
    nsuffix <- 1
    while(file.exists(s)){
      s <- paste0(saving_folder, "_", nsuffix)
      nsuffix <- nsuffix + 1
    }
    saving_folder <- paste0(s, "/")
  }else{ # Else, just add "/"
    saving_folder <- paste0(saving_folder, "/")
  }
  
  suffix <- ifelse(exists("nsuffix"), paste0("_", nsuffix - 1), "") # Create report suffix
  
  dir.create(saving_folder)
  
  rmarkdown::render(
    "script_space_part.Rmd", params = list(
      file = file,
      thr = thr,
      selcrit = selcrit,
      durations = durations,
      codeloc = codeloc,
      saving_folder = saving_folder,
      use_default_offset = use_default_offset,
      use_habcov = use_habcov,
      use_hour = use_hour,
      hour_slices = hour_slices,
      spatial_scales = spatial_scales,
      hour_offset = hour_offset
    ),
    output_file = paste0("space_part_", ncodeloc, "_thr", thr, habcov, offset, hour, "_", selcrit, suffix, ".html"),
    output_dir = saving_folder
  )
}

# ----------------------------------------------------------------------------------------------

# ---------------------------------------- Function call ---------------------------------------
render_report(file = paste0(folder, fname),  sfolder = sfolder,
              thr = th, selcrit = sel,
              durations = d, 
              use_default_offset = use_default_offset, use_habcov = habcov,
              use_hour = use_hour, hour_slices = hour_slices,
              spatial_scales = spatial_scales, hour_offset = hour_offset)
          
# --------------------------------------- Threshold choice -------------------------------------
source("../functions/data_processing.R")
source("../functions/cooc_classes.R")

# To be run before to choose threshold
k <- read.csv(paste0(folder, fname),
                header = TRUE, sep = ",",  colClasses = rep("character", 29))
# k <- read_files(patterns = c("MAD", "PLN"), path = folder)
k <- clean_table(k)
filter_rare_spp(k, filter = FALSE) + geom_hline(yintercept = th)

rm(k)
# ----------------------------------------------------------------------------------------------