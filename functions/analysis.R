########################################################
# Functions used for final network analysis and tweaks
########################################################

library(RColorBrewer)
library(dplyr)

# read_RDSlist <- function(files, path, codeloclen = 3){
#   # Reads RDS files (files) into an object list.
#   # Path: path to RDS (optional).
#   # codeloclen: length of codeloc id to remove from list name 'defaults to 3).
#   l <- vector("list", length(files))
#   for(f in 1:length(files)){
#     if(missing(path)){
#       x <- readRDS(files[f])
#     }
#     else{
#       x <- readRDS(paste0(path, files[f]))
#     }
#     l[[f]] <- x
#   }
#   # Naming list
#   fnames <- sapply(files, FUN = function(x) substr(x, 1, nchar(x) - 4)) # Delete .rds extension
#   
#   fnames <- sapply(fnames, FUN = function(x) substr(x[1], codeloclen + 7, nchar(x[1]))) # substrinf XXX_cooc_
#   
#   names(l) <- fnames
#   return(l)
# }

# read_RDS_in_folder <- function(folder, codeloclen = 3){
#   # Read only RDS files in the folder 'folder'.
#   # folder is the folder to read in.
#   # Codeloclen is the length of the region identifier in the file name 
#   #     (typically KAR => 3).
#   # Returns an RDS objects list.
#   
#   files <- list.files(path = folder)
#   files <- files[grepl(".rds", files, fixed = TRUE)]
#   
#   grlist <- read_RDSlist(files = files, path = folder, codeloclen = codeloclen) # Read all RDS
#   
#   return(grlist)
# }

get_same_space_scale <- function(nklist, scale, clust_scale){
  # Get the same spatial for a list of cooc_network objects.
  # nklist is a list of cooc_network obkects.
  # scale is the spatial scale we want to get.
  # clust scale is the cluster scale to use, if scale = clust; namely coarse of fine.
  
  # Get list of nks of the same scale
  lscale <- lapply(nklist, function(x){x$space_scale[x$space_scale$type == scale]})
  
  # Discard list of length zero
  cond <- sapply(seq_along(lscale), function(i){length(lscale[[i]]) != 0})
  relevant <- names(lscale[cond]) # Names of correct scale nks, and which are not empty
  
  # Nk list of same scale
  r <- nklist[relevant]
  
  if(scale == 'clust'){ # Particular case clust :
    lscale.relevant <- lscale[relevant] # get spatial scale for which list length is not zero
    nclust <- unlist(sapply(lscale.relevant, function(x){x$ngroups})) # extract the clusters counts
    
    if(clust_scale == 'fine'){ # we want the smaller clusters => more clusters
      nclust.desired <- max(nclust)
    }else if(clust_scale == 'coarse'){ # we want the bigger ones
      nclust.desired <- min(nclust)
    }else{
      stop("Cluster scale badly specified.")
    }
    
    cond2 <- sapply(lscale.relevant, function(x){x$ngroups == nclust.desired}) # Select scales with correct clusters count
    relevant2 <- names(lscale.relevant[cond2]) # get their names
    res <- r[relevant2] # senect the networks
    
  }else{
    res <- r
  }
  return(res)
}

get_guild_colours <- function(traits, level = "guild"){
  # Input: traits, a df with columns (at least) specif_guild3 and guild_realm.
  # Level: whether to return colours for guilds or guild_realm
  
  # Output: a named list of hexadecimal colours with level names.
  
  # N guilds in each class
  if(level == "guild"){
    col.summary <- traits %>% group_by(guild_realm) %>% 
      summarise(count_guild = length(unique(specif_guild3)))
    
    # Get n
    nherbi <- col.summary$count_guild[col.summary$guild_realm == "Herbivore"]
    ncarni <- col.summary$count_guild[col.summary$guild_realm == "Carnivore"]
    nomni <- col.summary$count_guild[col.summary$guild_realm == "Omnivore"]
    
  }else if(level == "diet"){
    col.summary <- traits %>% group_by(guild_realm) %>% 
      summarise(count_diet = length(unique(guild_realm)))
    
    # Get n
    nherbi <- col.summary$count_diet[col.summary$guild_realm == "Herbivore"]
    ncarni <- col.summary$count_diet[col.summary$guild_realm == "Carnivore"]
    nomni <- col.summary$count_diet[col.summary$guild_realm == "Omnivore"]
  }
  
  # Colors for each guild
  if(length(nherbi) != 0){
    greens <- c("#85d482", "#177f3f", "#136934", 
                "#0fa721", "#80e560","#59ab5f")
    cols.herbi <- greens[1:nherbi]
    
    if(level == "guild"){
      names(cols.herbi) <- unique(traits$specif_guild3[traits$guild_realm == "Herbivore"])
    }else{
      names(cols.herbi) <- "Herbivore"
    }
  }
  if(length(ncarni) != 0){
    reds <- c("#900a0a", "#d4271e", "#9a3e3e")
    cols.carni <- reds[1:ncarni]
    
    if(level == "guild"){
      names(cols.carni) <- unique(traits$specif_guild3[traits$guild_realm == "Carnivore"])
    }else{
      names(cols.carni) <- "Carnivore"
    }
  }
  if(length(nomni) != 0){
    cols.omni <- "#6495ed" # manual because one class only
    if(level == "guild"){
      names(cols.omni) <- unique(traits$specif_guild3[traits$guild_realm == "Omnivore"])
    }else{
      names(cols.omni) <- "Omnivore"
    }
  }
  
  if(exists("cols.omni") & exists("cols.herbi") & exists("cols.carni")){
    colours <- c(cols.herbi, cols.carni, cols.omni)
  }else{
    colours <- c(cols.herbi, cols.carni)
  }
  
  return(colours)
}

add_guilds_and_reorder <- function(adj, 
                                   guilds = c("Medium carnivores", "Omnivores", "Small carnivores",
                                              "Water-dependent grazers", 
                                              "Small nonsocial browsers", "Large carnivores",  
                                              "Medium-sized social mixed diet", "Extra-large browsers",
                                              "Large browsers", "Nonruminants")){
  # Reorder a matrix and add missing rows if needed to match all guiilds.
  # adj: an adjacency matrix
  # guilds: all col/rownames that must be present.
  
  # Find missing cols/rows
  not_present <- guilds[!(guilds %in% colnames(adj))]
  
  if(length(not_present) != 0){
    # Create rows of NAs
    rows_to_add <- matrix(rep(NA, ncol(adj)*length(not_present)), nrow = length(not_present))
    rownames(rows_to_add) <- not_present
    
    m.rows <- rbind(adj, rows_to_add)
    
    cols_to_add <- matrix(rep(NA, length(not_present)*length(guilds)), ncol = length(not_present))
    colnames(cols_to_add) <- not_present
    
    m.cols <- cbind(m.rows, cols_to_add)
  }else{
    m.cols <- adj
  }
  
  
  res <- m.cols[sort(rownames(m.cols)), sort(colnames(m.cols))]
  
  return(res)
}