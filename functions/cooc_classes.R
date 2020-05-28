########################################################
# Functions for custom classes
########################################################

# Define coocc_object class
cooc_object <- function(data, rm_bycatch, pres_abs, threshold = NA, 
                        round_hr = NA, timestep = NA, space_scale = NA,
                        hour = NA){
  co <- structure(list(data = data, 
                       rm_bycatch = rm_bycatch,
                       pres_abs = pres_abs,
                       threshold = threshold,
                       round_hr = round_hr,
                       timestep = timestep,
                       space_scale = space_scale,
                       hour = hour),
                  class = "cooc_object")
  return(co)
}

# Define cooc_table class
cooc_table <- function(data, rm_bycatch, pres_abs, threshold = NA, 
                       round_hr = NA, timestep = NA, space_scale = NA, hour = NA){
  ct <- cooc_object(data = data, rm_bycatch = rm_bycatch, pres_abs = pres_abs, 
                    threshold = NA, round_hr = NA, timestep = NA, space_scale = NA, 
                    hour = NA)
  class(ct) <- append(class(ct),"cooc_table")
  
  return(ct)
}

pln_to_tbl_graph <- function(pln_nk){
  # Transforms a PLN graph to a tbl_graph.
  # Also deletes confusing attributes added by PLN.
  
  nk <- plot(pln_nk, plot = FALSE)
  nk <- nk %>% as_tbl_graph()
  
  nk <- nk %>% activate(edges) %>% dplyr::select(-c(width, color))
  
  nk <- nk %>% activate(nodes) %>% dplyr::select(-c(label.cex, label, size, 
                                                    label.color, 
                                                    frame.color))
  return(nk)
}

cooc_network <- function(cooc_data, pln_nk){
  # Constructor foc cooc_nk object.
  # cooc_data: a cooc_data object with data slot transformed by PLN (has an Abundance table)
  # pln_nk: a PLNnetwork object iferred by getBestModel on the data of cooc_data.
  
  cn <- cooc_data
  cn$pln_network <- pln_nk
  nk <- pln_to_tbl_graph(pln_nk)
  
  abundance <- colSums(cooc_data$data$Abundance)
  ab <- data.frame(abundance) %>% rownames_to_column("name")
  
  nk <- nk %>% activate(nodes) %>% left_join(ab, by = "name")
  
  cn$network <- nk
  class(cn) <- c("cooc_obj" ,"cooc_network")
  
  return(cn)
}

# Summary method for cooc_object
summary.cooc_object <- function(cooc_obj){
  class <- class(cooc_obj)[length(class(cooc_obj))]
  print(paste("A", class, "object."))
  if(!is.na(cooc_obj$pres_abs) & cooc_obj$pres_abs){
    print("Data was passed to presence_absence.")
  }
  if(!is.na(cooc_obj$rm_bycatch)){
    print("Bycatches species were discarded.")
  }
  if(!is.na(cooc_obj$threshold)){
    print(paste("Species with less than",cooc_obj$threshold, "individuals were filtered."))
  }
  if(!is.na(cooc_obj$timestep)){
    msg_agg <- paste("Aggregation time:",cooc_obj$timestep)
    if(!is.na(cooc_obj$round_hr)){
      print(paste(msg_agg, "starting at", cooc_obj$round_hr, "h."))
    }
    else{
      print(msg_agg)
    }
  }
  else{
    print("No time aggregation.")
  }
  if(length(cooc_obj$space_scale) > 1){
    if(!(is.na(cooc_obj$space_scale)[1])){
      print(paste("Aggregation scale:", paste(unlist(cooc_obj$space_scale), collapse = " ")))
    }
    else{
      print("No space aggregation.")
    }
  }
  if(!is.na(cooc_obj$hour)){
    print(paste("Aggregated by hour:", cooc_obj$hour))
  }
}

summary.cooc_table <- function(ct){
  summary.cooc_object(ct)
  print("Head of data:")
  print(head(ct$data))
}

summary.cooc_network <- function(cnk){
  summary.cooc_object(cnk)
  print(head(cnk$network))
}

set_network <- function(cooc_nk, nk){
  UseMethod("set_network")
}

set_network.cooc_network <- function(cooc_nk, nk){
  if(!(is(nk, 'igraph') | is(nk, 'tbl_graph') | is(nk, 'PLNnetwork'))){
    stop("nk must be a network.")
  }
  cooc_nk$network <- nk
  return(cooc_nk)
}