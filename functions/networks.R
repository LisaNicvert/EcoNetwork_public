########################################################
# Functions used for network inference and plotting
########################################################

library(igraph)
library(PLNmodels)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(econetwork)


prepare_mydata <- function(cooc_tab, ncols, addsite = FALSE, cov, offset, 
                           merge_cov_by = 'site_ID', merge_offset_by = 'cluster_ID'){
  # Prepare data for PLN analysis.
  # Input : an cooc_tab with expanded df in slot data.
  # ncols : number of covariates columns.
  # addsite: whether the site should be added to the covariates table
  # (default : false)
  # cov: dataframe of covariates with names as the site IDs in df and the covariate(s) as columns.
  # This assumes the covariates are the first colums.
  # Returns a cooc_table with an abundance dataframe in slot data.
  
  if(!is(cooc_tab, "cooc_table")){
    stop("df must be a S3 cooc_table.")
  }
  
  df <- cooc_tab$data
  
  df.df <- as.data.frame(df)
  
  if(!missing(cov)){
    df.cov <- df.df %>% mutate(cov_ID = !!rlang::sym(merge_cov_by))
    
    cov.cov <- cov %>% rownames_to_column('cov_ID')
    
    discarded <- cov.cov$cov_ID[!(cov.cov$cov_ID %in% df.cov$cov_ID)]
    NA.clust <- unique(df.cov$site_ID[!(df.cov$cov_ID %in% cov.cov$cov_ID)])
    
    if(length(discarded) > 0){
      message(paste("Warning : site(s)", paste(discarded, collapse =", "), "is/are in the covariates but not in the data and will be discarded."))
    }
    if(length(NA.clust) > 0){
      message(paste("Warning : site(s)", paste(NA.clust, collapse =", "), "has/have no associate covariates and will be discarded."))
    }
    
    # Join covariates and data, thus removing unused offsets and not offsetted data
    df.cov <- cov.cov %>% inner_join(df.cov, by = 'cov_ID') %>% droplevels()
    
    covariates <- df.cov[,(1:(ncol(cov.cov)+ncols))] # Extract covariates (n first rows)
    covariates <- covariates %>% dplyr::select(-cov_ID)

    df.df <- df.cov[,-(1:ncol(cov.cov))] # Extract data
  }
  else{
    covariates <- df.df[,1:ncols]
  }
  
  if(!missing(offset)){
    df.offset <- df.df %>% mutate(offset_ID = !!rlang::sym(merge_offset_by))
    
    off <- offset %>% rownames_to_column("offset_ID")
    colnames(off)[2] <- "Freq"
    
    discarded.offset <- off$offset_ID[!(off$offset_ID %in% df.offset$offset_ID)]
    NA.clust.offset <- unique(df.offset$offset_ID[!(df.offset$offset_ID %in% off$offset_ID)])
    
    if(length(discarded.offset) > 0){
      message(paste("Warning : site(s)", paste(discarded.offset, collapse =", "), "is/are in the offsets but not in the data and will be discarded."))
    }
    if(length(NA.clust.offset) > 0){
      message(paste("Warning : site(s)", paste(NA.clust.offset, collapse =", "), "has/have no associate offsets and will be discarded."))
    }
    
    # Join offset and data, thus removing unused offsets and not offsetted data
    df.offset <- off %>% inner_join(df.offset, by = "offset_ID")  %>% droplevels()
    
    offs <- df.offset[,1:ncol(off)] # Extract offsets (n first rows)
    offs <- offs %>% dplyr::select(-offset_ID) # Remove offset ID
    offset_f <- offs$Freq
    names(offset_f) <- rownames(offs)
    
    df.df <- df.offset[,-(1:ncol(off))] # Extract data
  }
  
  ab <- df.df[,-(1:ncols)]

  if(addsite){
    covariates$site <- substr(covariates$site_ID,1,3) # Add site name
  }
  
  if(missing(offset)){
    res <- prepare_data(counts = ab,
                        covariates = covariates)
  }
  else{
    res <- prepare_data(counts = ab,
                        covariates = covariates,
                        offset = offset_f)
  }
  
  res_total <- cooc_tab
  res_total$data <- res
  
  return(res_total)
}

plot_network <- function(nk, nodes_colour, base_layout, plot, s, labcex,
                         edge_width_coeff, check_overlap, geom_edges, label, legend,
                         colourname){
  UseMethod("plot_network")
} 

plot_nk <- function(nk, nodes_colour, base_layout, 
                    plot, s, labcex, edge_width_coeff,
                    check_overlap, geom_edges, label, legend,
                    colourname){
  # ############ INTENDED FOR PRIVATE USE ##########################
  # nk is a cooc_network object.
  # nodes_colour: named list with names corresponding to nodes colour scheme.
  #     node colour attribute must be present in the df.
  #     value can be either a factor (or coercible), or an hexadecimal color name.
  # plot: boolean (if true then the nk is plotted)
  # s and labcex: resp. size and label multiplication coefficients
  #     for abundance
  # edge_width_coeff: coeff to apply to edge width
  # check_overlap: boolean (whether to optimise labels positions)
  # Returns ggraph plot or tbl_graph object.
  # label: named one column df for labels (rownames correspond to network nodes names)
  
  
  ink <- nk
  
  # Replace network attributes
  # --- Vertices attributes proportional to (relative) abundances (except if too small)
  if("avg_abundance" %in% names(vertex_attr(ink))){ # Use avg abundance if present
    # Plotting avg abundance of a group
    ink <- ink %>% activate(nodes) %>% 
      mutate(spp_norm = avg_abundance/max(avg_abundance))
    
  }else if("abundance" %in% names(vertex_attr(ink))){ # Use abundance else
    ink <- ink %>% activate(nodes) %>% 
      mutate(spp_norm = abundance/max(abundance))
  }else{
    ink <- ink %>% activate(nodes) %>% 
      mutate(spp_norm = 1)
  }
  
  if("avg_weight" %in% names(edge_attr(ink))){ # Use avg weight if present
    # Plotting avg weight of a group
    ink <- ink %>% activate(edges) %>% 
      mutate(width = avg_weight)
    
  }else if("weight" %in% names(edge_attr(ink))){ # Use abundance else
    ink <- ink %>% activate(edges) %>% 
      mutate(width = abs(weight))
  }else{
    ink <- ink %>% activate(edges) %>% 
      mutate(width = 1)
  }
  
  # Add nodes labels
  if(missing(label)){
    ink <- ink %>% activate(nodes) %>% 
      mutate(label = name)
  }else{
    lab <- label %>% rename("label" = 1) 
    ink <- ink %>% activate(nodes) %>% 
      left_join(lab %>% rownames_to_column("name"), by = "name")
  }
  
  # Edge width range
  w <- ink %>% activate(edges) %>% dplyr::select(width) %>% as.data.frame()
  width <- w$width
  edge_width_range = c(min(width), max(width))*edge_width_coeff
  
  # Node size range
  ns <- ink %>% activate(nodes) %>% dplyr::select(spp_norm) %>% as.data.frame()
  nodesize <- abs(ns$spp_norm)
  node_size_range <- c(ifelse(min(nodesize) < .5, .5, min(nodesize)), max(nodesize))*10*s
  
  if(!is.null(base_layout)){
    ink.layout <- create_layout(ink, layout = get_my_layout, base_layout = base_layout)
  }
  else{
    ink.layout <- create_layout(ink, 'circle')
  }
  
  if(missing(geom_edges)){
    if("avg_abundance" %in% names(vertex_attr(ink))){ # Means the graph was aggregated at
      # Some point so defalts to using fam
      geom_edges = 'fan'
    }else{
      geom_edges = 'diagonal'
    }
      
  }
  if(plot){
    gg <- ggraph(ink.layout)
    
    if(geom_edges == 'fan'){
      gg <- gg +
        geom_edge_fan(aes(colour = factor(ifelse(weight > 0, "positive", "negative")),
                                   width = abs(width)), alpha = .5, show.legend = FALSE)
    }else{
      gg <- gg + 
        geom_edge_diagonal(aes(edge_colour = factor(ifelse(weight > 0, "positive", "negative")),
                                        width = width), alpha = .5, show.legend = FALSE)
    }
    
    gg <- gg +
      geom_edge_loop0(aes(colour = factor(ifelse(weight > 0,"positive", "negative")),
                          width = abs(width)), alpha = .5, show.legend = FALSE) + 
      scale_edge_colour_manual(values = c("positive" = "#f8766d", "negative" = "#00bfc4"), 
                               guide = FALSE) +
      scale_edge_width(range = edge_width_range, guide = FALSE)
    
   if(!is.null(nodes_colour)){ # Use manual scale
     gg <- gg + geom_node_point(aes(size = spp_norm, colour = !!rlang::sym(colourname)), show.legend = legend) +
      scale_colour_manual(name = "Functional guilds", values = nodes_colour)
   }else{
     gg <- gg + geom_node_point(aes(size = spp_norm), colour = "darkgrey", show.legend = legend)
   }
   gg <- gg + scale_size_continuous(range = node_size_range, guide = FALSE) +
     geom_node_text(aes(label = label), size = labcex,
                              repel = check_overlap, show.legend = FALSE) +
      theme_void()
    
    return(gg)
  }
  else{
    return(ink.layout)
  }
}

plot_network.cooc_network <- function(nk, nodes_colour = NULL, base_layout = NULL, 
                                      plot = TRUE, s = 1, labcex = 4,
                                      edge_width_coeff = 20,
                                      check_overlap = TRUE, geom_edges, label,
                                      legend = FALSE, colourname = "specif_guild3"){
  # nk is a cooc_network object.
  # nodes_colour: dataframe with named rows, with names corresponding
  # to the nodes (any name that is not in the netwotk will be discarded)
  # plot: boolean (if true then the nk is plotted)
  # s and labcex: resp. size and label multiplication coefficients
  #     for abundance
  # edge_witdh_range: edge width range to use for plotting (controls edge width)
  # check_overlap: boolean (whether to optimise labels positions)
  # Returns ggraph plot or tbl_graph object.
  ink <- nk$network
  
  res <- plot_nk(nk = ink, nodes_colour = nodes_colour, 
                 base_layout = base_layout, 
                 plot = plot, s = s, labcex = labcex,
                 edge_width_coeff = edge_width_coeff,
                 check_overlap = check_overlap,
                 geom_edges = geom_edges, label = label,
                 legend = legend,
                 colourname = colourname)
  return(res)
}

plot_network.igraph <- function(nk, nodes_colour = NULL, base_layout = NULL, 
                                plot = TRUE, s = 1, labcex = 4,
                                edge_width_coeff = 20,
                                check_overlap = TRUE,
                                geom_edges, label,
                                legend = FALSE, colourname = "specif_guild3"){
  # nk is a cooc_network object.
  # nodes_colour: dataframe with named rows, with names corresponding
  # to the nodes (any name that is not in the netwotk will be discarded)
  # plot: boolean (if true then the nk is plotted)
  # s and labcex: resp. size and label multiplication coefficients
  #     for abundance
  # edge_witdh_range: edge width range to use for plotting (controls edge width)
  # check_overlap: boolean (whether to optimise labels positions)
  # Returns ggraph plot or tbl_graph object.
  
  if(is(nk, 'tbl_graph')){
    ink <- nk
  }else{
    ink <- as_tbl_graph(nk)
  }
  
  res <- plot_nk(nk = ink, nodes_colour = nodes_colour, 
                base_layout = base_layout, 
                plot = plot, s = s, labcex = labcex,
                edge_width_coeff = edge_width_coeff,
                check_overlap = check_overlap,
                geom_edges = geom_edges, 
                label = label, legend = legend,
                colourname = colourname)
  return(res)
}

spp_degree <- function(networks){
  # Returns a df of the degree distributions of several networks.
  # networks is a named list.
  Degree <- degree(networks[[1]])
  Dataset <- rep(names(networks)[1], length(Degree))
  for(i in 2:length(networks)){
    deg <- degree(networks[[i]])
    Degree <- c(Degree, deg)
    Dataset <- c(Dataset, rep(names(networks)[i], length(deg)))
  }
  Species <- names(Degree)
  df <- data.frame(Dataset, Degree, Species)
  return(df)
}

adjmatrix_from_count <- function(df){
  # This function creates an adjacency matrix from count data.
  # df is a count data table (only counts and not sites / times IDs.)
  countmat <- as.matrix(df)
  countmat[countmat > 0] <- 1
  
  adjmat <- t(countmat) %*% countmat
  return(adjmat)
}

incimat_from_count <- function(df, traits, groups = c("Carnivore", "Herbivore")){
  # This function creates an incidence matrix from count data (bipartite nk).
  # df is a count data table (only counts and not sites / times IDs.)
  # traits is a df with a column spp_mysites" and a column "guild_realm".
  # groups is a vector with the names of the 2 groups to keep (1st will be rows, 2nd cols)
  
  adjmat <- adjmatrix_from_count(df)
  traits.sub <- traits %>% filter(spp_mysites %in% colnames(adjmat) & 
                                    guild_realm %in% groups)
  incimat <- adjmat[as.character(traits.sub$spp_mysites[traits.sub$guild_realm == groups[1]]),]
  incimat <- incimat[,as.character(traits.sub$spp_mysites[traits.sub$guild_realm == groups[2]])]
  return(incimat)
}

plot_mat <- function(mat, reorder = TRUE, label = FALSE, labsize = 4, 
                     get_lower = FALSE, pal = 'Blues', legend = TRUE){
  # mat is a dist object as returned by distPairwise.
  # reorder: should the matriw be reordered by similarity?
  # label: whether to display labels.
  # labsize: label sizes
  # get_lower: whether to display lower triangular patrix only.
  # label: whether to display labels
  
  lim = c(-0.01,1.01)
  
  
  if(length(mat[is.na(mat)]) != 0){
    message("Warning: some distances are NAs and will be printed as -1.")
    mat[is.na(mat)] <- -1
  }
  
  if(reorder){
    hc <- hclust(mat, "centroid")
    ordered.hc <- reorder(hc, mat)
    mat2 <- as.matrix(mat)[ordered.hc$order, rev(ordered.hc$order)]
  }
  else{
    mat2 <- as.matrix(mat)
  }
  
  mid = .5
  
  if(get_lower){
    lower_tri <- get_lower_tri(mat2)
    melted_mat <- melt(lower_tri, na.rm = TRUE)
    
    melted_mat <- melted_mat %>% group_by(Var1) %>% mutate(n_Var2_for_Var1 = n())
    gg <- ggplot(data = melted_mat, aes(x = Var2, y = reorder(Var1, desc(n_Var2_for_Var1)),
                                        fill = value))
   
  }else{
    melted_mat <- melt(mat2, na.rm = TRUE)
    gg <- ggplot(data = melted_mat, aes(x = Var1, y = Var2, fill = value))
  }
  
  ncol <- 50
  cl <- colorRampPalette(brewer.pal(5, pal))(ncol)
  
  gg <- gg +
    geom_tile(color = "white") +
    scale_fill_gradientn(colours = cl, name = "Dissimilarity", limits = c(-0.01,1.01), breaks = c(0, .5, 1)) +
    theme_minimal()+ 
    coord_fixed() + theme(axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  if(label){
    gg <- gg + geom_text(aes(Var2, Var1, label = round(value,2)), size = labsize,
                         col = ifelse(melted_mat$value > .4, 'white', 'black'))
  }
  if(!legend){
    gg <- gg + theme(legend.position="none")
    
  }
  
  return(gg)
}


plot_mat2 <- function(mat, pal = 'Blues', title){
  # mat is a dist object as returned by distPairwise.
  # pal is the RColorBrewer palette to use.
  # title: title graph (optional)
  
  if(length(mat[is.na(mat)]) != 0){
    message("Warning: some distances are NAs.")
  }
  
  ncol <- 50
  cl <- colorRampPalette(brewer.pal(5, pal))(ncol)
  
  if(!missing(title)){
    heatmap(as.matrix(mat), symm = TRUE, col = cl, main = title,
            margins = c(0, 0),
            zlim = c(0,1))
  }else{
    heatmap(as.matrix(mat), symm = TRUE, col = cl,
            margins = c(0, 0),
            zlim = c(0,1))
  }
}

graph_all_nodes <- function(n, traits, by = "spp_mysites"){
  # ---- INTENDED FOR PRIVATE USE ---
  # /!\ will issue "Duplicate vertex names" error if one column in traits is not the same for all
  # (eg keeping species while nodes are indeed guilds)
  # Create a synthetic graph from a species list.
  # n is a nodes list.
  # traits is an (optional) dataframe of traits with one column identifying nodes.
  #   all nodes in traits absent from n will be discarded from the result.
  # by is a named vector specifying the names used to merge spp with traits
  #   (eg c("spp" = "spp_mysites")).
  # returns a graph with a single arbitrary edge and all nodes from spp.
  
  # Arbitrary edge
  from <- n[1]
  to <- n[2]
  edges <- data.frame(from, to)
  
  nodes <- as.data.frame(n) # A dtataframe with one column named n
  
  if(!(missing(traits))){ # if traits provided
    nodes <- left_join(nodes, traits, by = c("n" = by)) %>% distinct() # Add all columns in traits
  }

  g <- graph_from_data_frame(edges, vertices = nodes)
  res <- as_tbl_graph(x = g)
  
  return(res)
}


get_nodes_list <- function(X){
  UseMethod("get_nodes_list")
}

get_nodes_list.list <- function(X){
  # X: a list of graph objects.
  # Supported are igraphs, PLNnetwork or ccoc_network lists. 
  # Returns a character vector with unique nodes names.
  
  if(is(X[[1]], "igraph")){ # graphes are igraph objects
    nodes_table <- lapply(X, function(x) V(x)$name)
    nodes <- unique(unlist(nodes_table))
  }
  else if(is(X[[1]], "PLNfit")){ # graphes are PLNnetwork objects
    nodes_table <- lapply(X, function(x) V(plot(x, plot = FALSE))$name)
    nodes <- unique(unlist(nodes_table))
  }
  else if(is(X[[1]], "cooc_network")){ # graphes are PLNnetwork objects
    nodes_table <- lapply(X, function(x) V(x$network)$name)
    nodes <- unique(unlist(nodes_table))
  }
}

get_nodes_list.default <- function(X){
  # X: a graph object which is not a list.
  # Supported are igraphs, PLNnetwork or ccoc_network objects. 
  # Returns a character vector with unique nodes names.
  
  # --- Get the unique species list
  if(is(X, "igraph")){ # graphes are igraph objects
    nodes_table <- V(X)$name
    nodes <- unique(unlist(nodes_table))
  }
  else if(is(X, "PLNfit")){ # graphes are PLNnetwork objects
    nodes_table <- V(plot(X, plot = FALSE))$name
    nodes <- unique(unlist(nodes_table))
  }
  else if(is(X, "cooc_network")){ # graphes are PLNnetwork objects
    nodes_table <- V(X$network)$name
    nodes <- unique(unlist(nodes_table))
  }
  return(nodes)
}


common_layout <- function(nodes, traits, colname_for_nodes_in_traits = "spp_mysites",
                          layout_by = c('guild_realm', 'specif_guild3', 'weight')){
  # Define a common layout for species, given traits (optional)
  # nodes: a vector with unique nodes names (all the spp/groups present in the nks we want
  #     to plot.)
  # traits: optional dataframe with traits to use to produce the layout (basically all the traits you
  #   want in the final nodes attributes).
  # colname_for_nodes_in_traits is a string defining what is the colnames 
  #   identifying nodes in the (optional) traits table.
  # layout by: character vector with ordering colnames for the layout (represents the order around
  # the circle).
  # Returns a dataframe with x and y columns (coordinates), species' names and optional
  # factors of traits used for filtering.
  
  # --- Create a pseudo network to use a layout
  gr <- graph_all_nodes(n = nodes, traits = traits, 
                        by = colname_for_nodes_in_traits)
  
  # --- Create a layout for this graph
  if(!missing(traits)){ # if traits were provided
    gr.ordered <- gr %>% activate(nodes) %>% dplyr::arrange(!!! rlang::syms(layout_by))
  }
  else{ # just alphabetical order
    gr.ordered <- gr %>% activate(nodes) %>% dplyr::arrange(name)
  }
  
  layout <- create_layout(gr.ordered, layout = 'circle') %>% 
    dplyr::select(-c(.ggraph.orig_index, circular, .ggraph.index))
  
  return(layout)
}

get_my_layout <- function(my_network, base_layout){
  # --- Join layout to the nodes atttribute
  # Join
  my_network_tbl <- as_tbl_graph(my_network) %>% activate(nodes) %>%
    left_join(base_layout, by = c("name" = "name"))
  
  return(my_network_tbl)
}

prepare_list_for_econetwork <- function(l){
  # Prepare the list l for analysis with the econetwork package, namely:
  #   convert the networks to igraph objects
  #   convert links to absolute weight
  # l is a (best: named) list of PLNnetwork objects 
  if(is(l[[1]], "cooc_network")){
    l2 <- lapply(l, function(x) x$network)
  }else if(is(l[[1]], "PLNnetworkfit")){
    l2 <- lapply(l, function(x) plot(x, plot = FALSE))
  }else{
    l2 <- l
  }
  
  l2 <- lapply(l2, function(x) as_tbl_graph(x) %>% activate(edges) %>% mutate(weight = abs(weight)))
  
  return(l2)
}

prepare_abtable_for_econetwork <- function(l){
  # Prepare the list l to be an abundance table for econetwork.
  # l is a named list of networks that have been prepared for econetwork.
  # Returns a matrix n(spp) * n(sites).
  
  # --- Get specific abundances in each site
  res <- lapply(l, function(x){r <- x %>% activate(nodes) %>% as.data.frame() %>% 
    dplyr::select(name, abundance)
    re <- r$abundance
    names(re) <- r$name
    return(re)})
  Abundance <- unlist(res)
  
  # --- Get site and species ID
  spl <- strsplit(names(Abundance), "[.]")
  sites <- sapply(spl, function(x) x[1]) # corresponds to the names of the list
  spp <- sapply(spl, function(x) x[2]) # corresponds to the species' names
  
  # --- Coerce to dataframe
  # Group all columns in a df
  abundances <- as.data.frame(Abundance)
  abundances <- abundances %>% mutate(Site = sites,
                                      Species = spp)
  # wide form for species
  abundances.wide <- abundances %>% pivot_wider(names_from = Site, 
                                         values_from = Abundance)
  abundances.wide[is.na(abundances.wide)] <- 0
  abundances.wide <- abundances.wide %>% column_to_rownames('Species')
  
  countmat <- as.matrix(abundances.wide)
  return(countmat)
}

prepare_data_for_econetwork <- function(l){
  # Prepare the list l for analysis with the econetwork package, namely:
  #     Extracts networks only and absolute weights
  #     Prepares the abundance table from l.
  # l is a list of cooc_network objects.
  # Returns a list with the networks and the abundance table.
  
  l.prep <- prepare_list_for_econetwork(l)
  abtable <- prepare_abtable_for_econetwork(l.prep)
  
  return(list("networks" = l.prep,
              "abtable" = abtable))
}

group_nodes <- function(graph, group = "trophic_level", simplify = FALSE){
  # Get a graph summarising the relations as if all nodes belonging to one group 
  #   were a single group.
  # graph: a tbl_graph, igraph a cooc_network object.
  # group: name of the grouping variable (assumed to be present in the nodes of graph)
  # Simplify: return a multigraph?
  # Returns either a tbl_graph or e cooc_network, depending on input:
  #     Abundance of new nodes is the sum of the old nodes abundances.
  #     n_spp of the new nodes is the original unit richness in each node.
  #     Weight of new edges btw the same original nodes: sum
  #     Weight of new edges btw different original nodes: remain the same
  # Returns a multigraph (except if simplify is TRUE.)
  
  if(is(graph, 'cooc_network')){
    graph.g <- graph$network
  }else if(is(graph, 'igraph')){
    if(is(graph, 'tbl_graph')){
      graph.g <- graph
    }else{
      graph.g <- as_tbl_graph(graph)
    }
  }else{
    stop("Type not recognised.")
  }
  
  # --- Get simple nodes and edges df
  nodes <- graph.g %>% activate(nodes) %>% as.data.frame() %>% rownames_to_column("ID") %>%
    mutate(ID = as.integer(ID))
  edges <- graph.g %>% activate(edges) %>% as.data.frame()
  
  # Add the grouping variable to edges
  nodes.sub <- nodes %>% dplyr::select(ID, !! rlang::sym(group))
  edges <- edges %>% left_join(nodes.sub, by = c("from" =  "ID"))
  edges <- edges %>% rename(group_from = !! rlang::sym(group))
  
  edges <- edges %>% left_join(nodes.sub, by = c("to" =  "ID"))
  edges <- edges %>% rename(group_to = !! rlang::sym(group))
  
  # Summarise
  edges_new <- edges %>% group_by(from, to) %>% 
    summarise(n_edges = n(), # n edges btw the same spp
              weight = sum(weight), # weight of each node
              group_from = group_from,
              group_to = group_to)
  # Reorder/rename edge attributes for graph_from_df
  edges_new <- edges_new %>% rename(from_orig = from) %>% rename(to_orig = to) %>%
    rename(from = group_from) %>% rename(to = group_to) %>%
    dplyr::select(from, to, weight, n_edges, everything())
  
  nodes_new <- nodes %>% group_by(!! rlang::sym(group)) %>% 
    summarise(tot_abundance = sum(abundance), n_species = n()) %>% 
    mutate(avg_abundance = tot_abundance/n_species) %>%
    rename(name =  !! rlang::sym(group))
  
  igr <- graph_from_data_frame(d = edges_new, directed = FALSE, vertices = nodes_new)
  
  if(simplify){
    igr <- igraph::simplify(igr, remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb = "sum")
  }
  
  gr <- as_tbl_graph(igr)
  if(is(graph, 'cooc_network')){
    res <- graph
    res$network <- gr
  }else{
    res <- gr
  }
  
  return(res)
}

get_common_edges <- function(l, n = 'all'){
  UseMethod("get_common_edges")
} 

get_common_edges.list <- function(l, n = 'all'){
  if(is(l[[1]], "cooc_network")){
    l2 <- lapply(l, FUN = function(x) x$network)
  }else if(is(l[[1]], "igraph")){
    l2 <- l
  }
  res <- get_common_edges_general(l2, n = n)
  return(res)
}


get_common_edges_general <- function(l, n = 'all'){
  # --- INTENDED FOR PRIVATE USE ---
  # l is a list of igraph objects or tbl_graphs.
  # n is the minimal occurrence nb we want for this edge.
  # Returns a list of vertices shared by n graphs (defaults to all graphes).
  
  # Binarise graph
  l2 <- lapply(l, FUN = function(x) x %>% as_tbl_graph() %>% 
                 activate(edges) %>% mutate(weight = 1))
  
  # Select only the useful graph attributes
  l2 <- lapply(l2, FUN = function(x) {x %>% activate(edges) %>% 
      dplyr::select(weight, from, to) %>%
      activate(nodes) %>% dplyr::select(name)})
  
  # Compute graphs union
  union <- do.call(igraph::union,l2)
  
  # Compute edge count
  edg <- union %>% as_tbl_graph() %>% activate(edges) %>% as_tibble()
  edg$count = rowSums(edg[,3:ncol(edg)], na.rm = TRUE) 
  edg <- edg %>% select(from, to, count)
  
  # Merge count back in
  union <- union %>% as_tbl_graph() %>% activate(edges) %>% 
    left_join(edg, by = c("from", "to"))
  
  # Add edges names to vertices
  vertices <- as_tbl_graph(union) %>% activate(nodes) %>% as.data.frame() %>% 
    rownames_to_column("ID") %>%
    mutate(ID = as.integer(ID)) %>% dplyr::select(ID, name)
  
  res <- as_tbl_graph(union) %>% activate(edges) %>% left_join(vertices, by = c("from" =  "ID"))
  res <- res %>% rename(name_from = name)
  
  res <- res %>% left_join(vertices, by = c("to" =  "ID"))
  res <- res %>% rename(name_to = name)
  
  res.df <- res %>% activate(edges) %>% as.data.frame() %>% 
    dplyr::select(name_from, name_to, count) %>%
    mutate(edge_name = paste(name_from, name_to, sep = "--")) %>%
    dplyr::select(edge_name, everything()) %>% arrange(desc(count))
  
  if(n != 'all'){
    res.df <- res.df %>% dplyr::filter(count >= n)
  }
  
  return(res.df)
}

cormatrix <- function(df, thr, method = "spearman"){
  # Compute the complete correlation matrix from the expanded dataframe.
  # Thr is a threshold undre which correlations should be removed 
  # (if not present, nothing is deleted.)
  # method defaults to spearman but others czn be supplied.
  df2 <- as.data.frame(df)
  df2 <- df[, ! names(df) %in% c("time_ID","site_ID")]
  
  res <- cor(df2, method = method)
  if(!missing(thr)){
    res[abs(res) < thr] <- 0
  }
  res[is.na(res)] <- 0
  return(res)
}

# Get upper triangle of the correlation matrix
# From http://www.sthda.com/english/wi  ki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat, diag = TRUE)]<- NA
  return(cormat)
}

# Get lower triangle of the correlation matrix
# From http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat, diag = TRUE)]<- NA
  return(cormat)
}

plot_cormat <- function(cormat, lim = c(-1,1)){
  # lim sets the limits of the color scale
  upper_tri <- get_upper_tri(cormat)
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  if (lim == c(-1,1)){
    mid = 0
  }
  else{
    mid = (lim[2] - lim[1])/2
  }
  gg <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = mid, limit = lim, space = "Lab") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 7, hjust = 1))+
    coord_fixed() + theme(legend.position="none",
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  return(gg)
}

cormat_to_graph <- function(cormat, thr){
  # Transform a correlation matrix into a igraph object.
  # thr is a threshold under which the correlations (edges) should be removed 
  # (defaults to none)
  c2 <- c
  if(!missing(thr)){ # filter the matrix
    c2[abs(c)< thr] <- 0
  }
  
  gr <- graph_from_adjacency_matrix(c2, weight = TRUE, mode = "undirected",
                                    diag = FALSE)
  return(gr)
}

reorder_network_list <- function(network_list, flevels = c("cam", "clust", "site"),
                                 arrange = c("time", "scale", "n")){
  # Reorders a list of networks by time, then spatial scale
  # from finer to coarser.
  # Eg 2d_cam, 2d_clust1, 2d_clust2, 5d_cam...
  # network_list is a list of cooc_network objects.
  # flevels is an ordered (asc) factor list th give the spatial order.
  
  # Extract space and time from original list
  space <- sapply(network_list, function(x) x$space_scale)
  time <- sapply(network_list, function(x) as.numeric(x$timestep, "days"))
  
  # Format tables
  space <- matrix(unlist(space), ncol = 2, byrow = TRUE)
  rownames(space) <- names(time)
  colnames(space) <- c("n", "scale")
  space <- as.data.frame(space) %>% rownames_to_column("ID")
  time <- data.frame(time) %>% rownames_to_column("ID")
  
  # Merge time and space into df
  res <- inner_join(time, space, by = 'ID')
  
  # Order factor scale
  res$scale <- factor(res$scale, levels = flevels)
  
  # Convert n (factor) to numeric
  res$n <- as.numeric(levels(res$n))[res$n]
  
  # Get networks ordered list                 
  res <- res %>% arrange(!!! rlang::syms(arrange)) # desc n coz the less groups, 
  # the coarser the scale
  
  ordered_names <- res$ID
  lordered <- network_list[ordered_names]
  
  return(lordered)
}

add_node_attribute_to_list <- function(l, attribute){
  # Adds attributes to the nodes of the graphs in grlist by merging.
  # l is a list of cooc_networks or igraph or tbl_graph.
  # attribute is a df with n column and named rows that must match nodes names.
  # returns the list with added attributes.
  
  attr <- attribute %>% rownames_to_column("name")
  
  if(is(l[[1]], "cooc_network")){
    l.nks <- lapply(l, function(x) x$network %>% activate(nodes) %>%
                      left_join(attr, by = "name"))
    l.res <- lapply(seq_along(l), function(i) set_network(l[[i]],l.nks[[i]]))
    names(l.res) <- names(l)
    
  }else if(is(l[[1]], 'tbl_graph')){
    l.res <- lapply(l, function(x) x %>% activate(nodes) %>%
                      left_join(attr, by = "name"))
  }else if(is(l[[1]], 'igraph')){
    l.res <- lapply(l, function(x) as_tbl_graph(x) %>% activate(nodes) %>%
                      left_join(attr, by = "name"))
  }else{
    stop(paste('The type of l', paste(class(l[[1]]), sep = ","), 'is not recognised.'))
  }
  
  return(l.res)
}

my_graphs_union <- function(grlist, add = "sum"){
  # Performs the union of my graphs objects, wisely for attributes.
  # grlist is a list of igraph or cooc_network objects.
  # add: whether to sum ar average attributes
  
  # Returns a tbl_graph object.
  #     * New abundance: sum of abundances/average abundances.
  #     * New weight: sum of weights.
  #     Discards all other columns than name, abundance, from, to, weight.
  
  if(is(grlist[[1]], 'igraph')){
    l <- grlist
  }else if(is(grlist[[1]], 'cooc_network')){
    l <- lapply(grlist, function(x) as.igraph(x$network))
  }else{
    stop("grlist type is not recognised.")
  }
  
  # Discard useless attributes
  if(TRUE %in% grepl("avg_abundance",  names(vertex_attr(l[[1]])))){ # graph was merged
    l.filtered <- lapply(l, function(x){as_tbl_graph(x) %>% activate(nodes) %>%
        dplyr::select(name, avg_abundance)})
  }else{ # graph was not merged
    l.filtered <- lapply(l, function(x){as_tbl_graph(x) %>% activate(nodes) %>%
        dplyr::select(name, abundance)})
  }
  
  l.filtered <- lapply(l, function(x){as_tbl_graph(x) %>% activate(edges) %>%
      dplyr::select(from, to, weight)})
  
  # Graph union
  m <- do.call(igraph::union, l.filtered)
  
  m_tbl <- as_tbl_graph(m)
  
  if(TRUE %in% grepl("abundance",  names(vertex_attr(m_tbl)))){
    if(TRUE %in% grepl("avg_abundance",  names(vertex_attr(m_tbl)))){
      message("avg abundances detected in colnames, merging them.")
      # Subset df with those columns and nodes ID
      df <- m_tbl %>% activate(nodes) %>% 
        dplyr::select(name, contains("avg_abundance")) %>% as.data.frame()
    }else{
      message("abundances detected in colnames, merging them.")
      # Subset df with those columns and nodes ID
      df <- m_tbl %>% activate(nodes) %>% 
        dplyr::select(name, contains("abundance")) %>% as.data.frame()
    }
    
    if(add == "sum"){
      df$abundance <- apply(df[,-1], 1, sum, na.rm = TRUE) # Sum of abundance columns
    }else if(add == "mean"){
      df$abundance <- apply(df[,-1], 1, mean, na.rm = TRUE) # Mean of abundance columns
    }
    
    df <- df %>% dplyr::select(name, abundance)  # Select columns to be merged
    
    m_tbl <- m_tbl %>% activate(nodes) %>% select(!contains("abundance")) %>% 
      left_join(df, by = 'name')
  } 
  if(TRUE %in% grepl("weight",  names(edge_attr(m_tbl)))){
    message("weights detected in colnames, merging them.")
    df <- m_tbl %>% activate(edges) %>% 
      dplyr::select(from, to, contains("weight")) %>% as.data.frame()
    
    if(add == "sum"){
      df$weight <- apply(df[,-c(1,2)], 1, sum, na.rm = TRUE) # Sum of weight columns
    }else if(add == "mean"){
      df$weight <- apply(df[,-c(1,2)], 1, mean, na.rm = TRUE) # Mean of weight columns
    }
    
    df <- df %>% dplyr::select(from, to, weight)  # Select columns to be merged
    
    m_tbl <- m_tbl %>% activate(edges) %>% select(from, to) %>% 
      left_join(df, by = c('from', 'to'))
  } 
  return(m_tbl)
}