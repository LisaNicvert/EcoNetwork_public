########################################################
# Functions used for spatial clustering
########################################################

library(sf)
library(spdep)

library(ade4)
library(adespatial)
library(adegraphics)

library(dplyr)


coords_from_sfc <- function(sfcobj){
  # Returns a 2 column matrix x and y with coordinates extracted from the sfc object in entry.
  t <- as(sfcobj, "Spatial")
  mat.coord <- as.matrix(coordinates(t))
  colnames(mat.coord) <- c("x","y")
  return(mat.coord)
}

neighbours <- function(sfcobj, method = "dnearneigh", d1=0,
                       maxdist=3){
  # Find the neighbours of the points.
  # 2 methods are implemented: 
  #   dnearneigh is based on distance criteria (min and max)
  #   gabrielneigh is based on Gabriel graph.
  # Returns a nb object.
  if(method == "dnearneigh"){
    res <- dnearneigh(sfcobj, d1=d1, d2=maxdist)
  }
  else if(method == "gabrielneigh"){
    coords <- coords_from_sfc(sfcobj) # First we need the coordinated in matrix form
    res <- graph2nb(gabrielneigh(coords), sym = TRUE)
  }
  return(res)
}

SWM <- function(nbobj, sfcobj, weights = 'simple'){
  # Computes a spatial weighting matrix for the neighbours defined.
  # nbobj is the nb object we are interseted in
  # sfcobj is the corresponding sfc object or which the nd was calculated.
  # weights is a string that indicates which standardisation should be used:
  #   'simple' for row sums, 'dist' for distance weighting.
  # Returns a listw object.
  
  dist.all <- st_distance(sfcobj) # complete distance matrix

  dist <- nbdists(nbobj, st_geometry(sfcobj)) # compute distances between neighbours
  
  if(weights == 'dist'){
    w <- lapply(dist, function(x) 1 - x/as.numeric(max(dist.all)))
    res <- nb2listw(nbobj, glist = w)
  }
  else if(weights == 'simple'){
    res <- nb2listw(nbobj)
  }
  
  
  return(res)
}

prepare_for_MVA <- function(sfcobj, variables){
  # Prepare a sfc object for multivariate analysis.
  # sfcobj is the sfcobject to be prepared.
  # variables are a vector with columns names to extract for the multivariate analysis.
  for_mva <- sfcobj %>% dplyr::select(all_of(variables))
  
  # Prepare data for mv analysis
  rownames(for_mva) <- sfcobj$site_ID
  res <- st_drop_geometry(for_mva)  %>% droplevels
  
  return(res)
}

spatial_clustering <- function(dist, swm, coords){
  # Performs spatial clustering from
  #   dist: distance matrix in input (a dist-class dissimilarity matrix)
  #   swm: a listw weighted neighbouring list
  #   coords: the point coordinates (a n*2 matrix)
  
  neighbours2 <- listw2sn(swm)[,1:2]
  clust <- constr.hclust(dist,
                         coords = coords, 
                         links = neighbours2)
  return(clust)
}

spatial_clustering_from_sfc <- function(sfcobj, neighbour.method = "dnearneigh", d1=0,
                                        maxdist=3, variables, return.mva = TRUE,
                                        weights = 'simple', ...){
  # Performs spatial clustering on a sfc object.
  # sfcobj is the sfcobject for which we want to cluster the sites
  # neighbour.method is the method to use for finding neighbours (one of 
  #   dnearneigh and gabrielneigh)
  # d1 and maxdist are distance arguments
  # variables are the variables to include in the MVA
  # .. additional arguments for MVA (scannf and nf)
  
  
  # Neighbours
  nb <- neighbours(sfcobj, method = neighbour.method, d1=d1,
                   maxdist=maxdist)
  # SWM
  swm.listw <- SWM(nb, sfcobj) #, weights = weights)
  
  # Multivariate analysis
  df <- prepare_for_MVA(sfcobj, variables = variables) 
  hillsmith <- dudi.hillsmith(df, ...)
  dist <- dist.dudi(hillsmith)
  
  # Spatial clustering
  sfc.coords <- coords_from_sfc(sfcobj)
  clust <- spatial_clustering(dist = dist, 
                            swm = swm.listw, 
                            coords = sfc.coords)
  if(return.mva){
    res <- list(clust = clust,
                mva = hillsmith)
  }
  else{
    res <- clust
  }
  return(res)
}

plot_clust <- function(hclustobj, k, labels = TRUE){
  # Wrapper to plot spatial clustering nicely (with axes and points labels).
  # hclustobj is the object to plot.
  # k is the nb of clusters to show.
  # labels: if true, write site_ID on top of graphic.
  
  plot(hclustobj, links = TRUE, k = k, xlab = 'Longitude', ylab = 'Latitude')
  if(labels){
    text(hclustobj$coords[,2] ~ hclustobj$coords[,1], labels=hclustobj$labels, cex=0.7)
  }
}

cut_and_rename <- function(hclustobj, k, prefix){
  # Wrapper for cuttree functionn, that renames the clusters too in order to
  #   use them safely in dataframes (coz df doesn't like colnames beginning with figures).
  # hclustobj is the object to cut into clusters.
  # k is the desired clusters nb.
  # Prefix is an optional prefix to add before clusters (must begin with character).
  # Returns a named list with clusters and site_ID.
  
  r <- cutree(hclustobj, k=k)
  if(missing(prefix)){
    res <- paste0("C",r)
  }else{
    res <- paste0(prefix,r)
  }
  
  names(res) <- names(r)
  return(res)
}