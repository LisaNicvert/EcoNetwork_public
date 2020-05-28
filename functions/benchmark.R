########################################################
# Functions for model benchmark (synthetic data)
########################################################

library(igraph)
library(MASS)
library(PLNmodels)

generate_precmat <- function(n = 20, m = 3, u = .1, v = .3, 
                             gr, names){
  # Generate a precision matrix with the given parameters.
  # n: nodes count
  # gr (optional): the graph to use
  # m: density of network
  # v and u: parameters for resp. correlation and conditioning.
  # names (optional): colnames to give to the graph nodes
  
  if(missing(gr)){ # generate network
    gr <- sample_pa(n = n, m= m, 
                    directed=FALSE) # Network
  }
  
  # Create Omega with same sparsity pattern
  G <- as_adjacency_matrix(gr, sparse = FALSE)
  eig <- eigen(v*G, only.values = TRUE)$values
  Omega <- v*G + diag(abs(min(eig))+u,n,n)
  
  if(missing(names)){
    if(n <= 26){
      colnames(Omega) <- letters[1:n]
      rownames(Omega) <- letters[1:n]
    }
    else if(n > 26 & n <= 26*2){
      colnames(Omega) <- c(letters, LETTERS[1:(n-26)])
      rownames(Omega) <- colnames(Omega)
    }
    else{
      colnames(Omega) <- c(letters, LETTERS, seq(22*2, n))
      rownames(Omega) <- colnames(Omega)
    }
  }
  else{
    colnames(Omega) <- names
    rownames(Omega) <- names
  }
  return(Omega)
}

generate_obs <- function(nrow = 10, mean, Omega){
  # generate observations based on the given parameters.
  # nrow: number of observations to generate (replicates)
  # mean: the ordered mean vector (defaults to 0n)
  # Omega: precision matrix
  # Observations' colnames are Omega's colnames, rownames are 1:n.
  
  n <- nrow(Omega) # spp nb
  
  if(missing(mean)){ # mean defaults to zero
    mean <- rep(0, n)
  }
  
  Sigma <- solve(Omega)
  
  # Draw latent variables
  Z <- mvrnorm(n = nrow, mu = rep(0, n), Sigma = Sigma)
  colnames(Z) <- colnames(Omega)
  rownames(Z) <- seq(1, nrow)
  
  # Generate obs from Zi
  obs <- matrix(rpois(n = nrow*n, lambda = exp(Z + mean)), nrow = nrow)
  colnames(obs) <- colnames(Omega)
  rownames(obs) <- seq(1, nrow)
  
  return(list(lambda = exp(Z + mean),
              obs = obs))
}

sensi_speci <-function(n = 20, m = 3, u = .1, v=.3, gr, names,
                       mean, size = seq(10,100, by = 10), nrep = 10,
                       change = "size", new_graph = FALSE,
                       crit = "BIC"){
  # Compute FDR and FOR, and True/False Positives/Negatives evolution for different
  # sample sizes or means, depending on parameter "change".
  # n: number of nodes.
  
  # m: number of edges to attach to each step for Barabasi-Albert.
  # v and u: parameters for resp. correlation and conditioning.
  
  # mean: mean vector for latent Normal multivariate layer (dimension n if size
  #     changes; else, dimension n*whatever).
  
  # size: the different sample size to test (seq if size changes, 
  # else integer.).
  # Must be > 1.
  
  # nrep: repetitions number for each sample size (redraw Zi).
  
  # change: what changes between iterations. Possible values:
  #     mean, size.
  # new_graph: if TRUE, a new precision matrix is generated each time.
  
  # crit: the selection criterion to use with PLNnetwork ('BIC' or 'StARS').
  
  if(missing(mean)){
    if(change == "size"){
      mean <- rep(0,n) # initialises mean vector
    }
    else if(change == "mean"){ # mean as a matrix
      mean <- matrix(unlist(lapply(seq(0,5), FUN = rep, n)), ncol = n,
                  byrow = TRUE)
    }
  }
  # Initialisation of results vectors
  all.sensi <- vector()
  all.speci <- vector()
  
  if(change == "size"){
    jlim <- length(size) # set main loop size
    mean.now <- mean # mean cst
  }
  else if(change == "mean"){
    jlim = nrow(mean) # set main loop size
    nrow <- size # size cst
  }
  
  if(!new_graph){
    # --- Simulate a cov matrix
    if(missing(gr)){ # generate network
      Omega <- generate_precmat(n = n, m = m, u = u, v = v, 
                                gr = gr)
    }
    else{
      Omega <- generate_precmat(n = n, m = m, u = u, v = v)
    }
  }
  
  for(j in 1:jlim){ # loop for sample size / mean
    if(change == "size"){
      nrow <- size[j]
    }
    else if(change == "mean"){
      mean.now <- mean[j,]
    }
    
    # Initialisation for the same sample size / mean
    sensitivities <- vector()
    specificities <- vector()
    
    for(i in 1:nrep){
      rep <- i
      
      if(change == "size"){
        print(paste("----------------- nrep =", i, "| nrow =", nrow, "-----------------"))
      }
      else if(change == "mean"){
        print(paste("----------------- nrep =", i, "| mean vector =", paste0("(",paste(mean[j,],collapse = ","),")"), "-----------------"))
      }
      
      if(new_graph){
        # --- Simulate a cov matrix
        if(missing(gr)){ # generate network
          Omega <- generate_precmat(n = n, m = m, u = u, v = v, 
                                    gr = gr)
        }
        else{
          Omega <- generate_precmat(n = n, m = m, u = u, v = v)
        }
      }
      
      # Generate observations
      rand <- generate_obs(nrow = nrow, mean = mean.now, Omega = Omega)

      # Infer network
      cov <- rep(0, nrow) # initialise covariates
      # names(cov) <- rownames(rand$obs) # match rownames to disable warnings
      
      
      d.prep <- prepare_data(counts = rand$obs, covariates = cov)
      
      nks <- PLNnetwork(Abundance ~ 1, data = d.prep)
      nk <- getBestModel(nks, crit = crit)
      
      # Estimate for Omega
      Omega.est <- nk$model_par$Omega

      # Estimate for partial correlations
      Pcor.est <- Omega.est
      Pcor <- Omega 
      for(k in 1:nrow(Omega)){ # nrow(Omega) = nrow(Omega.est) so not important
        for (l in 1:ncol(Omega)){ # same for ncol
          # Real partial correlation
          Pcor[k,l] <- -Omega[k,l]/sqrt(Omega[k,k]*Omega[l,l])
          # Estimate partial correlation
          Pcor.est[k,l] <- -Omega.est[k,l]/sqrt(Omega.est[k,k]*Omega.est[l,l])
        }
      }
      
      # False omission rate
      # -- For Pcor
      TN <- length(Pcor.est[Pcor.est == 0 & Pcor == 0])
      FN <- length(Pcor.est[Pcor.est == 0 & Pcor != 0])
      FOR <- FN/(TN + FN)
      
      # False discovery rate
      # -- For Pcor
      FP <- length(Pcor.est[Pcor.est != 0 & Pcor == 0])
      TP <- length(Pcor.est[Pcor.est != 0 & Pcor != 0])
      FDR <- FP/(TP + FP)
      
      if(change == "size"){
        fac <- nrow
        r <- cbind(fac, rep, TP, FP, TN, FN, 
                     FDR, FOR)
      }
      else if(change == "mean"){ # mean of means vector used here
        fac <- mean(mean)
        r <- cbind(fac, rep, TP, FP, TN, FN, 
                     FDR, FOR)
      }
      # Update
      if(i == 1){
        res.rep <- r
      }
      else{
        res.rep <- rbind(res.rep, r)
      }
    }
    
    # Update
    if(j == 1){
      res <- res.rep
    }
    else{
      res <- rbind(res, res.rep)
    }
  }
  res <- as.data.frame(res)
  
  if(change == "size"){
    colnames(res)[colnames(res) == "fac"] <- "nrow"
  }
  else if(change == "mean"){
    colnames(res)[colnames(res) == "fac"] <- "mean"
  }
  
  return(res)
}

repeat_simul <- function(n, nrep, nrow, mean, Omega, crit){
  # Repeat the simulation nrep times with the same setting.
  for(i in 1:nrep){
    message(paste0("------------------------------ Repetition ", i,"/",nrep))
    # Generate observations
    rand <- generate_obs(nrow = nrow, mean = mean, Omega = Omega)
    nrowfinal <- nrow(rand$obs[apply(rand$obs, 1, sum) != 0,])
    
    # Infer network
    cov <- rep(0, nrow) # initialise covariates
    # names(cov) <- rownames(rand$obs) # match rownames to disable warnings
    
    d.prep <- prepare_data(counts = rand$obs, covariates = cov)
    
    nks <- PLNnetwork(Abundance ~ 1, data = d.prep)
    nk <- getBestModel(nks, crit = crit)
    
    # Estimate for Omega
    Omega.est <- nk$model_par$Omega
    
    
    # -- Get half
    Omega.est.tri <- Omega.est[lower.tri(Omega.est, diag = FALSE)]
    Omega.tri <- Omega[lower.tri(Omega, diag = FALSE)]
    
    TN <- length(Omega.est.tri[Omega.est.tri == 0 & Omega.tri == 0])
    FN <- length(Omega.est.tri[Omega.est.tri == 0 & Omega.tri != 0])
    FP <- length(Omega.est.tri[Omega.est.tri != 0 & Omega.tri == 0])
    TP <- length(Omega.est.tri[Omega.est.tri != 0 & Omega.tri != 0])
    
    # False negative rate
    FNR <- FN/(TP + FN)
    
    if(i == 1){
      res <- c(nrowfinal, TN, FN, FP, TP, FNR)
    }else{
      r <- c(nrowfinal, TN, FN, FP, TP, FNR)
      res <- rbind(res, r)
    }
  }
  res <- as.data.frame(res)
  colnames(res) <- c("nrow", "TN", "FN", "FP", "TP", "FNR")
  rownames(res) <- NULL
  return(res)
}