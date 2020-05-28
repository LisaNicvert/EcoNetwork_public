source("functions/analysis.R")
source("functions/networks.R")

null_model <- function(counts, nrep = 20, crit = "StARS", cov, sppnames, 
                       use_cov = FALSE){
  nspp <- ncol(counts)
  
  a <- array(dim = c(nspp, nspp, nrep),
             dimnames = list(sppnames,
                             sppnames,
                             NULL))
  
  for(i in 1:nrep){
    message(paste0("------------------------------------ Repetition ", i, "/", nrep))
    
    # Permutation of original data
    M.permut <- apply(counts, 2, sample, size = nrow(counts), replace = FALSE)
    rownames(M.permut) <- rownames(counts)
    
    # Prepare for PLN
    prep.permut <- prepare_data(counts = M.permut, 
                                covariates = cov)
    
    # Infer nk
    if(use_cov){
      nks <- PLNnetwork(Abundance ~ 1 + cov, data = prep.permut)
    }else{
      nks <- PLNnetwork(Abundance ~ 1, data = prep.permut)
    }
    
    if(sum(nks$coefficient_path()$Coeff) != 0){ # if there are potential edges
      nk <- getBestModel(nks, crit = crit)
      
      if(nk$density != 0){ # if there are edges in the final graph
        # Adjacency matrix
        nk.adj <- as_adjacency_matrix(plot(nk, plot = FALSE), 
                                      type = "lower", attr = "weight",
                                      sparse = FALSE)
      }
    }else{ # Empty adjacency matrix
      nk.adj <- matrix(data = rep(0, nspp*nspp), ncol = nspp, nrow = nspp)
    }
    
    # Store in array
    a[,,i] <- nk.adj
  }
  return(a)
}



# Test null model (permute columns)

# --------- Read original data
dat <- readRDS("null_model/SAM_cooc_10d_cam.rds")
counts <- dat$data$Abundance
cov <- dat$data$cov_fac

ref.adj <- as_adjacency_matrix(dat$network,
                              type = "lower", attr = "weight",
                              sparse = FALSE)
sppnames <- colnames(ref.adj)

# --------- Parameters
crit <- "BIC"
nrep <- 10

# # Uncomment to run simulations
# a <- null_model(counts = counts, nrep = nrep, crit = crit, sppnames = sppnames, cov = cov,
#                 use_cov = FALSE)
# saveRDS(a, "null_model/nullmodel_simul.rds")

a <- readRDS("null_model/nullmodel_simul.rds")


# --------- Get proportion of positive edges in simulations
simul.signif.count <- apply(a, c(1,2), FUN = function(x){length(x[x != 0])/nrep})


# --------- Gather reference and null array in an array
dim.adj <- nrow(ref.adj)
test <- array(data = c(ref.adj, simul.signif.count), dim = c(dim.adj, dim.adj, 2),
              dimnames = list(colnames(ref.adj),
                              colnames(ref.adj), NULL))

# Proba for initially detected edges to appear due to chance
prob <- test[,,2][test[,,1] != 0]
prob

names <- which(test[,,1] != 0, arr.ind = TRUE)
names

edgenames <- c()
for(i in 1:nrow(names)){
  e <- c(colnames(test[,,1])[names[i,1]], colnames(test[,,1])[names[i,2]])
  edgenames <- rbind(edgenames, e)
}

edgenames <- cbind(edgenames, prob)
colnames(edgenames) <- c("from", "to", "Frequency")

edgenames <- as.data.frame(edgenames)
edgenames$Frequency <- as.numeric(as.character(edgenames$Frequency))

abd <- dat$network %>% activate(nodes) %>% as.data.frame()

# ----------------------
# Plot graph
# ----------------------
traits <- read.csv("data/traits/traits.csv")

traits.sub <- traits %>% 
  dplyr::select(spp_mysites, guild_realm, specif_guild3, weight)

labels <- traits %>% dplyr::select(spp_mysites, label)

# ------------------- Graph from adjmatrix
ref.gr <- graph_from_adjacency_matrix(ref.adj, mode = 'lower', weighted = TRUE)
ref.gr <- as_tbl_graph(ref.gr)

ref.gr <- ref.gr %>% activate(edges) %>% mutate(from_name = .N()$name[from]) %>%
  mutate(to_name = .N()$name[to]) %>%
  left_join(edgenames, by = c("from_name" = "to",
                              "to_name" = "from"))

ref.gr <- ref.gr %>% activate(nodes) %>% left_join(abd, by = "name") %>%
  left_join(labels, by = c("name" = "spp_mysites")) 

# ----------------- Plot
spp <- ref.gr %>% activate(nodes) %>% as.data.frame()
spp <- spp$name

clayout <- common_layout(nodes = spp, traits = traits.sub, 
                         colname_for_nodes_in_traits = 'spp_mysites')

# Filter to get species in the site
cols <- traits %>% filter(spp_mysites %in% spp)
cols.c <- get_guild_colours(cols)

gr.layout <- create_layout(ref.gr, layout = get_my_layout, base_layout = clayout)

pal <- colorRampPalette(rev(brewer.pal(11, "RdYlGn")))

ggraph(gr.layout) + 
  geom_edge_diagonal(aes(width = abs(weight),
                         col = Frequency), alpha = .7) +
  geom_node_point(aes(size = abundance), colour = "cornflowerblue", show.legend = FALSE) +
  geom_node_text(aes(label = label), repel = TRUE) +
  scale_edge_width(range = c(1, 5), guide = FALSE) +
  scale_size_continuous(range = c(5, 10)) +
  # scale_colour_manual(name = "Functional guilds", values = cols.c) +
  scale_edge_colour_gradientn(colours = pal(100), limits = c(0,.20), breaks = c(0, .10, .20)) + theme_void() + 
  theme(legend.text=element_text(size=10))
