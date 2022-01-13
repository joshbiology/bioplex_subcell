## ---------------------------
##
## Script name: cache_deepwalk.R
##
## Purpose of script: Assemble interactomes, and feed them through
##    DeepWalk.
##
## Author: Joshua Pan
##
## Date Created: 2022-01-09
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------



library(ProjectTemplate); load.project()
library(ggraph)
library(igraph)
library(tidygraph)
source('./munge/munge_interactomes.R')
source("./munge/munge_subcell.R")


#Step 1: Assemble networks
#Karate
karate <- make_graph("Zachary")
karate_df <- as_tibble(as_edgelist(karate))

#Play islands
islands <- play_islands(5, 10, 0.8, 3) %>% 
  mutate(community = as.factor(group_infomap()))

islands_df <- islands %>% activate(edges) %>% as_tibble()

#Simple PPI #####
#Autophagy
autophagy

#DUB
dub

#Gygi####
bioplex3

hct

bioplex2

rep_net

#Assemble into list
networks <- list(Karate = karate_df,
     Islands = islands_df,
     Autophagy = autophagy,
     DUB = dub,
     Bioplex3 = bioplex3,
     HCT = hct,
     Bioplex2 = bioplex2,
     Rep_Net = rep_net)

#Filter only to network edges that connect two genes with Gringas labels
networks_filtered <- c(networks[1:2], map(networks[3:8], ~filter(., Gene.x %in% bioid_loc$entrezgene & Gene.y %in% bioid_loc$entrezgene)))

map(names(networks_filtered), function(x) {
  pth <- file.path(".", "output", "edgelist", paste(x, "edgelist", sep = "."))
  networks_filtered[[x]] %>% select(1:2) %>%  na.omit() %>% write_tsv(path = pth, col_names = F)
})


# Filter by subcell -------------------------------------------------------
#When running yourself, please be careful of paths.
today <- Sys.Date()
out_path <- file.path(".", "data", "interim", "deepwalk_embeddings", today)
dir.create(out_path, recursive = T)

#Run bash script to 
cat("#!/bin/bash", map_chr(names(networks_filtered), function(x) {
  in_pth <- file.path(".", "output", "edgelist", paste(x, "edgelist", sep = "."))
  out_pth <- file.path(out_path, paste(x, "deepwalk", "txt", sep = "."))
  bash <- paste("deepwalk --format edgelist --representation-size 64 --walk-length 40 --number-walks 16 --workers 8 --input", in_pth, "--output", out_pth)
}), sep = "\n", file = file.path(".", "scripts", "deepwalk_batch.sh"))

#Make sure you have deepwalk installed!
#https://pypi.org/project/deepwalk/

#Import embeddings
deepwalks <- list.files("./data/interim/deepwalk_embeddings", pattern = "deepwalk", full.names = T) %>% 
  map(~read_delim(., delim = " ", skip = 1, col_names = F) %>%  rename(name = X1))  


names(deepwalks) <- list.files("./data/interim", pattern = "deepwalk", full.names = F) %>% word(1, sep = "\\.")

deepwalk_mat <- deepwalks %>% 
  map(~as.data.frame(.) %>% set_rownames(.$name) %>% select(-name) %>% as.matrix())

ProjectTemplate::cache("deepwalk_mat")




# Redo with FULL networks -------------------------------------------------
#For these embeddings, filter to network edges that connect at least one gene
#with Gringas labels.
#The idea is to minimize uninformative edges before embedding.

networks_filtered <- c(networks[1:2], map(networks[3:8], ~filter(., Gene.x %in% bioid_loc$entrezgene | Gene.y %in% bioid_loc$entrezgene)))

networks %>% map(igraph::graph_from_data_frame) %>% map_dbl(~igraph::V(.) %>% length)

networks_filtered %>% map(~igraph::graph_from_data_frame(.,directed = F)) %>% map_dbl(~igraph::V(.) %>% length)


map(names(networks_filtered), function(x) {
  pth <- file.path(".", "output", "edgelist", paste(x,"full", "edgelist", sep = "."))
  networks_filtered[[x]] %>% select(1:2) %>%  na.omit() %>% write_tsv(path = pth, col_names = F)
})


#Generate bash script 
cat("#!/bin/bash", map_chr(names(networks_filtered), function(x) {
  in_pth <- file.path(".", "output", "edgelist", paste(x, "full","edgelist", sep = "."))
  out_pth <- file.path(out_path, paste(x,"full", "deepwalk", "txt", sep = "."))
  bash <- paste("deepwalk --format edgelist --representation-size 128 --walk-length 40 --number-walks 16 --workers 8 --input", in_pth, "--output", out_pth)
}), sep = "\n", file = "./scripts/deepwalk_full_batch.sh")



# Cache -------------------------------------------------------------------

#Hard coded to import the September 2021 run by Joshua Pan

frozen_path <- "./data/interim/deepwalk_embeddings/2021-09-16/"

#Import embeddings
deepwalks <- list.files(frozen_path, pattern = "full.deepwalk", full.names = T) %>% 
  map(~read_delim(., delim = " ", skip = 1, col_names = F) %>%  rename(name = X1))  


names(deepwalks) <- list.files(frozen_path, pattern = "full.deepwalk", full.names = F) %>% word(1, sep = "\\.")


deepwalk_mat_full <- deepwalks %>% 
  map(~as.data.frame(.) %>% set_rownames(.$name) %>% select(-name) %>% as.matrix())

ProjectTemplate::cache("deepwalk_mat_full")


