#comparing class probabilities trained on location vs. compartment


library(ProjectTemplate); load.project()
library(ggraph)
library(igraph)
library(tidygraph)
source('./munge/interaction/munge_interactomes.R')
source("./munge/gene/munge_genesets.R")

# Import labels -----------------------------------------------------------

source("./munge/gene/munge_subcell.R")

bioid_loc <- bioid_loc %>% 
  mutate(NMF_Rank = as.factor(NMF_Rank) %>% forcats::fct_explicit_na()) %>% 
  unique()

nmf_rank_to_name <- bioid_loc %>% 
  select(NMF_Rank, Location) %>% 
  unique() %>% 
  pull(Location)


names(nmf_rank_to_name) <- bioid_loc %>% 
  select(NMF_Rank, Location) %>% 
  unique() %>% 
  pull(NMF_Rank)

View(nmf_rank_to_name %>% enframe)

bioid_nmf_mixing <- read_excel("./data/raw/gringas/stdenis_rank_profile.xlsx", sheet = 2, .name_repair ="universal") %>% 
  mutate(entrezgene = convert_genes_mygeneinfo(gene, to = "entrezgene")) 


bioid_mixing_mat <- bioid_nmf_mixing %>% 
  group_by(entrezgene) %>% 
  filter(row_number() == 1, !is.na(entrezgene)) %>% 
  ungroup() %>% 
  select(-gene) %>% 
  column_to_rownames("entrezgene") %>% 
  as.matrix()

bioid_loc <- bioid_loc %>% filter(entrezgene %in% rownames(bioid_mixing_mat))


# Visualizations on mixing matrix ---------------------------------------------

all_genes <- bioid_loc %>% arrange(NMF_Rank) %>% pull(entrezgene) %>% na.omit()
mito_genes <- bioid_loc %>% filter(NMF_Rank == 4) %>% pull(entrezgene) %>% na.omit()


pheatmap::pheatmap(bioid_mixing_mat[intersect(all_genes, rownames(bioid_mixing_mat)),] %>% 
                     set_colnames(nmf_rank_to_name[as.character(1:20)]), cluster_cols = F, cluster_rows = F, show_rownames = F, angle_col = 45)

bioid_mixing_mat %>% rowSums() %>% hist()

#Localization correlation
cor(bioid_mixing_mat) %>% set_colnames(nmf_rank_to_name[as.character(1:20)]) %>% 
  set_rownames(nmf_rank_to_name[as.character(1:20)]) %>% pheatmap::pheatmap()


manual_mapping <- list(Nucleus = c("8", "2", "9",  "10"),
                       Cyto = c("16", "12", "5", "17"),
                       Mito = c("4", "18", "13"),
                       ER = c("3", "6", "15"),
                       Traffic = c("1", "11", "14", "7", "20"),
                       Misc = "19")


localization_to_compartment <- enframe(manual_mapping) %>% unnest() %>% pull("name")
names(localization_to_compartment) <-  enframe(manual_mapping) %>% unnest() %>% pull("value")


bioid_loc <- bioid_loc %>%
  mutate(Compartment = factor(localization_to_compartment[bioid_loc$NMF_Rank %>% as.character] %>% as.character(), levels = names(manual_mapping)))

bioid_nmf_mixing_compartment <- aaply(bioid_mixing_mat, 1, function(x) map_dbl(manual_mapping, function(y)sum(x[as.numeric(y)])))
  
# Assemble networks -------------------------------------------------------


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
networks <- list(Autophagy = autophagy,
                 DUB = dub,
                 Bioplex3 = bioplex3,
                 HCT = hct,
                 Bioplex2 = bioplex2,
                 Rep_Net = rep_net)

# Filter by subcell -------------------------------------------------------
networks_filtered <- c(map(networks, ~filter(., Gene.x %in% bioid_loc$entrezgene & Gene.y %in% bioid_loc$entrezgene)))



# Import embeddings -------------------------------------------------------
load("./cache/deepwalk_mat.RData")


#Split data

deepwalk_location_df <- map(deepwalk_mat[names(networks_filtered)], function(x) inner_join(bioid_loc %>% select(NMF_Rank, entrezgene),
                                                                                  as_tibble(x %>% scale(), rownames = "entrezgene")) %>% 
                     column_to_rownames("entrezgene"))

deepwalk_location_df_soft_labels <- map(deepwalk_location_df, function(x) bioid_mixing_mat[rownames(x),])

#Freeze for both location and compartment
train_indices <- map(deepwalk_location_df, ~caret::createDataPartition(.$NMF_Rank, p = .75, 
                                                              list = FALSE, 
                                                              times = 1))

deepwalk_location_split <- map2(deepwalk_location_df, train_indices, function(x, y) {list(train = x[y,], test = x[-y,])})

deepwalk_location_soft_labels_split <- map2(deepwalk_location_df_soft_labels, train_indices, function(x, y) {list(train = x[y,], test = x[-y,])})


# Repeat for compartment

deepwalk_compartment_df <- map(deepwalk_mat[names(networks_filtered)], function(x) inner_join(bioid_loc %>% select(Compartment, entrezgene),
                                                                                           as_tibble(x %>% scale(), rownames = "entrezgene")) %>% 
                              column_to_rownames("entrezgene"))

deepwalk_compartment_df_soft_labels <- map(deepwalk_compartment_df, function(x) bioid_mixing_mat[rownames(x),])

deepwalk_compartment_split <- map2(deepwalk_compartment_df, train_indices, function(x, y) {list(train = x[y,], test = x[-y,])})

deepwalk_compartment_soft_labels_split <- map2(deepwalk_compartment_df_soft_labels, train_indices, function(x, y) {list(train = x[y,], test = x[-y,])})


# XGBoost module ----------------------------------------------------------


library(doParallel)
cl <- makePSOCKcluster(6)
registerDoParallel(cl)



grid_default <- expand.grid(
  nrounds = 200,
  eta = c(0.05, 0.3),
  max_depth = c( 4, 6),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

train_control <- caret::trainControl(
  method = "cv",
  number = 3,
  verboseIter = TRUE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

#Classification by location
bioplex3_location_xgboost <- caret::train(
  x = deepwalk_location_split$Bioplex3$train[-1],
  y = deepwalk_location_split$Bioplex3$train[[1]],
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE,
  objective='multi:softprob',
  num_class = 20
)

#saveRDS(bioplex3_location_xgboost, "./bioplex3_location_xgboost.rds")
#bioplex3_location_xgboost <- readRDS("./bioplex3_location_xgboost.rds")

#Classification by Compartment
bioplex3_compartment_xgboost <- caret::train(
  x = deepwalk_split$Bioplex3$train[-1],
  y = deepwalk_split$Bioplex3$train[[1]],
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE,
  objective='multi:softprob',
  num_class = 6
)

#saveRDS(bioplex3_compartment_xgboost, "./bioplex3_compartment_xgboost.rds")
#bioplex3_compartment_xgboost <- readRDS("./bioplex3_compartment_xgboost.rds")


# Generate class prob ----------------------------------------------------------

location_class_probabilities <- predict(bioplex3_location_xgboost$finalModel, newdata = deepwalk_location_split$Bioplex3$test[-1] %>% as.matrix()) 
location_class_probabilities <-  matrix(location_class_probabilities, nrow = 20) %>% t()
rownames(location_class_probabilities) <- rownames(deepwalk_location_split$Bioplex3$test)


#Compute class probabilities

compartment_class_probabilities <- predict(bioplex3_compartment_xgboost$finalModel, newdata = deepwalk_compartment_split$Bioplex3$test[-1] %>% as.matrix()) 
compartment_class_probabilities <-  matrix(compartment_class_probabilities, nrow = 6) %>% t()
rownames(compartment_class_probabilities) <- rownames(deepwalk_compartment_split$Bioplex3$test)
colnames(compartment_class_probabilities) <- names(manual_mapping)

#True classes

location_true_class_prob <- t(bioid_mixing_mat) %*% diag(1/rowSums(bioid_mixing_mat)) %>% t() %>% 
  set_rownames(rownames(bioid_mixing_mat))

compartment_true_class_prob <- t(bioid_nmf_mixing_compartment) %*% diag(1/rowSums(bioid_nmf_mixing_compartment)) %>% t() %>% 
  set_rownames(rownames(bioid_nmf_mixing_compartment))


#Evaluate log-loss


compartment_log_loss <- map_dbl(deepwalk_compartment_split$Bioplex3$test %>% rownames, ~ModelMetrics::logLoss(compartment_class_probabilities[.,],
                                                                          compartment_true_class_prob[.,]))


location_log_loss <- map_dbl(deepwalk_location_split$Bioplex3$test %>% rownames, ~ModelMetrics::logLoss(location_class_probabilities[.,],
                                                                          location_true_class_prob[.,]))

log_loss_df <- tibble(entrezgene = deepwalk_location_split$Bioplex3$test %>% rownames,
                      symbol = convert_genes(entrezgene),
                      Location_Log_Loss = location_log_loss,
                      Compartment_Log_Loss = compartment_log_loss)

qplot(compartment_log_loss, location_log_loss)


# Generate plots ----------------------------------------------------------
focus_gene <- "54205"

embedding_report_plot <- function(focus_gene) {
  
  
  
  # Network plots ------------------------------------------------
  
  require("grid")
  require("ggplotify")
  require( 'gridExtra' )
  
  
  bioplex3_igraph <- graph_from_data_frame(networks_filtered$Bioplex3)
  bioplex3_secondary_neighbors <- ego(bioplex3_igraph, order = 2) %>% 
    set_names(V(bioplex3_igraph)$name)
  
  tmp_graph <- induced_subgraph(bioplex3_igraph, bioplex3_secondary_neighbors[[focus_gene]])
  
  tmp_graph_2 <- as_tbl_graph(tmp_graph) %>% 
    activate(nodes) %>% 
    left_join(bioid_loc, by = c("name" = "entrezgene")) %>% 
    mutate(Location = str_sub(Location, 1, 25))
  
  #Location plot
  g1 <- ggraph(tmp_graph_2, layout="focus", focus = 1) +
    geom_edge_fan(color = "gray", alpha = 0.5) +
    geom_node_point(aes(color =Location), size = 2) +
    geom_node_text(aes(label = symbol, color =Location), size = 3, nudge_y = 0.3)
  
  d1 <-  
    tableGrob(tmp_graph_2 %>% 
                activate(nodes) %>% 
                as_tibble() %>% 
                count(Location) %>% 
                arrange(-n) %>% 
                slice(1:6))
  
  #Compartment plot
  g2 <- ggraph(tmp_graph_2, layout="focus", focus = 1) +
    geom_edge_fan(color = "gray", alpha = 0.5) +
    geom_node_point(aes(color =Compartment), size = 2) +
    geom_node_text(aes(label = symbol, color =Compartment), size = 3, nudge_y = 0.3)
  
  d2 <- tableGrob(tmp_graph_2 %>% 
                    activate(nodes) %>% 
                    as_tibble() %>% 
                    count(Compartment) %>% 
                    arrange(-n))
  
  network_plot <- cowplot::plot_grid(d1, g1, d2 , g2, rel_widths = c(0.3, 0.7)) 
  
  
  # Class prob heatmaps plots ----------------------------------------------------------
  require(pheatmap)
  
  location_row_order <- manual_mapping %>% enframe("Compartment", "NMF_Rank") %>% unnest() %>% 
    mutate(Location = nmf_rank_to_name[NMF_Rank] %>% str_sub(1, 25))
  
  
  
  red_palette <- c("#FFFFFF", RColorBrewer::brewer.pal(n = 9, name =
                                                         "Reds"))[1:8]
  
  heat_1 <- rbind(location_true_class_prob[focus_gene,location_row_order$NMF_Rank %>% as.numeric()],
                  location_class_probabilities[focus_gene,location_row_order$NMF_Rank %>% as.numeric()]) %>% 
    t() %>% 
    set_rownames(location_row_order$Location) %>% 
    set_colnames(c("True", "Predicted")) %>% 
    pheatmap(color = colorRampPalette(red_palette)(20),
             cluster_cols = F,
             cluster_rows = F, display_numbers = T,
             number_format = "%.2f")
  
  heat_2 <- rbind(compartment_true_class_prob[focus_gene,],
                  compartment_class_probabilities[focus_gene,]) %>% 
    t() %>% 
    set_colnames(c("True", "Predicted")) %>% 
    pheatmap(color = colorRampPalette(red_palette)(20),
             cluster_cols = F,
             cluster_rows = F, display_numbers = T, 
             number_format = "%.2f")
  
  
  # Class prob heatmaps plots ----------------------------------------------------------
  
  location_grob <- grobTree(textGrob(paste("log loss = ", 
                                           log_loss_df %>% 
                                             filter(entrezgene == focus_gene) %>% 
                                             pull(Location_Log_Loss) %>% 
                                             round(2)
  ), x=0.9,  y=0.95, hjust = 1,
  gp=gpar(col="red", fontsize=13)))
  
  
  ll_1 <- log_loss_df %>% 
    ggplot(aes(Location_Log_Loss)) +
    geom_histogram() +
    geom_vline(data = log_loss_df %>% 
                 filter(entrezgene == focus_gene), 
               aes(xintercept = Location_Log_Loss),
               color = "red",
               linetype= "dashed") + 
    annotation_custom(location_grob)
  
  
  
  compartment_grob <- grobTree(textGrob(paste("log loss = ", 
                                              log_loss_df %>% 
                                                filter(entrezgene == focus_gene) %>% 
                                                pull(Compartment_Log_Loss) %>% 
                                                round(2)
  ), x=0.9,  y=0.95, hjust = 1,
  gp=gpar(col="red", fontsize=13)))
  
  ll_2 <- log_loss_df %>% 
    ggplot(aes(Compartment_Log_Loss)) +
    geom_histogram() +
    geom_vline(data = log_loss_df %>% 
                 filter(entrezgene == focus_gene), 
               aes(xintercept = Compartment_Log_Loss),
               color = "red",
               linetype= "dashed") + 
    annotation_custom(compartment_grob)
  
  
  pred_plot <- cowplot::plot_grid(heat_1$gtable, ll_1, heat_2$gtable, ll_2, rel_widths = c(0.6, 0.4))
  
  
  
  # Final combine  ----------------------------------------------------------
  final_plot <- cowplot::plot_grid(network_plot, pred_plot, scale = 0.9, 
                                   labels = c(paste(convert_genes(focus_gene) ,"2nd deg. neighbors"), "Embedding class probabilities"))
  
  out_name <- paste("./output/bioplex_subcell_predictions/", convert_genes(focus_gene), "_Location_Log_Loss=",
                    location_log_loss[which(rownames(deepwalk_split$Bioplex3$test) == focus_gene)] %>% round(2),".png", sep = "")
  
  cowplot::save_plot(out_name, final_plot, base_height = 12, base_width = 21)
  
}

map(rownames(deepwalk_split$Bioplex3$test), ~embedding_report_plot(.))
