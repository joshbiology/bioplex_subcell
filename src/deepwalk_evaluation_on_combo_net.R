## ---------------------------
##
## Script name: deepwalk_evaluation_on_combo_net

## Purpose of script: Train xgboost classifier using DeepWalk embeddings
##    as input.
##
## Author: Joshua Pan
##
## Date Created: 2021
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

load('./cache/pubmed_citations.RData')

# prepare networks and import embedding -------------------------------------------------------

combo_net_filtered <- combo_net %>% 
  filter(Detected > 1)%>% 
  filter(Gene.x %in% bioid_loc$entrezgene & Gene.y %in% bioid_loc$entrezgene)

#Import embeddings
deepwalk <- read_delim("./data/interim/deepwalk_embeddings/2022-01-09/Combo.full.deepwalk.txt",
                            delim = " ", 
                            skip = 1, 
                            col_names = F) %>%
  rename(name = X1) %>% 
  column_to_rownames("name") %>% 
  as.matrix()


#Split data - holding only those with bioid labels.
deepwalk_location_df <- inner_join(bioid_loc %>% select(NMF_Rank, entrezgene),
                                   as_tibble(deepwalk %>% scale(), rownames = "entrezgene")) %>% 
                              column_to_rownames("entrezgene")

deepwalk_location_df_soft_labels <- bioid_mixing_mat[rownames(deepwalk_location_df),]

#Freeze for both location and compartment
train_indices <- caret::createDataPartition(deepwalk_location_df$NMF_Rank, p = .75, 
                                                                       list = FALSE, 
                                                                       times = 1)

deepwalk_location_split <- list(train = deepwalk_location_df[train_indices,], test = deepwalk_location_df[-train_indices,])

deepwalk_location_soft_labels_split <- list(train = deepwalk_location_df_soft_labels[train_indices,], test = deepwalk_location_df_soft_labels[-train_indices,])


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
combo_location_xgboost <- caret::train(
  x = deepwalk_location_split$train[,-1],
  y = deepwalk_location_split$train$NMF_Rank,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE,
  objective='multi:softprob',
  num_class = 20
)

#saveRDS(combo_location_xgboost, "./combo_location_xgboost.Rds")


# Generate class prob ----------------------------------------------------------

location_class_probabilities <- predict(combo_location_xgboost$finalModel, newdata = deepwalk_location_split$test[,-1] %>% as.matrix) %>% as.matrix()
location_class_probabilities <-  matrix(location_class_probabilities, nrow = 20) %>% t()
rownames(location_class_probabilities) <- deepwalk_location_split$test %>% rownames()


#True classes

location_true_class_prob <- t(bioid_mixing_mat) %*% diag(1/rowSums(bioid_mixing_mat)) %>% t() %>% 
  set_rownames(rownames(bioid_mixing_mat))

#Evaluate log-loss
location_log_loss <- map_dbl(intersect(rownames(location_class_probabilities),
                                       rownames(location_true_class_prob)), ~ModelMetrics::logLoss(location_class_probabilities[.,],
                                                                                                        location_true_class_prob[.,]))

log_loss_df <- tibble(entrezgene = deepwalk_location_split$test %>% rownames(),
                      symbol = convert_genes(entrezgene),
                      Location_Log_Loss = location_log_loss)


# Agreement ---------------------------------------------------------------

combo_test_summary <- caret::confusionMatrix(deepwalk_location_split$test[,1], max.col(location_class_probabilities) %>% factor(levels = 1:20))

rank <- combo_test_summary$byClass[,"Precision"] %>% order()

overall_summary <- combo_test_summary$table  %>% 
  set_rownames(nmf_rank_to_name[as.character(1:20)]) %>% 
  set_colnames(nmf_rank_to_name[as.character(1:20)] %>% str_sub(1,20))

overall_summary[rev(rank),rev(rank)] %>% 
  pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, 
                     display_numbers = T, number_format = "%.0f",
                     ) 



# Generate plots ----------------------------------------------------------
focus_gene <- "55037"

embedding_report_plot <- function(focus_gene, save_plot = F) {
  
  
  
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
                    location_log_loss[which(rownames(deepwalk_location_split$Bioplex3$test) == focus_gene)] %>% round(2),".png", sep = "")
  
  if(save_plot) cowplot::save_plot(out_name, final_plot, base_height = 12, base_width = 21)
  
  else
    final_plot
}
embedding_report_plot("57506")
map(rownames(deepwalk_split$Bioplex3$test), ~embedding_report_plot(.))


# Aggregate info ----------------------------------------------------------

combo_net_meta

#Uniprot
uniprot_loc <- read_csv("./data/raw/uniprot_loc.csv", col_names = c("entrezgene", "Uniprot_Loc")) %>% 
  group_by(entrezgene) %>% 
  filter(row_number() == 1) %>% 
  ungroup()

#bioid_loc
bioid_loc %>% 
  select(entrezgene, BioID_Location = Location, BioID_Compartment = Compartment)

#pubmed citations
pubmed_citations

#subcell atlas unique gene_loc_unique
gene_loc_unique %>% select(-1) %>% rename_with(~paste("HPA", .x, sep= "_"), -entrezgene)

#Mitocarta
mitocarta %>% 
  select(entrezgene, `MitoCarta3.0_SubMitoLocalization`)


meta <-  
  list(pubmed_citations, combo_net_meta %>% mutate(entrezgene = as.character(GeneID)) %>% select(-GeneID),
                          bioid_loc %>% 
                            dplyr::select(entrezgene, BioID_Location = Location, BioID_Compartment = Compartment),
                          gene_loc_unique %>% dplyr::select(-1) %>% rename_with(~paste("HPA", .x, sep= "_"), -entrezgene),
                          mitocarta %>% 
                            dplyr::select(entrezgene, `MitoCarta3.0_SubMitoLocalization`),
                          uniprot_loc) %>% 
  reduce(full_join, by = "entrezgene")


# Predict unseen ----------------------------------------------------------


bioplex_loc_out <- predict(combo_location_xgboost$finalModel, newdata = deepwalk)

bioplex_loc_out <-  matrix(bioplex_loc_out, nrow = 20) %>% t()
bioplex_loc_out <- tibble(entrezgene = deepwalk %>% rownames(),
               symbol = convert_genes(entrezgene)) %>% 
  cbind(bioplex_loc_out %>% set_colnames(nmf_rank_to_name[as.character(1:20)])) 


tmp1 <- bioplex_loc_out %>% 
  left_join(meta)

write_tsv(tmp1, "./output/combo_localization_predictions_full.tsv")




