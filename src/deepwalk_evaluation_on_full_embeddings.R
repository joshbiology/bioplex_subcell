## ---------------------------
##
## Script name: deepwalk_evaluation_on_full_embeddings
##
## Purpose of script: Train xgboost classifier using DeepWalk embeddings
##    as input.
##
## Author: Joshua Pan
##
## Date Created: 2021rest
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

#Manually defining compartments based on correlation
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
load("./cache/deepwalk_mat_full.RData")


#Split data

deepwalk_location_df <- map(deepwalk_mat_full[names(networks)], function(x) inner_join(bioid_loc %>% select(NMF_Rank, entrezgene),
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

deepwalk_compartment_df <- map(deepwalk_mat_full[names(networks)], function(x) inner_join(bioid_loc %>% select(Compartment, entrezgene),
                                                                                              as_tibble(x %>% scale(), rownames = "entrezgene")) %>% 
                                 column_to_rownames("entrezgene"))

deepwalk_compartment_df_soft_labels <- map(deepwalk_compartment_df, function(x) bioid_mixing_mat[rownames(x),])

deepwalk_compartment_split <- map2(deepwalk_compartment_df, train_indices, function(x, y) {list(train = x[y,], test = x[-y,])})

deepwalk_compartment_soft_labels_split <- map2(deepwalk_compartment_df_soft_labels, train_indices, function(x, y) {list(train = x[y,], test = x[-y,])})



# Manual determination of test and train genes ----------------------------

#This is to test if there is consistency issues between the runs due to seeding from diff. test and train genes.

intersect_genes <- intersect(deepwalk_mat_full$Bioplex3 %>% rownames(), deepwalk_mat_full$HCT %>% rownames()) %>% 
  intersect(bioid_loc$entrezgene)

train_df <- bioid_loc %>% filter(entrezgene %in% intersect_genes)

#ignore class imbalance for now with compartments, prioritize localizations.

manual_split <- caret::createDataPartition(train_df$NMF_Rank, p = .75, 
                           list = FALSE, 
                           times = 1)

train_genes <- train_df$entrezgene[manual_split]
train_labels <- train_df$NMF_Rank[manual_split]
train_comp_labels <- train_df$Compartment[manual_split]
test_genes <-train_df$entrezgene[-manual_split]
test_labels <-train_df$NMF_Rank[-manual_split]
test_comp_labels <-train_df$Compartment[-manual_split]


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
  x = deepwalk_mat_full$Bioplex3[train_genes,],
  y = train_labels,
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
  x = deepwalk_mat_full$Bioplex3[train_genes,],
  y = train_comp_labels,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE,
  objective='multi:softprob',
  num_class = 6
)

#saveRDS(bioplex3_compartment_xgboost, "./bioplex3_compartment_xgboost.rds")
#bioplex3_compartment_xgboost <- readRDS("./bioplex3_compartment_xgboost.rds")


#HCT116

#Classification by Loc
hct_location_xgboost <- caret::train(
  x = deepwalk_mat_full$HCT[train_genes,],
  y = train_labels,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE,
  objective='multi:softprob',
  num_class = 20
)

#saveRDS(hct_location_xgboost, "./hct_location_xgboost.rds")

#Classification by Compartment
hct_compartment_xgboost <- caret::train(
  x = deepwalk_mat_full$HCT[train_genes,],
  y = train_comp_labels,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE,
  objective='multi:softprob',
  num_class = 6
)
#saveRDS(hct_compartment_xgboost, "./hct_compartment_xgboost.rds")

# Generate class prob ----------------------------------------------------------

location_class_probabilities <- predict(bioplex3_location_xgboost$finalModel, newdata = deepwalk_mat_full$Bioplex3[test_genes,] %>% as.matrix()) 
location_class_probabilities <-  matrix(location_class_probabilities, nrow = 20) %>% t()
rownames(location_class_probabilities) <- test_genes


#Compute class probabilities

compartment_class_probabilities <- predict(bioplex3_compartment_xgboost$finalModel, newdata = deepwalk_mat_full$Bioplex3[test_genes,] %>% as.matrix()) 
compartment_class_probabilities <-  matrix(compartment_class_probabilities, nrow = 6) %>% t()
rownames(compartment_class_probabilities) <- test_genes
colnames(compartment_class_probabilities) <- names(manual_mapping)

#True classes

location_true_class_prob <- t(bioid_mixing_mat) %*% diag(1/rowSums(bioid_mixing_mat)) %>% t() %>% 
  set_rownames(rownames(bioid_mixing_mat))

compartment_true_class_prob <- t(bioid_nmf_mixing_compartment) %*% diag(1/rowSums(bioid_nmf_mixing_compartment)) %>% t() %>% 
  set_rownames(rownames(bioid_nmf_mixing_compartment))


#Evaluate log-loss


compartment_log_loss <- map_dbl(test_genes, ~ModelMetrics::logLoss(compartment_class_probabilities[.,],
                                                                                                              compartment_true_class_prob[.,]))


location_log_loss <- map_dbl(test_genes, ~ModelMetrics::logLoss(location_class_probabilities[.,],
                                                                                                        location_true_class_prob[.,]))

log_loss_df <- tibble(entrezgene = test_genes,
                      symbol = convert_genes(entrezgene),
                      Location_Log_Loss = location_log_loss,
                      Compartment_Log_Loss = compartment_log_loss)

qplot(compartment_log_loss, location_log_loss)



# HCT116 Generate class prob ----------------------------------------------------------

hct_location_class_probabilities <- predict(hct_location_xgboost$finalModel, newdata = deepwalk_mat_full$HCT[test_genes,] %>% as.matrix()) 
hct_location_class_probabilities <-  matrix(hct_location_class_probabilities, nrow = 20) %>% t()
rownames(hct_location_class_probabilities) <- test_genes


#Compute class probabilities

hct_compartment_class_probabilities <- predict(hct_compartment_xgboost$finalModel, newdata = deepwalk_mat_full$HCT[test_genes,] %>% as.matrix()) 
hct_compartment_class_probabilities <-  matrix(hct_compartment_class_probabilities, nrow = 6) %>% t()
rownames(hct_compartment_class_probabilities) <- test_genes
colnames(hct_compartment_class_probabilities) <- names(manual_mapping)


#Evaluate log-loss


hct_compartment_log_loss <- map_dbl(test_genes, ~ModelMetrics::logLoss(hct_compartment_class_probabilities[.,],
                                                                                                              compartment_true_class_prob[.,]))


hct_location_log_loss <- map_dbl(test_genes, ~ModelMetrics::logLoss(hct_location_class_probabilities[.,],
                                                                                                        location_true_class_prob[.,]))

hct_log_loss_df <- tibble(entrezgene = test_genes,
                      symbol = convert_genes(entrezgene),
                      Location_Log_Loss = hct_location_log_loss,
                      Compartment_Log_Loss = hct_compartment_log_loss)

qplot(hct_compartment_log_loss, hct_location_log_loss)



# Agreement ---------------------------------------------------------------

bioplex3_test_summary <- caret::confusionMatrix(test_labels, max.col(location_class_probabilities) %>% factor(levels = 1:20))

hct_test_summary <- caret::confusionMatrix(test_labels, max.col(hct_location_class_probabilities) %>% factor(levels = 1:20))


overall_summary <- cbind(bioplex3_test_summary$byClass[,3], hct_test_summary$byClass[,3]) %>% 
  set_colnames(c("Bioplex3_PPV", "HCT_PPV"))  %>% 
  set_rownames(nmf_rank_to_name[as.character(1:20)]) %>% 
  magrittr::extract(order(bioplex3_test_summary$byClass[,3], decreasing = T),)


overall_summary %>% 
  pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, display_numbers = T)


#Location log loss
inner_join(log_loss_df, hct_log_loss_df %>% rename(HCT_Location_Log_Loss = Location_Log_Loss,  HCT_Compartment_Log_Loss=  Compartment_Log_Loss)) %>% 
  left_join(bioid_loc) %>% 
  ggplot(aes(Location_Log_Loss, HCT_Location_Log_Loss, color = Location)) +
  geom_point() +
  geom_abline(slope=1, intercept = 0)

inner_join(log_loss_df, hct_log_loss_df %>% rename(HCT_Location_Log_Loss = Location_Log_Loss,  HCT_Compartment_Log_Loss=  Compartment_Log_Loss)) %>% 
  left_join(bioid_loc) %>% 
  ggplot(aes(Compartment_Log_Loss, HCT_Compartment_Log_Loss, color = Compartment)) +
  geom_point() +
  geom_abline(slope=1, intercept = 0)

inner_join(log_loss_df, hct_log_loss_df %>% rename(HCT_Location_Log_Loss = Location_Log_Loss,  HCT_Compartment_Log_Loss=  Compartment_Log_Loss)) %>% 
  pivot_longer(names_to = "Dataset", values_to = "Log_Loss",c(Compartment_Log_Loss, HCT_Compartment_Log_Loss)) %>% 
  left_join(bioid_loc) %>% 
  ggplot(aes(Log_Loss)) +
  geom_histogram(aes(fill = Compartment)) +
  facet_grid(Dataset~Compartment)

hct_log_loss_df %>% 
  left_join(bioid_loc) %>% 
  ggplot(aes(Location_Log_Loss)) +
  geom_histogram(aes(fill = Location)) +
  facet_wrap(~Location)

hct_log_loss_df %>% 
  left_join(bioid_loc) %>% 
  ggplot(aes(Compartment_Log_Loss)) +
  geom_histogram(aes(fill = Compartment)) +
  facet_wrap(~Compartment)


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

bioplex3
# 
# library(UniProt.ws)
# up <- UniProt.ws(taxId=9606)
# keys <- union(deepwalk_mat_full$Bioplex3 %>% rownames(), deepwalk_mat_full$HCT %>% rownames())
# columns <- c("SUBCELLULAR-LOCATIONS")
# kt <- "ENTREZ_GENE"
# res <- select(up, keys, columns, kt)

#Uniprot
uniprot_loc <- read_csv("./output/uniprot_loc.csv", col_names = c("entrezgene", "Uniprot_Loc")) %>% 
  group_by(entrezgene) %>% 
  filter(row_number() == 1) %>% 
  ungroup()

#BIOPLEX3 BAIT status
bait_bioplex3 %>% 
  select(entrezgene = Gene) %>% 
  mutate(Is_293T_Bait = T)

#HCT BAIT status
bait_hct %>% 
  select(entrezgene = Gene) %>% 
  mutate(Is_HCT_Bait = T)

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
  list(pubmed_citations,
                          bait_bioplex3 %>% 
                            dplyr::select(entrezgene = Gene) %>% 
                            mutate(Is_293T_Bait = T),
                          bait_hct %>% 
                            dplyr::select(entrezgene = Gene) %>% 
                            mutate(Is_HCT_Bait = T),
                          bioid_loc %>% 
                            dplyr::select(entrezgene, BioID_Location = Location, BioID_Compartment = Compartment),
                          gene_loc_unique %>% dplyr::select(-1) %>% rename_with(~paste("HPA", .x, sep= "_"), -entrezgene),
                          mitocarta %>% 
                            dplyr::select(entrezgene, `MitoCarta3.0_SubMitoLocalization`),
                          uniprot_loc) %>% 
  reduce(full_join, by = "entrezgene")


# Predict unseen ----------------------------------------------------------


bioplex_loc_out <- predict(bioplex3_location_xgboost$finalModel, newdata = deepwalk_mat_full$Bioplex3) 

bioplex_loc_out <-  matrix(bioplex_loc_out, nrow = 20) %>% t()
bioplex_loc_out <- tibble(entrezgene = deepwalk_mat_full$Bioplex3 %>% rownames(),
               symbol = convert_genes(entrezgene)) %>% 
  cbind(bioplex_loc_out %>% set_colnames(nmf_rank_to_name[as.character(1:20)])) 

hct_loc_out <-  predict(hct_location_xgboost$finalModel, newdata = deepwalk_mat_full$HCT) 
hct_loc_out <-  matrix(hct_loc_out, nrow = 20) %>% t()

hct_loc_out <- tibble(entrezgene = deepwalk_mat_full$HCT %>% rownames(),
               symbol = convert_genes(entrezgene)) %>% 
  cbind(hct_loc_out %>% set_colnames(nmf_rank_to_name[as.character(1:20)])) 

tmp1 <- inner_join(bioplex_loc_out %>% rename_with(~paste("293T", .x, sep = "_"), -c(entrezgene, symbol)), 
                  hct_loc_out%>% rename_with(~paste("HCT116", .x, sep = "_"), -c(entrezgene, symbol))) %>% 
  left_join(meta)

write_tsv(tmp1, "/Users/joshpan/Dropbox (Partners HealthCare)/Bioplex Data Science/output/localization_predictions_full.tsv")



bioplex3_comp_out <- predict(bioplex3_compartment_xgboost$finalModel, newdata = deepwalk_mat_full$Bioplex3)
bioplex3_comp_out <-  matrix(bioplex3_comp_out, nrow = 6) %>% t()

bioplex3_comp_out <- tibble(entrezgene = deepwalk_mat_full$Bioplex3 %>% rownames(),
               symbol = convert_genes(entrezgene)) %>% 
  cbind(bioplex3_comp_out %>% set_colnames(names(manual_mapping)))


hct_comp_out <- predict(bioplex3_compartment_xgboost$finalModel, newdata = deepwalk_mat_full$HCT)
hct_comp_out <-  matrix(hct_comp_out, nrow = 6) %>% t()

hct_comp_out <- tibble(entrezgene = deepwalk_mat_full$HCT %>% rownames(),
                            symbol = convert_genes(entrezgene)) %>% 
  cbind(hct_comp_out %>% set_colnames(names(manual_mapping)))




tmp2 <- inner_join(bioplex3_comp_out %>% rename_with(~paste("293T", .x, sep = "_"), -c(entrezgene, symbol)), 
                  hct_comp_out%>% rename_with(~paste("HCT116", .x, sep = "_"), -c(entrezgene, symbol))) %>% 
  left_join(meta)



write_tsv(tmp2, "/Users/joshpan/Dropbox (Partners HealthCare)/Bioplex Data Science/output/compartment_predictions_full.tsv")



# Agreement 2 -------------------------------------------------------------

loc_compare <- rbind(bioplex_loc_out %>% pivot_longer(names_to = "Location", values_to = "Prob", -c(entrezgene, symbol)) %>% mutate(Dataset = "HEK293T"),
                     hct_loc_out %>% pivot_longer(names_to = "Location", values_to = "Prob", -c(entrezgene, symbol)) %>% mutate(Dataset = "HCT116")) %>% 
  mutate(Training = case_when(entrezgene %in% train_genes ~ "Training",
                              entrezgene %in% test_genes ~ "Held_Out_Test",
                              TRUE ~ "Unlabeled") %>% factor(levels = c("Training", "Held_Out_Test","Unlabeled")))

#labeled data
train_df %>% 
  left_join(loc_compare %>% filter(Training != "Unlabeled")) %>% 
  pivot_wider(names_from = "Dataset", values_from = Prob) %>% 
  filter(!is.na(Training)) %>% 
  mutate(Location = factor(Location, levels= nmf_rank_to_name[order(bioplex3_test_summary$byClass[,3], decreasing = T) %>% as.character()])) %>% 
  ggplot(aes(HEK293T, HCT116, color = Location)) +
  geom_point() +
  facet_wrap(~Location+Training, ncol = 8)


#Unlabeled data
loc_compare %>% filter(Training == "Unlabeled") %>% 
  pivot_wider(names_from = "Dataset", values_from = Prob) %>% 
  filter(!is.na(Training)) %>% 
  mutate(Location = factor(Location, levels= nmf_rank_to_name[order(bioplex3_test_summary$byClass[,3], decreasing = T) %>% as.character()])) %>% 
  ggplot(aes(HEK293T, HCT116, color = Location)) +
  geom_point() +
  geom_vline(xintercept = 0.4, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.4, color = "red", linetype = "dashed") +
  facet_wrap(~Location+Training, ncol = 8)



loc_compare %>% 
  pivot_wider(names_from = "Dataset", values_from = Prob) %>% 
  ggplot(aes(HEK293T, HCT116, color = Location)) +
  geom_point() +
  facet_wrap(~Location+Training, ncol = 9)


comp_compare <- rbind(bioplex3_comp_out   %>% pivot_longer(names_to = "Compartment", values_to = "Prob", -c(entrezgene, symbol)) %>% mutate(Dataset = "HEK293T"),
                      hct_comp_out %>% pivot_longer(names_to = "Compartment", values_to = "Prob", -c(entrezgene, symbol)) %>% mutate(Dataset = "HCT116"))

comp_compare %>% 
  pivot_wider(names_from = "Dataset", values_from = Prob) %>% 
  ggplot(aes(HEK293T, HCT116, color = Compartment)) +
  geom_point() +
  facet_wrap(~Compartment)


