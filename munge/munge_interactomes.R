#01-Interactomes

#TODO: Fix the bug where cached is constantly rerunning


##### Gygi #####
##Version 1: July 2018
# bioplex3_raw <- read_tsv('./data/raw/gygi/interactomeNetwork_BioPlex_BP3.tsv')
# hct_raw <- read_tsv('./data/raw/gygi/interactomeNetwork_BioPlex_HCT_5k.tsv')
# replicated_raw <- read_tsv('./data/raw/gygi/replicatedEdgeList.tsv')

##Version 2: December 2018
#Clusters
bioplex3_mcl <- read_tsv('./data/raw/gygi/bioplex/v2/BioPlex3_MCL_Clusters_Dec2018.tsv', col_types = "ccc") %>% 
  dplyr::select(-1)
hct_mcl <- read_tsv('./data/raw/gygi/bioplex/v2/BioPlex_HCT_MCL_Clusters_Dec2018.tsv', col_types = "ccc") %>% 
  dplyr::select(-1)
bioplex2_mcl <- read_csv('./data/raw/NIHMS867766-supplement-supp_table6.csv', col_types = "ccc-") %>% 
  dplyr::select(-3) %>% 
  set_colnames(c("clusterID", "geneID"))

# #CCLE
# #TODO: Get non-excel version of this data
# ProjectTemplate::cache("ccle_prot", {
#   ccle_prot <- read_csv('./data/raw/gygi/ccle/20171122_ccle_protein_expression_draft9.csv') %>% 
#     gather(Cell_Line, Value, 49:426) %>% 
#     dplyr::select(HUGO = Gene_Symbol, Description, Uniprot = Uniprot_Acc, Cell_Line, Value) %>% 
#     de_excelify_genes(HUGO) %>% 
#     left_join(hugo_to_entrez_df(.$HUGO, cached = T)) %>% 
#     group_by(Gene, Cell_Line) %>% 
#     slice(1) %>% 
#     ungroup() %>% 
#     na.omit() %>% 
#     dplyr::select(Gene, Cell_Line, Value) %>% 
#     dependr::df.to.mat() %>% 
#     t()
# })


#Interaction networks
bioplex3 <- read_tsv('./data/raw/gygi/bioplex/v2/BioPlex3_InteractionList_Dec2018.tsv', col_types = "ccccn") %>% 
  dplyr::rename(Gene.x = `GeneID A`, Gene.y = `GeneID B`, HUGO.x = `Symbol A`, HUGO.y = `Symbol B`)

hct <- read_tsv('./data/raw/gygi/bioplex/v2/BioPlex_HCT_InteractionList_Dec2018.tsv', col_types = "ccccn") %>% 
  dplyr::rename(Gene.x = `GeneID A`, Gene.y = `GeneID B`, HUGO.x = `Symbol A`, HUGO.y = `Symbol B`)

rep_net <- read_tsv("./data/raw/gygi/bioplex/v2/bioPlex_replicated_network_Dec_2018.tsv", col_types = "cc") %>% 
  dplyr::rename(Gene.x = GeneA, Gene.y = GeneB)

bioplex2 <- read_tsv('data/raw/nature22366s2.tsv', col_types = "ccccccnnn") %>% 
  dplyr::select(1:2, 5:6, 9) %>% 
  set_names(c('Gene.x', 'Gene.y', 'HUGO.x', 'HUGO.y', 'Score')) 

bioplex_master <- full_join(bioplex3 %>% rename(Bioplex3 = Score),
                            hct %>% rename(HCT116 = Score)) %>% 
  left_join(rep_net %>% 
              rbind(rep_net %>% rename(Gene.x = Gene.y, Gene.y = Gene.x)) %>% 
              mutate(Replicated = 1)) %>% 
  mutate(Dataset = case_when(Replicated == 1 ~ "Replicated",
                             !is.na(Bioplex3) & is.na(HCT116) ~ "Bioplex3",
                             !is.na(HCT116) & is.na(Bioplex3) ~ "HCT116",
                             TRUE ~ "Double_Obs"))

#Write out interactomes for ingestion in Cytoscape
community_interactome <- split(bioplex3_mcl$geneID, bioplex3_mcl$clusterID) %>%
  llply(function(x)combn(x, 2) %>% t() %>% as_tibble()) %>%
  enframe(name = "gene_set") %>%
  unnest() %>%
  dplyr::rename(Gene.x = V1, Gene.y = V2) %>%
  inner_join(bioplex3)

# write_tsv(community_interactome, "./output/files/community_interactome.tsv")

#Bait selection
bait_bioplex3 <- read_tsv("./data/raw/gygi/bioplex/v2/baitList_BioPlex3.tsv", col_types = "cc") %>% 
  rename(HUGO = Symbol, Gene = GeneID)

bait_hct <- read_tsv("./data/raw/gygi/bioplex/v2/baitList_BioPlexHCT.tsv", col_types = "cc") %>% 
  rename(HUGO = Symbol, Gene = GeneID)

bait_bioplex2 <- read_csv("./data/raw/NIHMS867766-supplement-supp_table1.csv", col_types = "ccc-") %>% 
  set_colnames(c("HUGO", "Gene", "Num_Interactions"))

bait_meta <- list("Bioplex3" = bait_bioplex3,
                  "HCT116" = bait_hct,
                  "Bioplex2" = bait_bioplex2) %>% 
  enframe("Dataset") %>% 
  unnest() %>% 
  select(-Num_Interactions)


# #Prey abundance TMT data, Ed Huttin, 2-22-2019
# ProjectTemplate::cache("prey_abundance", {
#   prey_abundance <- read_tsv("./data/raw/gygi/bioplex/v2/PQ13551_11CellLines_Reps1-3_Scaled.tsv") %>% 
#     mutate(Uniprot = map_chr(`Protein Id`, ~ word(., 2, sep = "\\|"))) %>% 
#     left_join(uniprot_to_entrez_df(.$Uniprot)) %>% 
#     select(Gene, 9:(ncol(.)-2)) %>% 
#     gather(Key, Abundance, 2:ncol(.)) %>% 
#     separate(Key, into = c("Rep", "Code"), sep = "~") %>% 
#     mutate(Rep = word(Rep, 2)) %>% 
#     filter(!is.na(Code))
#   
#   cell_line_mapping <- tibble(Code = unique(prey_abundance$Code), 
#                               Cell_Line = c("RKO","A549","U87","HCT116","HEK293T","HeLa","MCF7","U2OS","SUM159","PANC1","Jurkat"))
#   
#   prey_abundance <- prey_abundance %>% 
#     left_join(cell_line_mapping) %>% 
#     group_by(Gene, Cell_Line) %>% 
#     summarize(Mean_Abundance = mean(Abundance)) %>% 
#     ungroup() %>% 
#     spread(Cell_Line, Mean_Abundance) %>% 
#     left_join(entrez_to_hugo_df(.$Gene, cached = T))
# })


#####Lists of complexes####
#CORUM
corum <- read_tsv("./data/raw/corum_human_core_complexes.tsv",
         col_types = "icc") %>%
  distinct(Complex, Gene, .keep_all = T) %>%
  add_count(Complex) %>% 
  dplyr::rename(symbol = Gene) %>%
  mutate(entrezgene = convert_genes(symbol, "symbol", "entrez_id")) %>% 
  unique()

corum_list <- split(corum$entrezgene, corum$Complex)


#hu.MAP
humap_list <- scan("./data/raw/humap_clusters.txt", what="", sep="\n")  %>%
  str_split(pattern = " ") %>%
  set_names(1:length(.))

#Bioplex
bioplex_list <- split(bioplex3_mcl$geneID, bioplex3_mcl$clusterID) 
bioplex2_list <- split(bioplex2_mcl$geneID, bioplex2_mcl$clusterID)

#HCT
hct_list <- split(hct_mcl$geneID, hct_mcl$clusterID) 

##### Other Interactomes #####
humap <- read_tsv("./data/raw/humap_ints.txt", col_names = c("Gene.x", "Gene.y", "Confidence"), col_types = "ccn")


#Harper
autophagy <- read_excel("./data/raw/gygi/other/behrends_2010_nature.xls", .name_repair = "universal") %>% 
  mutate(entrezgene = convert_genes(Bait, "symbol", "entrez_id")) %>% 
  mutate(Gene.x = entrezgene, Gene.y = GeneID) %>% 
  select(Gene.x, Gene.y, everything())

dub <- read_excel("./data/raw/gygi/other/sowa_2009_cell.xls", .name_repair = "universal") %>% 
  mutate(entrezgene = convert_genes(BAIT, "symbol", "entrez_id")) %>% 
  mutate(Gene.x = entrezgene, Gene.y = GeneID) %>% 
  select(Gene.x, Gene.y, everything())

#Biogrid
biogrid_human <- read_tsv("./data/raw/BIOGRID-ALL-3.5.176.tab2.txt") %>% 
  filter(`Organism Interactor A` == "9606", `Organism Interactor B` == "9606", `Experimental System` != "Two-hybrid") %>% 
  rename(entrezgene.x = "Entrez Gene Interactor A",
         entrezgene.y = "Entrez Gene Interactor B",
         symbol.x = "Official Symbol Interactor A",
         symbol.y = "Official Symbol Interactor B")

# 
# #Anna Malovannaya 2011 Cell
# malov_expt_meta <- read_csv("./data/raw/Malovannaya_2011_ab_meta.csv", skip = 11) %>% 
#   rename(Bait = `Intended Antigen`)
# malov_gene_meta <- read_csv("./data/raw/Malovannaya_2011_gene_meta.csv", skip = 5)
# 
# malov_data <- read_csv("./data/raw/Malovannaya_2011_all_ms.csv", skip = 3) %>% 
#   gather(Antibody, Count, 2:ncol(.)) %>% 
#   mutate(Antibody = str_sub(Antibody, start = 5)) %>% 
#   left_join(malov_expt_meta) %>% 
#   rename(Gene.y = GeneID) %>% 
#   left_join(malov_gene_meta %>% select(1:2), by = c("Bait" = "GeneSymbol")) %>% 
#   rename(Gene.x = GeneID) %>% 
#   mutate(Gene.x = as.character(Gene.x),
#          Gene.y = as.character(Gene.y))