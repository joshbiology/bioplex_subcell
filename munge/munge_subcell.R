#Subcellular localization


# HPA ---------------------------------------------------------------------


#hierarchy from this mapping: https://www.proteinatlas.org/humancell
loc_hier <- read_csv("./data/raw/subcell_hierarchy.csv")

#The Human Protein Atlas version 18 and Ensembl version 88.38 (https://www.proteinatlas.org/about/download)
#reliability scores described here: https://www.proteinatlas.org/about/assays+annotation#ifre
#https://stackoverflow.com/questions/15347282/split-delimited-strings-in-a-column-and-insert-as-new-rows

gene_loc <- read_tsv("./data/raw/subcellular_location.tsv") %>%
  dplyr::select(1, 4:7) %>%
  dplyr::rename(ensembl = 'Gene') %>%
  gather(Reliability, Cell_Structure, 2:ncol(.)) %>%
  filter(!is.na(Cell_Structure)) %>%
  separate_rows(Cell_Structure, sep = ";") %>%
  left_join(loc_hier) %>% 
  mutate(entrezgene = convert_genes(ensembl, from = "ensembl_gene_id", to = "entrez_id")) %>% 
  na.omit()

gene_loc_unique <- gene_loc %>%
  mutate(Reliability = factor(Reliability, levels = c("Enhanced",  "Supported", "Approved", "Uncertain"))) %>% 
  arrange(Reliability) %>% 
  group_by(entrezgene) %>% 
  filter(row_number() == 1) %>% #this is alphabetical at the moment...
  ungroup()


# HumanCellMap (bioid) -----------------------------------------------------------

#BioID gene map: https://humancellmap.org
bioid_loc <- read_tsv("./data/raw/humancellmap/preys-latest.txt") %>% 
  rename(entrezgene = Entrez, Location = `MMF localization`, NMF_Rank = `NMF rank`) %>%
  mutate(entrezgene = as.character(entrezgene))

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


# Hela fractionation ------------------------------------------------------


#Hela fractionated proteomes: http://mapofthecell.biochem.mpg.de
# hela_frac <- read_excel("./data/raw/itzhak_2016_elife_subcell.xlsx", sheet = 2) %>% 
#   mutate(entrezgene = convert_genes_mygeneinfo(`Gene name`, to = "entrezgene")) %>% 
#   select(entrezgene, everything())


# Mitocarta ---------------------------------------------------------------

mitocarta <- read_excel("./data/raw/mootha/Human.MitoCarta3.0.xls", sheet =2) %>% 
  rename(entrezgene = HumanGeneID) %>%
  mutate(entrezgene = as.character(entrezgene))


export_flag_subcell <- F
if (export_flag_subcell) {
  write_tsv(gene_loc, "./output/localization/gene_loc.tsv")
  write_tsv(bioid_loc, "./output/localization/bioid_loc.tsv")
  write_tsv(hela_frac, "./output/localization/hela_frac.tsv")
}

# Preprocess Gringras data ------------------------------------------------


# Visualizations on mixing matrix ---------------------------------------------

all_genes <- bioid_loc %>% arrange(NMF_Rank) %>% pull(entrezgene) %>% na.omit()
mito_genes <- bioid_loc %>% filter(NMF_Rank == 4) %>% pull(entrezgene) %>% na.omit()


#pheatmap::pheatmap(bioid_mixing_mat[intersect(all_genes, rownames(bioid_mixing_mat)),] %>% 
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

