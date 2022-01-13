#Pubmed_Citations
#Source: wget -N -P /Users/joshpan/gene_function_dev/data/raw/ ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz
#Updated 10/10/2019

ProjectTemplate::cache("pubmed_citations", { 
  read_tsv("./data/raw/gene2pubmed.txt") %>% 
    filter(`#tax_id` == 9606) %>%
    dplyr::rename(entrezgene = GeneID) %>% 
    dplyr::mutate(entrezgene = as.character(entrezgene)) %>% 
    dplyr::count(entrezgene) %>% 
    rename(Pubmed_Count = n)})
