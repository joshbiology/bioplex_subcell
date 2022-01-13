#Ingest and plot Laura's quality control data

library(tidyverse)
library(readxl)
xl_data <- './data/U2OS_IF_QC.xlsx'
tab_names <- excel_sheets(path = xl_data)
list_all <- lapply(tab_names, function(x) read_excel(path = xl_data, sheet = x, col_types = c(rep("guess", 4), "date", rep("guess", 5))))

qc <- enframe(list_all) %>% unnest() %>% select(-name) %>% tibble(.name_repair = "universal")

colnames(qc)

qc <- rename(qc, Is_High_Quality = IF.Quality..suitable.for.validation.) 

qc <- qc %>% 
  mutate(Localization = str_trim(Localization) %>% tolower()) 


#Summary of U2OS001-007

qc %>% 
  filter(Is_High_Quality %in% c('Yes', 'No')) %>% 
  group_by(Plate) %>% 
  count(Is_High_Quality) %>% 
  ggplot(aes(x = Plate, y = n, fill =  Is_High_Quality)) +
  geom_bar(position="stack", stat = "identity")

#Stats
qc %>% count(IF.Quality..suitable.for.validation.)

#Plot types.
qc %>% 
  filter(Is_High_Quality == "Yes") %>% 
  count(Localization) %>% View()

write_tsv(qc %>% 
            filter(Is_High_Quality == "Yes"), "./out.tsv")


convert_genes_mygeneinfo <- function(genes, to = c("entrezgene", "symbol")) {
  to <- match.arg(to)
  
  dat <- mygene::queryMany(genes, 
                           scopes = "entrezgene,symbol,alias, ensembl.gene", 
                           fields = "entrezgene,name,symbol,taxid,type_of_gene", 
                           species="human", 
                           size = 1,
                           entrezonly = T) %>% 
    as_tibble()
  
  key <- pull(dat, query)
  
  value <- pull(dat, to)
  
  dict <- value %>% magrittr::set_names(key)
  
  return(dict[genes] %>% as.character())
}

#Rename genes
source("/Users/joshpan/gene_fn/lib/gene_id_conversion.R")
library(mygene)
qc <- qc %>% 
  mutate(entrezgene = convert_genes_mygeneinfo(Bait.Symbol))

write_tsv(qc, "./laura_annotations.tsv")


pred <- read_excel('./data/localization_predictions_full.xlsx')

pred <- mutate(pred, entrezgene = as.character(entrezgene))

pred %>% 
  inner_join(qc %>% filter(Is_High_Quality == "Yes")) %>% 
  write_tsv("./test.tsv")

View(pred %>% 
       inner_join(qc %>% filter(Is_High_Quality == "Yes")))

meta <- pred %>% 
       inner_join(qc %>% filter(Is_High_Quality == "Yes")) %>% 
  mutate(Pubmed_Count = as.numeric(Pubmed_Count)) %>% 
  arrange(Pubmed_Count)

pred_293 <- pred %>% column_to_rownames('entrezgene') %>% 
  dplyr::select(starts_with("293T")) %>% 
  as.matrix()

categories <- (colnames(pred_293)) %>% str_sub(6, -1)


pred_hct <- pred %>% column_to_rownames('entrezgene') %>% 
  dplyr::select(starts_with("HCT")) %>% 
  as.matrix()

colnames(pred_293) <- categories
colnames(pred_hct) <- categories


library(cowplot)
#Make text report file:


#Make heatmap 293T


plot_subcell_report <- function(index) {
  
  
  report <- c(header = paste("************I. Gene information************"),
              entrez = sprintf('entrez: %s', meta$entrezgene[index]),
              hugo = sprintf('symbol: %s', meta$symbol[index]),
              pubmed = sprintf('pubmed count: %d', meta$Pubmed_Count[index]),
              sep1 = paste("************II. Bioplex annotations************"),
              laura = sprintf('coverslip prediction: %s', str_wrap(meta$Localization[index], width = 40)),
              plate = sprintf('plate ID: %s', meta$Plate[index]),
              exp = sprintf('expression: %s', meta$Bait.Expression..choose..present..high.expression..no.expression.[index]),
              comments = sprintf('comments: %s', meta$Loc..Comments[index]),
              date = sprintf('date: %s', meta$Date[index]),
              sep2 = paste("************III. Published annotations************"),
              uniprot = (sprintf('uniprot: %s', str_wrap(meta$Uniprot_Loc[index] %>% str_sub(1,40), width = 40))),
              gingras = sprintf('Gingras prox. labeling: %s', meta$BioID_Location[index]),
              hpa = sprintf('Human prot. atlas: %s', meta$HPA_Cell_Structure[index]),
              hpa = sprintf('HPA reliability: %s', meta$HPA_Reliability[index]))
  
  report_g <- ggplot() + theme_void() + xlim(c(0,1)) + draw_text(report, x = rep(0,length(report)), y = rev(seq_len(length(report))), hjust = 0)
  
  
  df_293 <- pred_293[meta$entrezgene[index],] %>% 
    enframe() %>% 
    mutate(name = factor(name, levels = rev(categories)))
  
  g_293 <- df_293 %>% 
    ggplot(aes(1, name, fill = value)) + 
    geom_tile() +
    geom_text(aes(label = round(value, 1))) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1)) +
    theme_minimal() +
    scale_y_discrete(labels = function(x) str_wrap(rev(categories), width = 30)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "bottom") +
    labs(fill = "293T")
  
  
  
  df_hct <- pred_hct[meta$entrezgene[index],] %>% 
    enframe() %>% 
    mutate(name = factor(name, levels = rev(categories)))
  
  g_hct <- df_hct %>% 
    ggplot(aes(1, name, fill = value)) + 
    geom_tile() +
    geom_text(aes(label = round(value, 1))) +
    
    scale_fill_gradient(low = "white", high = "purple", limits = c(0,1)) + 
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "bottom")  +
    labs(fill = "HCT")
  
  
  cowplot::plot_grid(report_g, g_293, g_hct, nrow = 1, rel_widths = c(1, 0.6, 0.3))
  
  
  ggsave(file.path(".", "output", "pred", paste(index, "_", meta$symbol[index], ".pdf", sep = "")), width = 10, height = 6)
  
  
}

map(seq_len(nrow(meta)), plot_subcell_report)
