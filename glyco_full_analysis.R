# set up workspace
assign('last.warning', NULL, envir = baseenv())
rm(list = ls(all.names = TRUE))
sessionInfo()
closeAllConnections()
gc() 

#install.packages('seqinr')
#install.packages('readxl')

setwd('~/Dropbox/Broad-Research/ChemBio/R_Drpbx/sleep/glyco/')

library(tidyverse)
library(readxl)
library(seqinr)


## IDENTIFIER MAPPINGS

uniprotid_geneid_map_flybase <-read_tsv('fbgn_NAseq_Uniprot_fb_2020_04.tsv', skip = 4) %>%
  dplyr::select(1,3,6) %>%
  dplyr::rename(fly_gene_symbol=1, fbgn_id=2, uniprot_id=3) %>%
  filter(!is.na(uniprot_id))

uniprotid_geneid_map_uniprot <-read_tsv('DROME_7227_idmapping_selected.tab', col_names = F) %>% 
  dplyr::select(1,2) %>% 
  dplyr::rename(uniprot_id = 1, uniprot_name = 2)

## CAUTION !!
# different genes with same uniprot IDs
uniprotid_geneid_map_flybase %>% group_by(uniprot_id) %>% filter(n()>1)
# different uniprot IDs with same gene symbols
uniprotid_geneid_map_flybase %>% group_by(fly_gene_symbol) %>% filter(n()>1)
# fix different genes with same uniprot IDs
uniprotid_geneid_map_flybase <- uniprotid_geneid_map_flybase %>% 
  group_by(uniprot_id) %>% 
  summarize_all(list(~toString(unique(.)))) 

## PROTEOME

tbl_proteome_1 <-read_excel('prot_matrix_RNAi1_RNAi2.xlsx', sheet = 'prot_matrix_RNAi1')
tbl_proteome_2 <-read_excel('prot_matrix_RNAi1_RNAi2.xlsx', sheet = 'prot_matrix_RNAi2') 

tbl_proteome <- full_join(tbl_proteome_1, tbl_proteome_2) %>% 
  dplyr::rename(uniprot_id = 1)


tbl_proteome %>% group_by(uniprot_id) %>% filter(n()>1)

tbl_proteome <- left_join(left_join(tbl_proteome, uniprotid_geneid_map_flybase), uniprotid_geneid_map_uniprot) %>% 
  dplyr::select(12,1,14,13, 2:11) %>% 
  dplyr::rename('Ctrl-1' = 'WT-1', 
         'Ctrl-2' = 'WT-2', 
         'Ctrl-3' = 'WT-3',
         'Ctrl-4' = 'WT-4',
         
         'Alg10_RNAi_#1-1' = 'RNAi #1-1',
         'Alg10_RNAi_#1-2' = 'RNAi #1-2',
         'Alg10_RNAi_#1-3' = 'RNAi #1-3',
         
         'Alg10_RNAi_#2-1' = 'RNAi #2-1',
         'Alg10_RNAi_#2-2' = 'RNAi #2-2',
         'Alg10_RNAi_#2-3' = 'RNAi #2-3')  

tbl_proteome <- tbl_proteome %>% 
  add_column(ttest_pval = NA, .after = 4) %>% 
  add_column(fc_median = NA, .after = 5) %>% 
  add_column(fc_mean = NA, .after = 6) %>% 
  add_column(p_fdr = NA, .after = 7) 



idx_ctrl_cols <- which(str_detect(names(tbl_proteome), 'Ctrl'))
idx_alg10_cols <- which(str_detect(names(tbl_proteome), 'Alg10'))

for(i in 1:nrow(tbl_proteome)){
  tbl_proteome$ttest_pval[i] <- {t.test(tbl_proteome[i,idx_alg10_cols],tbl_proteome[i,idx_ctrl_cols])}$p.val
  tbl_proteome$fc_median[i] <- median(as.numeric(tbl_proteome[i,idx_alg10_cols]),na.rm=T)/median(as.numeric(tbl_proteome[i,idx_ctrl_cols]), na.rm=T)
  tbl_proteome$fc_mean[i] <- mean(as.numeric(tbl_proteome[i,idx_alg10_cols]),na.rm=T)/mean(as.numeric(tbl_proteome[i,idx_ctrl_cols]), na.rm=T)
}


tbl_proteome$p_fdr <- p.adjust(tbl_proteome$ttest_pval, method='fdr')

tbl_proteome <- tbl_proteome %>%
  arrange(ttest_pval)


## GLYCOPROTEOME - GENE LEVEL
tbl_glycoproteome_1 <-read_excel('prot_matrix_WT4_RNAi1_RNAi2.xlsx', sheet = 'prot_matrix_WT4_RNAi1') %>% 
  dplyr::select(1:8) %>% 
  dplyr::rename(uniprot_id = 1)
tbl_glycoproteome_2 <-read_excel('prot_matrix_WT4_RNAi1_RNAi2.xlsx', sheet = 'prot_matrix_WT4_RNAi2') %>% 
  dplyr::select(1:8) %>% 
  dplyr::rename(uniprot_id = 1)

# WT- 126, 127N, 127C, 128N;Alg10 RNAi #1- 128C, 129N, 129C;
# WT- 126, 127N, 127C, 128N;Alg10 RNAi #2; 130N, 130C, 131

tbl_glycoproteome <- full_join(tbl_glycoproteome_1, tbl_glycoproteome_2) %>% 
  dplyr::rename('Ctrl-1' = 'F_126', 
         'Ctrl-2' = 'F_127N', 
         'Ctrl-3' = 'F_127C',
         'Ctrl-4' = 'F_128N',
         
         'Alg10_RNAi_#1-1' = 'F_128C',
         'Alg10_RNAi_#1-2' = 'F_129N',
         'Alg10_RNAi_#1-3' = 'F_129C',
         
         'Alg10_RNAi_#2-1' = 'F_130N',
         'Alg10_RNAi_#2-2' = 'F_130C',
         'Alg10_RNAi_#2-3' = 'F_131N')

tbl_glycoproteome <- left_join(left_join(tbl_glycoproteome, uniprotid_geneid_map_flybase), 
                               uniprotid_geneid_map_uniprot) 

tbl_glycoproteome <- tbl_glycoproteome %>% 
  dplyr::select(1,12,13,14,2:11) %>% 
  add_column(ttest_pval = NA, .after = 4) %>% 
  add_column(fc_median = NA, .after = 5) %>% 
  add_column(fc_mean = NA, .after = 6) %>% 
  add_column(p_fdr = NA, .after = 7) 


idx_ctrl_cols <- which(str_detect(names(tbl_glycoproteome), 'Ctrl'))
idx_alg10_cols <- which(str_detect(names(tbl_glycoproteome), 'Alg10'))


for(i in 1:nrow(tbl_glycoproteome)){
  tbl_glycoproteome$ttest_pval[i] <- {t.test(tbl_glycoproteome[i,idx_alg10_cols],tbl_glycoproteome[i,idx_ctrl_cols])}$p.val
  tbl_glycoproteome$fc_median[i] <- median(as.numeric(tbl_glycoproteome[i,idx_alg10_cols]),na.rm=T) / 
    median(as.numeric(tbl_glycoproteome[i,idx_ctrl_cols]),na.rm=T)
  tbl_glycoproteome$fc_mean[i] <- mean(as.numeric(tbl_glycoproteome[i,idx_alg10_cols]),na.rm=T) / 
    mean(as.numeric(tbl_glycoproteome[i,idx_ctrl_cols]),na.rm=T)
}

tbl_glycoproteome$p_fdr <- p.adjust(tbl_glycoproteome$ttest_pval, method='fdr')



## GLYCOPROTEOME - SITE LEVEL
tbl_glycoproteome_sitelevel_1 <-read_excel('siteLevel_matrix_WT4_RNAi1_RNAi2.xlsx', sheet = 'RNAi1') %>% 
  dplyr::select(1:9) %>% 
  dplyr::rename(uniprot_id = 1, annotated_sequence = 2)
tbl_glycoproteome_sitelevel_2 <-read_excel('siteLevel_matrix_WT4_RNAi1_RNAi2.xlsx', sheet = 'RNAi2') %>% 
  dplyr::select(1:9) %>% 
  dplyr::rename(uniprot_id = 1, annotated_sequence = 2)

tbl_glycoproteome_sitelevel <- full_join(tbl_glycoproteome_sitelevel_1, tbl_glycoproteome_sitelevel_2) %>% 
  dplyr::rename('Ctrl-1' = 'F_126', 
         'Ctrl-2' = 'F_127N', 
         'Ctrl-3' = 'F_127C',
         'Ctrl-4' = 'F_128N',
         
         'Alg10_RNAi_#1-1' = 'F_128C',
         'Alg10_RNAi_#1-2' = 'F_129N',
         'Alg10_RNAi_#1-3' = 'F_129C',
         
         'Alg10_RNAi_#2-1' = 'F_130N',
         'Alg10_RNAi_#2-2' = 'F_130C',
         'Alg10_RNAi_#2-3' = 'F_131N')

tbl_glycoproteome_sitelevel <- left_join(left_join(tbl_glycoproteome_sitelevel, uniprotid_geneid_map_flybase), 
                                   uniprotid_geneid_map_uniprot) 

tbl_glycoproteome_sitelevel <- tbl_glycoproteome_sitelevel %>% 
  dplyr::select(13,1,15,14, 2:12) %>%
  add_column(ttest_pval = NA, .after = 5) %>%
  add_column(fc_median = NA, .after = 6) %>%
  add_column(fc_mean = NA, .after = 7) %>%
  add_column(p_fdr = NA, .after = 8)


idx_ctrl_cols <- which(str_detect(names(tbl_glycoproteome_sitelevel), 'Ctrl'))
idx_alg10_cols <- which(str_detect(names(tbl_glycoproteome_sitelevel), 'Alg10'))


for(i in 1:nrow(tbl_glycoproteome_sitelevel)){
  tbl_glycoproteome_sitelevel$ttest_pval[i] <- {t.test(tbl_glycoproteome_sitelevel[i,idx_alg10_cols],
                                                       tbl_glycoproteome_sitelevel[i,idx_ctrl_cols])}$p.val
  tbl_glycoproteome_sitelevel$fc_median[i] <- median(as.numeric(tbl_glycoproteome_sitelevel[i,idx_alg10_cols]),na.rm=T)/median(as.numeric(tbl_glycoproteome_sitelevel[i,idx_ctrl_cols]), na.rm=T)
  tbl_glycoproteome_sitelevel$fc_mean[i] <- mean(as.numeric(tbl_glycoproteome_sitelevel[i,idx_alg10_cols]),na.rm=T)/mean(as.numeric(tbl_glycoproteome_sitelevel[i,idx_ctrl_cols]), na.rm=T)
}

tbl_glycoproteome_sitelevel$p_fdr <- p.adjust(tbl_glycoproteome_sitelevel$ttest_pval, method='fdr')

tbl_glycoproteome_sitelevel_hits <- tbl_glycoproteome_sitelevel %>% 
  filter(p_fdr<0.25)


## JOINT ANALYSES

tbl_prot_join <- tbl_proteome %>%
  setNames(paste0('prot_', names(.))) %>%
  dplyr::rename('fly_gene_symbol' = 1, 'uniprot_id' = 2, 'uniprot_name' = 3, 'fbgn_id' = 4)

tbl_glycoprot_join <- tbl_glycoproteome %>%
  setNames(paste0('glycoprot_', names(.))) %>%
  dplyr::rename('uniprot_id' = 1, 'fly_gene_symbol' = 2, 'fbgn_id' = 3, 'uniprot_name' = 4)

tbl_glycoprot_sitelevel_join <- tbl_glycoproteome_sitelevel %>%
  setNames(paste0('glycoprotsitelevel_', names(.))) %>%
  dplyr::rename('fly_gene_symbol' = 1, 'uniprot_id' = 2, 'uniprot_name' = 3, 'fbgn_id' = 4,
                'annotated_sequence' = 5)


tbl_joined_sitelevel <- full_join(full_join(tbl_glycoprot_sitelevel_join, tbl_glycoprot_join, by='uniprot_id'), 
                                  tbl_glycoprot_join, by='uniprot_id') %>% 
  arrange(glycoprotsitelevel_ttest_pval) %>% 
  dplyr::select(2,1,4,3,5:9,23:26,10:19,25:34) %>% 
  dplyr::rename('uniprot_id' = 1, 'fly_gene_symbol' = 2, 'fbgn_id' = 3, 'uniprot_name' = 4) %>% 
  unite('peptide_identifier',fly_gene_symbol,annotated_sequence, remove = FALSE)

tbl_joined_genelevel <- full_join(tbl_glycoprot_join, tbl_prot_join, by='uniprot_id') %>%  
  arrange(glycoprot_ttest_pval) %>% 
  dplyr::select(1:8,22:25,9:18,26:35) %>% 
  dplyr::rename('uniprot_id' = 1, 'fly_gene_symbol' = 2, 'fbgn_id' = 3, 'uniprot_name' = 4)






 


## PLOTS

# glycoprot gene level - mean
fdr_thresh <- 0.25
fc_thresh <- 2^0.5

tbl_cat <- tbl_joined_genelevel %>%
  dplyr::select(2) %>% 
  add_column(hit_category = 'NONE', .after = 1)
tbl_cat$hit_category[tbl_joined_genelevel$glycoprot_p_fdr<fdr_thresh]  <- 'FDR_SIGNIFICANT'
tbl_cat$hit_category[tbl_joined_genelevel$glycoprot_p_fdr<fdr_thresh & 
                       tbl_joined_genelevel$glycoprot_fc_mean >= fc_thresh]  <- 'UP'
tbl_cat$hit_category[tbl_joined_genelevel$glycoprot_p_fdr<fdr_thresh & 
                       tbl_joined_genelevel$glycoprot_fc_mean <= 1/fc_thresh]  <- 'DOWN'


point_labels <- tbl_joined_genelevel$fly_gene_symbol
point_labels[tbl_cat$hit_category == 'NONE'] <- ''
point_labels[tbl_cat$hit_category == 'FDR_SIGNIFICANT'] <- ''

hitlist_genelevel_fcmean <- tbl_cat %>% 
  filter(hit_category == 'UP' | hit_category == 'DOWN') %>% 
  arrange(hit_category) %>% 
  dplyr::rename(hit_category_genelevel_fcmean = 2)

p_glycoprot_genelevel_mean <- ggplot(tbl_joined_genelevel, aes(x = log2(glycoprot_fc_mean), 
                                                           y = -log10(glycoprot_ttest_pval))) +
  geom_point(aes(color = tbl_cat$hit_category)) + 
  geom_text(aes(label=point_labels),hjust=1, vjust=0.35) +
  geom_vline(xintercept = log2(fc_thresh), linetype = 'dashed') +
  geom_vline(xintercept = -log2(fc_thresh), linetype = 'dashed') +
  geom_hline(yintercept = - log10(0.05), linetype = 'dashed') +
  xlim(-2, 2) +
  ylim(0, 3.5) + 
  ggtitle('Gene Level') +
  scale_fill_brewer(palette='Set1')

p_glycoprot_genelevel_mean


# glycoprot gene level - median
fdr_thresh <- 0.25
fc_thresh <- 2^0.5

tbl_cat <- tbl_joined_genelevel %>%
  dplyr::select(2) %>% 
  add_column(hit_category = 'NONE', .after = 1)
tbl_cat$hit_category[tbl_joined_genelevel$glycoprot_p_fdr<fdr_thresh]  <- 'FDR_SIGNIFICANT'
tbl_cat$hit_category[tbl_joined_genelevel$glycoprot_p_fdr<fdr_thresh & 
                       tbl_joined_genelevel$glycoprot_fc_median >= fc_thresh]  <- 'UP'
tbl_cat$hit_category[tbl_joined_genelevel$glycoprot_p_fdr<fdr_thresh & 
                       tbl_joined_genelevel$glycoprot_fc_median <= 1/fc_thresh]  <- 'DOWN'


point_labels <- tbl_joined_genelevel$fly_gene_symbol
point_labels[tbl_cat$hit_category == 'NONE'] <- ''
point_labels[tbl_cat$hit_category == 'FDR_SIGNIFICANT'] <- ''

hitlist_genelevel_fcmedian <- tbl_cat %>% 
  filter(hit_category == 'UP' | hit_category == 'DOWN') %>% 
  arrange(hit_category) %>% 
  dplyr::rename(hit_category_genelevel_fcmedian = 2)

p_glycoprot_genelevel_median <- ggplot(tbl_joined_genelevel, aes(x = log2(glycoprot_fc_median), 
                                                             y = -log10(glycoprot_ttest_pval))) +
  geom_point(aes(color = tbl_cat$hit_category)) + 
  geom_text(aes(label=point_labels),hjust=1, vjust=0.35) +
  geom_vline(xintercept = log2(fc_thresh), linetype = 'dashed') +
  geom_vline(xintercept = -log2(fc_thresh), linetype = 'dashed') +
  geom_hline(yintercept = - log10(0.05), linetype = 'dashed') +
  xlim(-2, 2) +
  ylim(0, 3.5) + 
  ggtitle('Gene Level') +
  scale_fill_brewer(palette='Set1')

p_glycoprot_genelevel_median




# glycoprot site level - mean
fdr_thresh <- 0.25
fc_thresh <- 2^0.5

tbl_cat <- tbl_joined_sitelevel %>%
  dplyr::select(2) %>% 
  add_column(hit_category = 'NONE', .after = 1)
tbl_cat$hit_category[tbl_joined_sitelevel$glycoprotsitelevel_p_fdr<fdr_thresh]  <- 'FDR_SIGNIFICANT'
tbl_cat$hit_category[tbl_joined_sitelevel$glycoprotsitelevel_p_fdr<fdr_thresh & 
                       tbl_joined_sitelevel$glycoprotsitelevel_fc_mean >= fc_thresh]  <- 'UP'
tbl_cat$hit_category[tbl_joined_sitelevel$glycoprotsitelevel_p_fdr<fdr_thresh & 
                       tbl_joined_sitelevel$glycoprotsitelevel_fc_mean <= 1/fc_thresh]  <- 'DOWN'

point_labels <- tbl_joined_sitelevel$peptide_identifier
point_labels[tbl_cat$hit_category == 'NONE'] <- ''
point_labels[tbl_cat$hit_category == 'FDR_SIGNIFICANT'] <- ''


hitlist_sitelevel_fcmean <- tbl_cat %>% 
  filter(hit_category == 'UP' | hit_category == 'DOWN') %>% 
  arrange(hit_category) %>% 
  dplyr::rename(hit_category_sitelevel_fcmean = 2)

p_glycoprot_sitelevel_mean <- ggplot(tbl_joined_sitelevel, aes(x = log2(glycoprotsitelevel_fc_mean), y = -log10(glycoprotsitelevel_ttest_pval))) +
  geom_point(aes(color = tbl_cat$hit_category)) + 
  geom_text(aes(label=point_labels),hjust=1, vjust=0.35) +
  geom_vline(xintercept = log2(fc_thresh), linetype = 'dashed') +
  geom_vline(xintercept = -log2(fc_thresh), linetype = 'dashed') +
  geom_hline(yintercept = - log10(0.05), linetype = 'dashed') +
  xlim(-3, 3) +
  ylim(0, 4) + 
  ggtitle('site level') +
  scale_fill_brewer(palette='Set1')

p_glycoprot_sitelevel_mean


# glycoprot site level - median
fdr_thresh <- 0.25
fc_thresh <- 2^0.5

tbl_cat <- tbl_joined_sitelevel %>%
  dplyr::select(2) %>% 
  add_column(hit_category = 'NONE', .after = 1)
tbl_cat$hit_category[tbl_joined_sitelevel$glycoprotsitelevel_p_fdr<fdr_thresh]  <- 'FDR_SIGNIFICANT'
tbl_cat$hit_category[tbl_joined_sitelevel$glycoprotsitelevel_p_fdr<fdr_thresh & 
                       tbl_joined_sitelevel$glycoprotsitelevel_fc_median >= fc_thresh]  <- 'UP'
tbl_cat$hit_category[tbl_joined_sitelevel$glycoprotsitelevel_p_fdr<fdr_thresh & 
                       tbl_joined_sitelevel$glycoprotsitelevel_fc_median <= 1/fc_thresh]  <- 'DOWN'

point_labels <- tbl_joined_sitelevel$peptide_identifier
point_labels[tbl_cat$hit_category == 'NONE'] <- ''
point_labels[tbl_cat$hit_category == 'FDR_SIGNIFICANT'] <- ''

hitlist_sitelevel_fcmedian <- tbl_cat %>% 
  filter(hit_category == 'UP' | hit_category == 'DOWN') %>% 
  arrange(hit_category) %>% 
  dplyr::rename(hit_category_sitelevel_fcmedian = 2)

p_glycoprot_sitelevel_median <- ggplot(tbl_joined_sitelevel, aes(x = log2(glycoprotsitelevel_fc_median), y = -log10(glycoprotsitelevel_ttest_pval))) +
  geom_point(aes(color = tbl_cat$hit_category)) + 
  geom_text(aes(label=point_labels),hjust=1, vjust=0.35) +
  geom_vline(xintercept = log2(fc_thresh), linetype = 'dashed') +
  geom_vline(xintercept = -log2(fc_thresh), linetype = 'dashed') +
  geom_hline(yintercept = - log10(0.05), linetype = 'dashed') +
  xlim(-3, 3) +
  ylim(0, 4) + 
  ggtitle('site level') +
  scale_fill_brewer(palette='Set1')

p_glycoprot_sitelevel_median


## HITLIST
full_hitlist_gene <- full_join(hitlist_genelevel_fcmean, hitlist_genelevel_fcmedian)
full_hitlist_site <- full_join(hitlist_sitelevel_fcmean, hitlist_sitelevel_fcmedian)
full_hitlist_site <- full_hitlist_site %>% 
  separate(peptide_identifier, into = c('fly_gene_symbol','annotated_sequence'), remove = FALSE, sep = '_')

full_hitlist <- full_join(full_hitlist_site, full_hitlist_gene, by = 'fly_gene_symbol')

load('~/Dropbox/Broad-Research/ChemBio/R_Drpbx/future/tbl_diopt.Rdata')
tbl_diopt <- tbl_diopt %>% filter(searched_species == 'Fly')
tbl_diopt <- tbl_diopt %>%
  group_by(matched_gene) %>%
  summarize_all(list(~toString(unique(.)))) %>% 
  dplyr::rename(fly_gene_symbol = matched_gene)

full_hitlist <- left_join(full_hitlist, tbl_diopt, by = 'fly_gene_symbol') %>% 
  arrange(hit_category_sitelevel_fcmedian) %>% 
  unite(combined_hit_categories, hit_category_sitelevel_fcmean, 
        hit_category_sitelevel_fcmedian, 
        hit_category_genelevel_fcmean, 
        hit_category_genelevel_fcmedian,
        remove = FALSE, sep = '_')  %>%
  add_column(sortable_category_score = NA, .after = 'combined_hit_categories')
  

full_hitlist$sortable_category_score <- str_count(full_hitlist$combined_hit_categories, 'DOWN') - str_count(full_hitlist$combined_hit_categories, 'UP')

full_hitlist <- full_hitlist %>% 
  arrange(desc(sortable_category_score)) %>% 
  dplyr::select(-c(combined_hit_categories,sortable_category_score))

write_tsv(full_hitlist, '~/Desktop/full_hitlist.tsv')

full_hitlist_withrawdata <- left_join(left_join(full_hitlist, tbl_joined_genelevel, by='fly_gene_symbol'), 
                                  tbl_joined_sitelevel, by = 'peptide_identifier')
## CLEANUP

rm(list=setdiff(ls(), c('tbl_joined_genelevel','tbl_joined_sitelevel','full_hitlist')))

write_tsv(full_hitlist_withrawdata, 'full_hitlist_withrawdata.tsv')
write_tsv(full_hitlist, 'full_hitlist.tsv')
