library(dplyr)

data <- ipsc_family_variants

#extract fam3 variants
filtered_data <- data %>%
  filter(Sample %in% c("3_1_10", "3_3_1", "3_2_1", "3_4_2", "SG226"))

write.csv(filtered_data, "iPSC_variant_data_fam3.csv")

#3110_versus_321_331_342

iPSC_variant_data_3 <- iPSC_variant_data_fam3

immature_3110_versus <- im_3110_versus_321_331_342_deseq2_padj0.05

filtered_variant <- iPSC_variant_data_3 %>%
  filter(gene_id %in% immature_3110_versus$ensembl_gene_id)

write.csv(filtered_variant, "variant_fam3_3110_versus_331_321_342_ipsc.csv")

#filter vid only 3_1_10 and 3_2_1 have but not in 3_3_1 and 3_4_2

data <- variant_fam3_3110_versus_331_321_342_npc

sample_3110 <- data %>% 
  dplyr::filter(Sample == "3_1_10") %>% 
  dplyr::select(vid)

sample_321 <- data %>% 
  dplyr::filter(Sample == "3_2_1") %>% 
  dplyr::select(vid)

sample_331 <- data %>% 
  dplyr::filter(Sample == "3_3_1") %>% 
  dplyr::select(vid)

sample_342 <- data %>% 
  dplyr::filter(Sample == "3_4_2") %>% 
  dplyr::select(vid)



# Merge sample_331 and sample_342
merged_331_342 <- union(sample_331$vid, sample_342$vid)


# Find vids that meet the conditions
valid_vids <- intersect(sample_3110$vid, sample_321$vid) %>%
  setdiff(merged_331_342)


filtered_data <- data %>% dplyr::filter(vid %in% valid_vids)


write.csv(filtered_data, "variant_3110_321_immature.csv")
