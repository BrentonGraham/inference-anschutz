# Author: Brenton Graham
# Last Edited: 10/26/2022
# Description: Script to select and pre-process relevant data for longitudinal analysis project 

# Import packages
require(tidyverse)
require(phyloseq)
require(microbiome)

# Import full multiomics set
multiomics_set <- paste("/Users/bgraham/Research/883/Data/Merged Data/883MergedData.20221026.csv", sep="/") %>% 
  read_delim(delim=",")

# Enforce inclusion criteria
# 3 samples were each run on 2 separate proteomics batches to control for batching effect. 
# Since proteomics data is not used in this analysis, dropping these replicates results in 3 dropped observations.
# Only sputum samples were kept. This will be necessary if we want to bring in the 16S data.
inclusion_data <- multiomics_set %>%
  dplyr::select(-proteomics_batch, -contains("seq.")) %>% # Remove proteomics data
  unique() %>%                                     # 6 observations lost here
  filter(sample_type_16S == "Sputum") %>%          # Remove saliva samples
  mutate_at(vars(psa_muc:burk), as.numeric) %>%
  mutate(
    pseudo_aer = psa_muc + psa_nonmuc,
    staph_aur = staph + mrsa,
    achromobacter = axylos + achromo
  ) %>%
  mutate(bmi = weight / ((height * 0.01)^2)) %>%
  dplyr::select(sid, time, sid_first, ex_num, gender, genotype, age, bmi, pe_pastyr, 
                hosp_pastyr, admit_novirus, negative, pseudo_aer, staph_aur, 
                achromobacter, hflu, smalto, burk, pe_pastyr, fev1_pred, pes_total, 
                crp, elastase, contains("Bacteria/")) %>%
  # Binary variables
  mutate(pseudo_aer = ifelse(is.na(pseudo_aer), NA, ifelse(pseudo_aer == 0, 0, 1)),
         pseudo_aer = as.numeric(ifelse(negative == 1, 0, pseudo_aer)),
         staph_aur = ifelse(is.na(staph_aur), NA, ifelse(staph_aur == 0, 0, 1)),
         staph_aur = as.numeric(ifelse(negative == 1, 0, staph_aur)),
         achromobacter = ifelse(is.na(achromobacter), NA, ifelse(achromobacter == 0, 0, 1)),
         achromobacter = as.numeric(ifelse(negative == 1, 0, achromobacter)),
         hflu = as.numeric(ifelse(negative == 1, 0, ifelse(hflu == "NO RESULTS", NA, hflu))),
         smalto = as.numeric(ifelse(negative == 1, 0, ifelse(smalto == "NO RESULTS", NA, smalto))),
         burk = as.numeric(ifelse(negative == 1, 0, burk)),
         admit_novirus = ifelse(admit_novirus == "NOT DONE", NA, admit_novirus),
         negative = ifelse(negative == "NO RESULTS", NA, negative)) %>%
  arrange(sid_first, sid, time)
# Number of IDs in the unique set: 46
inclusion_data %>% dplyr::select(sid) %>% unique() %>% nrow()

# Select samples with data at all three time points
eligible_sid <- inclusion_data %>%
  dplyr::select(sid) %>%
  group_by(sid) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n == 3) %>%
  dplyr::select(sid) %>%
  pull()
# Number of IDs with three time points: 23
length(eligible_sid)

# Baseline measure of pathogen count and co-infection
baseline_pathogens <- inclusion_data %>%
  filter(sid %in% eligible_sid) %>%
  filter(time == 1) %>%
  mutate(
    admit.cf_pathogens = pseudo_aer + staph_aur + achromobacter + hflu + smalto + burk,
    admit.cf_bin = ifelse(admit.cf_pathogens == 0, 0, 1),
    admit.virus = 1 - as.numeric(admit_novirus), # Negate no-virus to virus
    admit.coinfect = as.numeric(admit.virus) + as.numeric(admit.cf_bin),
    admit.coinfect = factor(case_when(
      admit.coinfect <= 1 ~ 0,
      admit.coinfect > 1 ~ 1
    ))
  ) %>%
  dplyr::select(sid, admit.virus, admit.cf_pathogens, admit.coinfect) %>%
  arrange(sid)

# 16S microbiome data
microbiome_data <- inclusion_data %>%
  filter(sid %in% eligible_sid) %>%
  dplyr::select(sid_first, sid, time, contains("Bacteria/")) %>%
  mutate(total_reads = dplyr::select(., contains("Bacteria/")) %>% rowSums(na.rm = TRUE)) 

# CLR transformation
otu_table <- microbiome_data %>% dplyr::select(contains("Bacteria")) %>%
  t() %>% otu_table(taxa_are_rows=T)
metadata <- microbiome_data %>% dplyr::select(-contains("Bacteria"))
otu_names <- otu_table %>% as.data.frame() %>% rownames()
tax_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu_tax_split <- sapply(otu_names, FUN=function(x) strsplit(x, "/"))
tax_table <- plyr::ldply(otu_tax_split, rbind)[-1] %>%
  set_rownames(otu_names) %>%
  set_colnames(tax_levels) %>%
  as.matrix()

# Phyloseq object
physeq <- phyloseq(otu_table(otu_table), tax_table(tax_table), sample_data(metadata))
microbiomeCLR <- physeq %>%
  microbiome::transform("clr")

# Final microbiome set with CLR transform
CLR <- abundances(microbiomeCLR) %>% t() %>% as.data.frame() %>%
  merge(physeq %>% sample_data(), by = 'row.names') %>%
  column_to_rownames("Row.names") %>%
  dplyr::select(sid_first, sid, time, 
                contains("Pseudomonas aeruginosa"), contains("Staphylococcus aureus"),
                contains("Achromobacter"), contains("Haemophilus"), contains("Stenotrophomonas"), 
                contains("Bacteria/Proteobacteria/Betaproteobacteria/Burkholderiales/Burkholderiaceae/Burkholderia"))

# Remove ineligible SIDs from inclusion set
longitudinal_data <- inclusion_data %>%
  dplyr::select(-contains("Bacteria/")) %>%
  filter(sid %in% eligible_sid) %>%
  merge(CLR) %>%
  merge(baseline_pathogens) %>%
  mutate(
    admit.virus = factor(admit.virus),
    admit.cf_pathogens = factor(admit.cf_pathogens),
    admit.coinfect = factor(admit.coinfect)
  ) %>%
  relocate(contains("admit."), .before = pe_pastyr) %>%
  dplyr::select(-admit_novirus, -negative) %>%
  # Rename culture variables
  rename(pseudo_aer.cult = pseudo_aer, staph_aur.cult = staph_aur, 
         achromobacter.cult = achromobacter, hflu.cult = hflu, 
         smalto.cult = smalto, burk.cult = burk) %>%
  # Rename 16S variables
  rename(
    pseudo_aer.seq = "Bacteria/Proteobacteria/Gammaproteobacteria/Pseudomonadales/Pseudomonadaceae/Pseudomonas/Pseudomonas aeruginosa",
    staph_aur.seq = "Bacteria/Firmicutes/Bacilli/Bacillales/Staphylococcaceae/Staphylococcus/Staphylococcus aureus",
    achromobacter.seq = "Bacteria/Proteobacteria/Betaproteobacteria/Burkholderiales/Alcaligenaceae/Achromobacter",
    hflu.seq = "Bacteria/Proteobacteria/Gammaproteobacteria/Pasteurellales/Pasteurellaceae/Haemophilus",
    smalto.seq = "Bacteria/Proteobacteria/Gammaproteobacteria/Xanthomonadales/Xanthomonadaceae/Stenotrophomonas",
    burk.seq = "Bacteria/Proteobacteria/Betaproteobacteria/Burkholderiales/Burkholderiaceae/Burkholderia"
  ) %>%
  arrange(sid_first, sid, time) %>%
  dplyr::select(sid, time, sid_first, ex_num, gender, genotype, age, contains("cult"),
                contains("seq"), contains("admit"), everything())

# Number of unique subjects: 21 (2 subjects with 2 PEx)
longitudinal_data %>% dplyr::select(sid_first) %>% unique() %>% nrow()

# Output data set
longitudinal_data %>% write.table("883LongitudinalData.csv", sep = ",", row.names = F)
