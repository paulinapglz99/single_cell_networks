#metadata analysis

#Libraries --- ---

pacman::p_load("tidyverse", 
               "ggplot2")

#Get data --- ---

clinical <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_clinical.csv")

biospecimen <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_biospecimen_metadata.csv")

assay <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_assay_scrnaSeq_metadata.csv")

Exp_1_atlas <- vroom::vroom("/datos/rosmap/single_cell/metadata/Experiment1/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv")

#Filter to obtain only sc assays --- ---

dim(assay)
table(assay$assay)
table(assay$dataContributionBatch)
table(assay$platform)
table(assay$platformLocation)

#Biospecimen

biospecimen <- biospecimen %>% 
  filter(nucleicAcidSource == "single nucleus") %>% 
  filter(exclude == "FALSE")

table(biospecimen$tissue)
table(biospecimen$organ)
table(biospecimen$BrodmannArea)
table(biospecimen$nucleicAcidSource)

#Demographics --- ---

individuals <- biospecimen$individualID

dim(clinical)
clinical <- clinical %>% 
  filter(clinical$individualID %in% individuals)
dim(clinical)
View(clinical)

#Covariate table --- ---

table(clinical$cogdx)

table(clinical$braaksc)

table(clinical$ceradsc)

table(clinical$educ)

#Age by dx

age <- as.numeric(gsub("\\+", "", metadata_DLFPC_ROSMAP$age_death))
age_mean <- mean(age, na.rm = TRUE)
age_sd <- sd(age, na.rm = TRUE)

age_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  mutate(age_death = as.numeric(gsub("\\+", "", age_death))) %>% 
  group_by(dicho_NIA_reagan) %>%
  summarise(
    mean_age = mean(age_death, na.rm = TRUE),
    sd_age = sd(age_death, na.rm = TRUE)
  )

table(metadata_DLFPC_ROSMAP$dicho_NIA_reagan)

#Schooling

mean_schooling <-  metadata_DLFPC_ROSMAP %>%
  group_by(dicho_NIA_reagan) %>% 
  summarise(
    mean_educ = mean(educ, na.rm = TRUE),
    sd_educ = sd(educ, na.rm = TRUE)
  )

#Individuals by sex

sex_dx_counts <- metadata_DLFPC_ROSMAP %>%
  count(dicho_NIA_reagan, msex)  # Contar el número de individuos por diagnóstico y sexo

#Post-mortem interval

pmi_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  group_by(dicho_NIA_reagan) %>%
  summarise(
    mean_pmi = mean(pmi, na.rm = TRUE),
    sd_pmi = sd(pmi, na.rm = TRUE)
  )

#APOE genotype

apoe_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  count(dicho_NIA_reagan, apoe_genotype) 

#Race 

race_counts <- metadata_DLFPC_ROSMAP %>%
  group_by(dicho_NIA_reagan) %>% 
  mutate(race = case_when(
    race == 1 ~ "White",
    race == 2 ~ "Black or African American",
    race == 3 ~ "American Indian or Alaska Native",
    race == 4 ~ "Native Hawaiian or Other Pacific Islander",
    race == 5 ~ "Asian",
    race == 6 ~ "Other",
    race == 7 ~ "Unknown",
    TRUE ~ "Missing"
  )) %>% 
  count(race)

#Study

study_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  count(dicho_NIA_reagan, Study) 

#END