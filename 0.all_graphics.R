#metadata analysis

#Libraries --- ---

pacman::p_load("tidyverse", 
               "ggplot2")

#Get data --- ---

clinical <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_clinical.csv")

biospecimen <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_biospecimen_metadata.csv")

assay <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_assay_scrnaSeq_metadata.csv")

Exp_1_atlas <- vroom::vroom("/datos/rosmap/single_cell/metadata/Experiment1/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv")

View(biospecimen)


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



#424 IDs 
filtered_individuals <- clinical$individualID

print(filtered_individuals)
write.csv(filtered_individuals, "filtered_individuals.csv", row.names = FALSE)



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
  mutate(age_death = as.numeric(gsub("\\+", "", age_death))) %>%  # gsup 
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


#Data stratification
metadata_ROSMAP <- clinical %>%
  mutate(is_AD = case_when(
    cogdx == 1 & (braaksc != 0 & (ceradsc == 1 | ceradsc ==2)) ~ "AD-NC_ASYM",
    cogdx == 1  & (ceradsc == 4 | ceradsc == 3 ) ~ "control",
    (cogdx %in% c(4, 5) & ceradsc == 1) ~ "r AD-NC_SYM",
    cogdx %in% c(2, 3) ~ "MCI",
    TRUE ~ NA_character_
  ))

































# 
ggplot(clinical, aes(x = factor(cogdx))) +
  geom_bar(fill = "steelblue",  alpha = 0.7) +
  labs(title = "   Cogdx", x = "Category", y = "Frecuency")+
  theme_minimal()
table(clinical$cogdx)

write.csv(freq_table, "cogdx_frequencies.csv", row.names = FALSE)
View(freq_table)  

#Brack score 
ggplot(freq_braaksc, aes(x = factor(braaksc), y = n)) +
  geom_col(fill = "orange", alpha = 0.7, width = 0.7) +
  scale_x_discrete(limits = as.character(0:6)) + 
  labs(title = "Braak Score", x = "Category", y = "Frequency") +
  theme_minimal()

table(clinical$braaksc)

freq_braaksc <- clinical %>%
  count(braaksc) %>%
  arrange(desc(n))
write.csv(freq_braaksc, "braaksc_frequencies.csv", row.names = FALSE)
View(freq_braaksc)

#CERAD score
print(freq_ceradsc)

ggplot(freq_ceradsc, aes(x = factor(ceradsc), y = n)) +
  geom_col(fill = "yellow", alpha = 0.7, width = 0.7) +  
  scale_x_discrete(limits = as.character(1:4)) + 
  labs(title = "CERAD Score", x = "Category", y = "Frequency") +
  theme_minimal()


freq_ceradsc <- clinical %>%
  count(ceradsc) %>%
  arrange(desc(n))
write.csv(freq_ceradsc, "ceradsc_frequencies.csv", row.names = FALSE)
View(freq_ceradsc)



#

ggplot(clinical, aes(x = factor(dcfdx_lv))) +
  geom_bar(fill = "pink", alpha = 0.7) +
  labs(title = " Dcfdx", x = "Category", y = "Frecuency")+
  theme_minimal()

freq_dcfdx_lv <- clinical %>%
  count(dcfdx_lv) %>%
  arrange(desc(n))
write.csv(freq_ceradsc, "dcfdx_frequencies.csv", row.names = FALSE)
View(freq_dcfdx_lv)


# Demographic overview

#Sex 
ggplot(sex_dx_counts, aes(x = factor(msex), y = n)) +
  geom_bar(stat = "identity", fill = "purple", alpha = 0.7) +
  labs(title = "Sex distribution", x = "Sexo", y = "Frecuencia") +
  theme_minimal()

sex_dx_counts <- clinical %>%
  count(msex) %>%
  arrange(desc(n))
write.csv(sex_dx_counts, "msex_frequencies.csv", row.names = FALSE)
View(sex_dx_counts)

#Race
ggplot(race_counts, aes(x = factor(race), y = n)) +
  geom_bar(stat = "identity", fill = "magenta", alpha = 0.7) +
  labs(title = "Race", x = "Race", y = "Frecuencia") +
  theme_minimal()
  
library(dplyr)

# ROS/MAP
study_summary_by_dx <- clinical %>%
  count(Study) %>%
  arrange(desc(n))
print(study_summary_by_dx)

ggplot(study_summary_by_dx, aes(x = factor(Study), y = n)) +
  geom_bar(stat = "identity", fill = "cyan4", alpha = 0.7) +
  labs(title = "Study", x = "Study", y = "Frecuency") +
  theme_minimal()

#Apoe-genotype
apoe_summary_by_dx <- clinical %>%
  count(apoe_genotype) %>%
  arrange(desc(n))
print(apoe_summary_by_dx)

ggplot(apoe_summary_by_dx, aes(x = factor(apoe_genotype), y = n)) +
  geom_bar(stat = "identity", fill = "brown1", alpha = 0.7) +
  labs(title = " APOE Genotype", x = "APOE genotype ", y = "Frecuencia") +
  theme_minimal()

write.csv(apoe_summary_by_dx, "apoe_genotype_frequencies.csv", row.names = FALSE)
View(apoe_summary_by_dx)

# make de Corplot 
install.packages("corrplot")
install.packages("ggcorrplot")
library(corrplot)
library(ggcorrplot)
library(tidyverse)

# Filter the columns of interest 
datos_cor <- clinical %>% 
  select(cogdx, braaksc, apoe_genotype, spanish,educ, msex, ceradsc)

datos_cor_numeric <- datos_cor %>% mutate_all(as.numeric)


cor_matrix <- cor(datos_cor_numeric, use = "pairwise.complete.obs")

corrplot(cor_matrix, method = "square", type = "upper", 
         tl.col = "red2", tl.srt = 45, 
         col = colorRampPalette(c("red4", "white", "#27408B"))(200),
         cl.pos = "b",            
         zlim = c(-1, 1),         
         cl.breaks = seq(-1, 1, by = 0.2),
         cl.cex = 0.7,            
         cl.align = "c")          


#libraries  ----- 
library(corrplot)
library(tidyverse)
# make corplot with numbers

B <- clinical %>% select(cogdx, braaksc, apoe_genotype, spanish,educ, msex, ceradsc)
B <- B %>% drop_na()

cor_matrix <- cor(B, use = "pairwise.complete.obs", method = "pearson")

corrplot(cor_matrix, method = "number", 
         number.cex = 0.8,  # Tamaño de los números
         col = colorRampPalette(c("red4", "white", "#27408B"))(200),  # Escala de color
         addCoef.col = "black",  # Números en negro
         cl.pos = "b",  # Barra de calor en la parte inferior
         tl.col = "red2", tl.srt = 45,  # Ajustes en etiquetas
         addgrid.col = "gray90")  # Líneas más suaves



# # data stratification

metadata_ROSMAP <- clinical %>%
mutate(is_AD = case_when(
  cogdx == 1 & (braaksc != 0 & (ceradsc == 1 | ceradsc ==2)) ~ "AD-NC_ASYM",
  cogdx == 1  & (ceradsc == 4 | ceradsc == 3 ) ~ "control",
  (cogdx %in% c(4, 5) & ceradsc == 1) ~ "r AD-NC_SYM",
  cogdx %in% c(2, 3) ~ "MCI",
  TRUE ~ NA_character_
))