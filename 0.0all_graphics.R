# Metadata analysis

# Load libraries
library(tidyverse)
library(ggplot2)
library(corrplot)
library(vroom)

# Load data
clinical <- vroom("/datos/rosmap/single_cell/metadata/ROSMAP_clinical.csv")
biospecimen <- vroom("/datos/rosmap/single_cell/metadata/ROSMAP_biospecimen_metadata.csv")
assay <- vroom("/datos/rosmap/single_cell/metadata/ROSMAP_assay_scrnaSeq_metadata.csv")

# Filter biospecimen to keep only "single nucleus" samples
biospecimen <- biospecimen %>%
  filter(nucleicAcidSource == "single nucleus", exclude == "FALSE")

# Filter clinical data to keep only selected individual IDs
individuals <- biospecimen$individualID
clinical <- clinical %>%
  filter(individualID %in% individuals)

# Save filtered individual IDs
write.csv(clinical$individualID, "filtered_individuals.csv", row.names = FALSE)

# Covariate tables
print(table(clinical$cogdx))
print(table(clinical$braaksc))
print(table(clinical$ceradsc))
print(table(clinical$educ))

# Bar plots
# Cogdx
ggplot(clinical, aes(x = factor(cogdx))) +
  geom_bar(fill = "steelblue", alpha = 0.7) +
  labs(title = "Cognitive Diagnosis (Cogdx)", x = "Category", y = "Frequency") +
  theme_minimal()

# Braak Score
freq_braaksc <- clinical %>% count(braaksc)
ggplot(freq_braaksc, aes(x = factor(braaksc), y = n)) +
  geom_col(fill = "orange", alpha = 0.7, width = 0.7) +
  labs(title = "Braak Score", x = "Category", y = "Frequency") +
  theme_minimal()

# CERAD Score
freq_ceradsc <- clinical %>% count(ceradsc)
ggplot(freq_ceradsc, aes(x = factor(ceradsc), y = n)) +
  geom_col(fill = "yellow", alpha = 0.7, width = 0.7) +
  labs(title = "CERAD Score", x = "Category", y = "Frequency") +
  theme_minimal()

# Sex distribution
freq_sex <- clinical %>% count(msex)
ggplot(freq_sex, aes(x = factor(msex), y = n)) +
  geom_col(fill = "purple", alpha = 0.7) +
  labs(title = "Sex Distribution", x = "Sex", y = "Frequency") +
  theme_minimal()

# APOE Genotype
freq_apoe <- clinical %>% count(apoe_genotype)
ggplot(freq_apoe, aes(x = factor(apoe_genotype), y = n)) +
  geom_col(fill = "brown1", alpha = 0.7) +
  labs(title = "APOE Genotype", x = "APOE Genotype", y = "Frequency") +
  theme_minimal()

# Correlation plot
datos_cor <- clinical %>%
  select(cogdx, braaksc, apoe_genotype, spanish, educ, msex, ceradsc) %>%
  mutate_all(as.numeric)

cor_matrix <- cor(datos_cor, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "square", type = "upper",
         tl.col = "red2", tl.srt = 45,
         col = colorRampPalette(c("red4", "white", "#27408B"))(200),
         cl.pos = "b", zlim = c(-1, 1))

# Stratify data
metadata_ROSMAP <- clinical %>%
  mutate(is_AD = case_when(
    cogdx == 1 & (braaksc != 0 & ceradsc %in% c(1, 2)) ~ "AD-NC_ASYM",
    cogdx == 1 & ceradsc %in% c(3, 4) ~ "control",
    cogdx %in% c(4, 5) & ceradsc == 1 ~ "AD-NC_SYM",
    cogdx %in% c(2, 3) ~ "MCI",
    TRUE ~ NA_character_
  ))

# save stratified metadata
#write.csv(metadata_ROSMAP, "metadata_ROSMAP_stratified.csv", row.names = FALSE)
