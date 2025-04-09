#metadata analysis

#Libraries --- ---
pacman::p_load("tidyverse", 
               "ggplot2", 
               "gt")
#Get data 
clinical <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_clinical.csv")
biospecimen <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_biospecimen_metadata.csv")
assay <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_assay_scrnaSeq_metadata.csv")
Exp_1_atlas <- vroom::vroom("/datos/rosmap/single_cell/metadata/Experiment1/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv")
demultiplex <- vroom::vroom("/datos/rosmap/single_cell/metadata/ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv")

#Filter to obtain only single nucleus assays
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

#Demographics
individuals <- biospecimen$individualID
dim(clinical)

clinical <- clinical %>% 
  filter(clinical$individualID %in% individuals)
dim(clinical)
 
filtered_individuals <- clinical$individualID

length(filtered_individuals)
#write.csv(filtered_individuals, "filtered_individuals.csv", row.names = FALSE)

#Covariate table 

table(clinical$cogdx)
table(clinical$braaksc)
table(clinical$ceradsc)
table(clinical$educ)

# Create Alzheimer's diagnosis classification variable
clinical <- clinical %>%
  mutate(is_AD = case_when(
    cogdx == 1 & (braaksc != 0 & (ceradsc == 1 | ceradsc == 2)) ~ "AD-NC_ASYM",
    cogdx == 1 & (ceradsc == 4 | ceradsc == 3) ~ "control",
    (cogdx %in% c(4, 5) & ceradsc == 1) ~ "AD-NC_SYM",
    cogdx %in% c(2, 3) ~ "MCI",
    TRUE ~ NA_character_
  ))

clinical <- clinical %>%
  mutate(age_death = as.numeric(gsub("\\+", "", age_death)))

table(clinical$is_AD)
#write.csv(clinical, "clinical_stratified.csv", row.names = FALSE)

# Create demographic summary table grouped by Alzheimer's diagnosis (is_AD)
tabla_resumen <- clinical %>%
  group_by(is_AD) %>%
  summarise(
    n = as.character(n()),
    Age = paste0(round(mean(age_death, na.rm = TRUE), 1), " ± ", round(sd(age_death, na.rm = TRUE), 2)),
    Education = paste0(round(mean(educ, na.rm = TRUE), 1), " ± ", round(sd(educ, na.rm = TRUE), 2)),
    Males = as.character(sum(msex == 1, na.rm = TRUE)),
    Females = as.character(sum(msex == 0, na.rm = TRUE))
  ) %>%
  rename(Categoria = is_AD)

# Count APOE genotypes per group and convert counts to character for merging
apoe_tabla <- clinical %>%
  mutate(apoe_genotype = ifelse(is.na(apoe_genotype), "Unknown", apoe_genotype)) %>%
  group_by(is_AD, apoe_genotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = is_AD, values_from = Count, values_fill = 0) %>%
  rename(Genotype = apoe_genotype) %>%
  mutate(across(-Genotype, as.character))  # Convert all but Genotype column to character
# Transpose summary table so group categories become columns
tabla_resumen_transpuesta <- tabla_resumen %>%
  pivot_longer(cols = -Categoria, names_to = "Variable", values_to = "Valor") %>%
  pivot_wider(names_from = Categoria, values_from = Valor)

# Combine demographic summary and APOE genotype counts into one final table
tabla_final <- bind_rows(
  tabla_resumen_transpuesta %>% rename(Genotype = Variable),
  apoe_tabla
)

print(tabla_final)
#write.csv(tabla_final, "tabla_clinica_resumen.csv", row.names = FALSE)

# tabla with gt
tabla_final_gt <- tabla_final %>%
  gt(rowname_col = "Genotype") %>%
  tab_header(
    title = md("**Table 1.** Demographic and genetic characteristics of cognitive groups")
  ) %>%

  # Group rows into meaningful sections
  tab_row_group(label = "APOE genotype", rows = Genotype %in% c("22", "23", "24", "33", "34", "44", "Unknown")) %>%
  tab_row_group(label = "Sex", rows = Genotype %in% c("Males", "Females")) %>%
  tab_row_group(label = "Years of education", rows = Genotype == "Education") %>%
  tab_row_group(label = "Age", rows = Genotype == "Age") %>%
  tab_row_group(label = "Sample size (n)", rows = Genotype == "n") %>%

  cols_align(align = "center", columns = everything()) %>%

  #style
  tab_options(
    table.border.top.width = px(2),
    table.border.bottom.width = px(2),
    heading.title.font.size = 14,
    heading.title.font.weight = "bold",
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    row_group.font.size = 12
  )

tabla_final_gt

#ANOVA 
clinical_test <- clinical %>%
  filter(!is.na(is_AD)) %>%
  mutate(
    msex = as.numeric(msex),
    Study_num = ifelse(Study == "ROS", 1, 0)
  )

# Function to run ANOVA + Tukey and return results
get_tukey_table <- function(variable, var_name = "Variable"){
  model <- aov(variable ~ is_AD, data = clinical_test)
  print(summary(model))
  tukey <- TukeyHSD(model)
  
  df <- as.data.frame(tukey$is_AD)
  df$Comparison <- rownames(df)
  df$Significant <- ifelse(df$`p adj` < 0.05, "Yes", "No")
  df$Variable <- var_name
  rownames(df) <- NULL
  df <- df[, c("Variable", "Comparison", "diff", "lwr", "upr", "p adj", "Significant")]
  return(df) } 

  # Run for each variable
table_age   <- get_tukey_table(clinical_test$age_death, "age_death")
table_educ  <- get_tukey_table(clinical_test$educ, "educ")
table_sex   <- get_tukey_table(clinical_test$msex, "Sex")
table_study <- get_tukey_table(clinical_test$Study_num, "Study (0=MAP, 1=ROS)")

# Combine all tables into one
final_tukey_table <- bind_rows(table_age, table_educ, table_sex, table_study)

#write.csv(final_tukey_table, "tukey_significance_results.csv", row.names = FALSE)

# Function to create significance matrix
get_significance_matrix <- function(var, var_label) {
  model <- aov(var ~ is_AD, data = clinical_test)
  tukey <- TukeyHSD(model)
  df <- as.data.frame(tukey$is_AD)
  df$Comparison <- rownames(df)
  df[[var_label]] <- ifelse(df$`p adj` < 0.05, "Yes", "No")
  df %>% select(Comparison, !!sym(var_label))
}

# Get individual matrices
sig_age   <- get_significance_matrix(clinical_test$age_death, "Age")
sig_educ  <- get_significance_matrix(clinical_test$educ, "Education")
sig_sex   <- get_significance_matrix(clinical_test$msex, "Sex")
sig_study <- get_significance_matrix(clinical_test$Study_num, "Study")

# Merge all matrices
significance_matrix <- reduce(
  list(sig_age, sig_educ, sig_sex, sig_study),
  full_join,
  by = "Comparison"
)

#write.csv(significance_matrix, "significance_matrix.csv", row.names = FALSE)

# Prepare for heatmap
significance_matrix_long <- significance_matrix %>% 
  mutate(across(-Comparison, ~ ifelse(. == "No", 0, ifelse(. == "Yes", 1, .)))) %>% 
  mutate(across(-Comparison, as.numeric)) %>% 
  pivot_longer(cols = -Comparison, names_to = "Variable", values_to = "Sig")

# Plot
significance_matrix_plot <- ggplot(significance_matrix_long, 
                                   aes(x = Variable, y = Comparison, fill = as.factor(Sig))) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  scale_fill_manual(values = c("0" = "#4F94CD", "1" = "indianred1")) +

  theme_classic()

significance_matrix_plot


#Add data to multiplexing info --- ---



