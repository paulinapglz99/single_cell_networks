#metadata analysis

#Libraries --- ---

pacman::p_load("tidyverse", 
               "ggplot2")

#Get data 

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

# Crear una nueva variable para la clasificación de la enfermedad de Alzheimer y otros estados cognitivos
clinical <- clinical %>%
  mutate(is_AD = case_when(
    cogdx == 1 & (braaksc != 0 & (ceradsc == 1 | ceradsc == 2)) ~ "AD-NC_ASYM",
    cogdx == 1 & (ceradsc == 4 | ceradsc == 3) ~ "control",
    (cogdx %in% c(4, 5) & ceradsc == 1) ~ "AD-NC_SYM",
    cogdx %in% c(2, 3) ~ "MCI",
    TRUE ~ NA_character_
  ))


table(clinical$is_AD)

write.csv(clinical, "clinical_stratified.csv", row.names = FALSE)


#Demographic and genetic summary table preparation
clinical <- clinical %>%
  mutate(age_death = as.numeric(gsub("\\+", "", age_death)))

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
write.csv(tabla_final, "tabla_clinica_resumen.csv", row.names = FALSE)

library(gt)


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

# ANOVA analysis

# Filter data with valid stratification
clinical_test <- clinical %>%
  filter(!is.na(is_AD)) 

# ANOVA for age at death
anova_age <- aov(age_death ~ is_AD, data = clinical_test)
summary(anova_age)

# ANOVA for education
anova_educ <- aov(educ ~ is_AD, data = clinical_test)
summary(anova_educ)

# If ANOVA shows significant differences, perform Tukey's HSD test to identify specific group differences
turkey <- TukeyHSD(anova_age)
TukeyHSD(anova_educ)

turkey.df <- as.data.frame(turkey$is_AD)
turkey.df$sig <- ifelse(turkey.df$`p adj` > 0.05, "Not sig", "Sig")
print(turkey.df)
turkey.df
# ANOVA for education
anova_educ <- aov(educ ~ is_AD, data = clinical_test)
summary(anova_educ)

# Tukey HSD for education
tukey_educ <- TukeyHSD(anova_educ)

# Display results in a table
tukey.df <- as.data.frame(tukey_educ$is_AD)
tukey.df$sig <- ifelse(tukey.df$`p adj` > 0.05, "Not sig", "Sig")

# Print results
print(tukey.df)

# --------------------------------------
# Load necessary libraries
pacman::p_load(tidyverse)

# Clean and prepare data
clinical_test <- clinical %>%
  filter(!is.na(is_AD)) %>%
  mutate(
    age_death = as.numeric(gsub("\\+", "", age_death)),
    msex = as.numeric(msex),
    Study_num = ifelse(Study == "ROS", 1, 0)
  )

# Function to run ANOVA + Tukey and return a table with significance column
get_tukey_table <- function(variable, var_name = "Variable") {
  model <- aov(variable ~ is_AD, data = clinical_test)
  tukey <- TukeyHSD(model)
  df <- as.data.frame(tukey$is_AD)
  df$Comparison <- rownames(df)
  df$Significant <- ifelse(df$`p adj` < 0.05, "Yes", "No")
  df$Variable <- var_name
  rownames(df) <- NULL
  df <- df[, c("Variable", "Comparison", "diff", "lwr", "upr", "p adj", "Significant")]
  return(df)
}

# Run for each variable
table_age   <- get_tukey_table(clinical_test$age_death, "age_death")
table_educ  <- get_tukey_table(clinical_test$educ, "educ")
table_sex   <- get_tukey_table(clinical_test$msex, "Sex")
table_study <- get_tukey_table(clinical_test$Study_num, "Study (0=MAP, 1=ROS)")

# Combine all tables into one
final_tukey_table <- bind_rows(table_age, table_educ, table_sex, table_study)

# View table
View(final_tukey_table)

# (Optional) Save as CSV
write.csv(final_tukey_table, "tukey_significance_results.csv", row.names = FALSE)

# --------------------------------------
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

# Merge all matrices by "Comparison"
significance_matrix <- reduce(
  list(sig_age, sig_educ, sig_sex, sig_study),
  full_join,
  by = "Comparison"
)

# View the combined matrix
View(significance_matrix)

# (Optional) Save as CSV
write.csv(significance_matrix, "significance_matrix.csv", row.names = FALSE)

# --------------------------------------

create_tukey_heatmap <- function(variable, variable_label) {
  model <- aov(variable ~ is_AD, data = clinical_test)
  tukey <- TukeyHSD(model)
  df <- as.data.frame(tukey$is_AD)
  
  # Create a column with the name of the comparisons
  df$Comparison <- rownames(df)
  
  # Split comparison into two groups - using separation from the last dash
  df <- df %>%
    separate(Comparison, into = c("Group1", "Group2"), sep = "(?<=\\w)-(?!.*-)", remove = FALSE) %>%
    mutate(across(c(Group1, Group2), str_trim)) %>%
    mutate(Significant = ifelse(`p adj` < 0.05, "Yes", "No"))
  
  # Create all possible combinations between groups
  groups <- sort(unique(c(df$Group1, df$Group2)))
  combinations <- expand.grid(Group1 = groups, Group2 = groups, stringsAsFactors = FALSE)
  
  # Join combinations with results
  matrix_long <- combinations %>%
    left_join(df, by = c("Group1", "Group2")) %>%
    mutate(Significant = case_when(
      Group1 == Group2 ~ "—",
      is.na(Significant) ~ "No",
      TRUE ~ Significant
    ))
  
  # Plot heatmap
  ggplot(matrix_long, aes(x = Group1, y = Group2, fill = Significant)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("Yes" = "tomato", "No" = "skyblue", "—" = "grey90")) +
    coord_fixed() +
    labs(
      title = paste("Significance between groups -", variable_label),
      x = "Group 1", y = "Group 2", fill = "Significant?"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

  
  
  # Generalized function to create heatmap from Tukey results
  
  library(tidyverse)
  
  create_tukey_heatmap <- function(variable, variable_label) {
    model <- aov(variable ~ is_AD, data = clinical_test)
    tukey <- TukeyHSD(model)
    df <- as.data.frame(tukey$is_AD)
    df$Comparison <- rownames(df)
    df$Significant <- ifelse(df$`p adj` < 0.05, "Yes", "No")
    
    # Split comparison into group pairs
    df <- df %>%
      separate(Comparison, into = c("Group1", "Group2"), sep = "(?<=\\w)-(?=\\w)") %>%
      mutate(across(c(Group1, Group2), str_trim))
    
    # All possible combinations
    groups <- sort(unique(c(df$Group1, df$Group2)))
    separate(Comparison, into = c("Group1", "Group2"), sep = "(?<=\\w)-(?!.*-)", remove = FALSE)
    
    #combinations <- expand.grid(Group1 = groups, Group2 = groups, stringsAsFactors = FALSE)
    
    # Join with Tukey results
    matrix_long <- combinations %>%
      left_join(df, by = c("Group1", "Group2")) %>%
      mutate(Significant = case_when(
        Group1 == Group2 ~ "—",         # diagonal
        is.na(Significant) ~ "No",      # fill empty with No
        TRUE ~ Significant
      ))
    
    # Plot heatmap
    ggplot(matrix_long, aes(x = Group1, y = Group2, fill = Significant)) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("Yes" = "tomato", "No" = "skyblue", "—" = "grey90")) +
      coord_fixed() +
      labs(
        title = paste("Significance between groups -", variable_label),
        x = "Group 1", y = "Group 2", fill = "Significant?"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Example: run for age
  create_tukey_heatmap(clinical_test$age_death, "Age at death")
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  # Heatmap visualization using geom_tile
  
  library(ggplot2)
  
  # Convert matrix to long format
  long_matrix <- matriz_edad %>%
    rownames_to_column("Group1") %>%
    pivot_longer(-Group1, names_to = "Group2", values_to = "Significant")
  
  # Heatmap plot
  ggplot(long_matrix, aes(x = Group1, y = Group2, fill = Significant)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = c("Yes" = "tomato", "No" = "skyblue", "—" = "grey90")
    ) +
    coord_fixed() +  # square tiles
    labs(
      title = "Significance between groups - Age at death",
      x = "Group 1",
      y = "Group 2",
      fill = "Significant?"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # --------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #ANOVA 
  # Filter data with valid estratification
  clinical_test <- clinical %>%
    filter(!is.na(is_AD)) 
  
  # ANOVA for age 
  anova_age <- aov(age_death ~ is_AD, data = clinical_test)
  summary(anova_age)
  
  # ANOVA for edu
  anova_educ <- aov(educ ~ is_AD, data = clinical_test)
  summary(anova_educ)
  
  ## If ANOVA shows significant differences, perform Tukey's HSD test to identify specific group differences
  turkey <- TukeyHSD(anova_age)
  TukeyHSD(anova_educ)
  turkey.df <- as.data.frame(turkey$is_AD)
  turkey.df$sig <- ifelse(turkey.df$`p adj` > 0.05, "No sig", "Sig")
  print(tukey.df)
  tukey.df
  
  
  
  
  # ANOVA para educación
  anova_educ <- aov(educ ~ is_AD, data = clinical_test)
  summary(anova_educ)
  
  # Tukey para educación
  tukey_educ <- TukeyHSD(anova_educ)
  
  # Visualizar resultados en tabla
  tukey.df <- as.data.frame(tukey_educ$is_AD)
  tukey.df$sig <- ifelse(tukey.df$`p adj` > 0.05, "No sig", "Sig")
  
  # Imprimir
  print(tukey.df)
  
  ###
  # Cargar librerías necesarias
  pacman::p_load(tidyverse)
  
  # Limpiar y preparar los datos
  clinical_test <- clinical %>%
    filter(!is.na(is_AD)) %>%
    mutate(
      age_death = as.numeric(gsub("\\+", "", age_death)),
      msex = as.numeric(msex),
      Study_num = ifelse(Study == "ROS", 1, 0)
    )
  
  # Función para ejecutar ANOVA + Tukey y devolver tabla con columna de significancia
  get_tukey_table <- function(variable, var_name = "Variable") {
    model <- aov(variable ~ is_AD, data = clinical_test)
    tukey <- TukeyHSD(model)
    df <- as.data.frame(tukey$is_AD)
    df$Comparación <- rownames(df)
    df$Significativo <- ifelse(df$`p adj` < 0.05, "Sí", "No")
    df$Variable <- var_name
    rownames(df) <- NULL
    df <- df[, c("Variable", "Comparación", "diff", "lwr", "upr", "p adj", "Significativo")]
    return(df)
  }
  
  # Ejecutar para cada variable
  tabla_age   <- get_tukey_table(clinical_test$age_death, "age_death")
  tabla_educ  <- get_tukey_table(clinical_test$educ, "educ")
  tabla_sex   <- get_tukey_table(clinical_test$msex, "Sex")
  tabla_study <- get_tukey_table(clinical_test$Study_num, "Study (0=MAP, 1=ROS)")
  
  # Unir todas las tablas en una sola
  tabla_final_tukey <- bind_rows(tabla_age, tabla_educ, tabla_sex, tabla_study)
  
  # Visualizar tabla
  View(tabla_final_tukey)
  
  # (Opcional) Guardar como CSV
  write.csv(tabla_final_tukey, "tukey_resultados_significancia.csv", row.names = FALSE)
  
  
  
  
  # Función para sacar matriz de significancia
  get_significance_matrix <- function(var, var_label) {
    model <- aov(var ~ is_AD, data = clinical_test)
    tukey <- TukeyHSD(model)
    df <- as.data.frame(tukey$is_AD)
    df$Comparación <- rownames(df)
    df[[var_label]] <- ifelse(df$`p adj` < 0.05, "Sí", "No")
    df %>% select(Comparación, !!sym(var_label))
  }
  
  # Obtener matrices individuales
  sig_age   <- get_significance_matrix(clinical_test$age_death, "Edad")
  sig_educ  <- get_significance_matrix(clinical_test$educ, "Educación")
  sig_sex   <- get_significance_matrix(clinical_test$msex, "Sexo")
  sig_study <- get_significance_matrix(clinical_test$Study_num, "Estudio")
  
  # Unir todas por "Comparación"
  matriz_significancia <- reduce(
    list(sig_age, sig_educ, sig_sex, sig_study),
    full_join,
    by = "Comparación"
  )
  
  # Visualizar como tabla
  View(matriz_significancia)
  
  # (Opcional) Guardar como CSV
  write.csv(matriz_significancia, "matriz_significancia.csv", row.names = FALSE)
  
  
  library(ggplot2)
  
  # Convertir la matriz a formato largo
  matriz_edad_long <- matriz_edad %>%
    rownames_to_column("Grupo1") %>%
    pivot_longer(-Grupo1, names_to = "Grupo2", values_to = "Significativo")
  
  # Gráfico tipo heatmap
  ggplot(matriz_edad_long, aes(x = Grupo1, y = Grupo2, fill = Significativo)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = c("Sí" = "tomato", "No" = "skyblue", "—" = "grey90")
    ) +
    coord_fixed() +  # para que los tiles sean cuadrados
    labs(
      title = "Significancia entre grupos - Edad al morir",
      x = "Grupo 1",
      y = "Grupo 2",
      fill = "¿Significativo?"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  
  
  
  
  ####
  
  library(tidyverse)
  
  # Función para construir heatmap de significancia entre grupos
  graficar_tukey_matrix <- function(variable, variable_label) {
    modelo <- aov(variable ~ is_AD, data = clinical_test)
    tukey <- TukeyHSD(modelo)
    df <- as.data.frame(tukey$is_AD)
    df$Comparación <- rownames(df)
    df$Significativo <- ifelse(df$`p adj` < 0.05, "Sí", "No")
    
    # Separar los grupos correctamente
    df <- df %>%
      separate(Comparación, into = c("Grupo1", "Grupo2"), sep = "(?<=\\w)-(?=\\w)") %>%
      mutate(across(c(Grupo1, Grupo2), str_trim))
    
    # Crear todas las combinaciones posibles
    grupos <- sort(unique(c(df$Grupo1, df$Grupo2)))
    combinaciones <- expand.grid(Grupo1 = grupos, Grupo2 = grupos, stringsAsFactors = FALSE)
    
    # Unir con los datos del Tukey
    matriz_long <- combinaciones %>%
      left_join(df, by = c("Grupo1", "Grupo2")) %>%
      mutate(Significativo = case_when(
        Grupo1 == Grupo2 ~ "—",                  # Diagonal
        is.na(Significativo) ~ "No",             # Si no hay valor, no significativo
        TRUE ~ Significativo
      ))
    
    # Graficar con geom_tile
    ggplot(matriz_long, aes(x = Grupo1, y = Grupo2, fill = Significativo)) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("Sí" = "tomato", "No" = "skyblue", "—" = "grey90")) +
      coord_fixed() +
      labs(
        title = paste("Significancia entre grupos -", variable_label),
        x = "Grupo 1", y = "Grupo 2", fill = "¿Significativo?"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Ejecutar para edad, por ejemplo
  graficar_tukey_matrix(clinical_test$age_death, "Edad al morir")
  
  



