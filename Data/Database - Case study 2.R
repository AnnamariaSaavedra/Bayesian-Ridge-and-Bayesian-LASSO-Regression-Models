# 0. Set seed

rm(list = ls()); set.seed(123)

# 1. Load necessary libraries

suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(fastDummies)))

# 2. Import database

# Household and Respondent Characteristics

Home <- read_delim("~/Downloads/Características y composición del hogar.CSV", 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

Home <- Home %>%
  select(DIRECTORIO, SECUENCIA_ENCUESTA, SECUENCIA_P, ORDEN, # Primary key that uniquely identifies each observation
         P6020, # Biological sex at birth
         P6040, # Age
         P6080, # Ethnic group
         P6087, # Father's highest level of educational attainment
         P6088) %>% # Mother's highest level of educational attainment
  # Create the primary key that uniquely identifies each observation
  unite("ID", c(DIRECTORIO, SECUENCIA_ENCUESTA, SECUENCIA_P, ORDEN), sep = "", remove = TRUE) 

# Educational attainment

Education <- read_delim("~/Downloads/Educación.CSV", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)

Education <- Education %>% 
  select(DIRECTORIO, SECUENCIA_ENCUESTA, SECUENCIA_P, ORDEN, # Primary key that uniquely identifies each observation
         P8587) %>% # Highest level of educational attainment
  # Create the primary key that uniquely identifies each observation
  unite("ID", c(DIRECTORIO, SECUENCIA_ENCUESTA, SECUENCIA_P, ORDEN), sep = "", remove = TRUE)

# Labor force

Labor <- read_delim("~/Downloads/Fuerza de trabajo.CSV", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

Labor <- Labor %>% 
  select(DIRECTORIO, SECUENCIA_ENCUESTA, SECUENCIA_P, ORDEN, # Primary key that uniquely identifies each observation
         P6435, # Target population
         P6450, # Employment contract
         P8624, # Gross Labor Income
         P8626, P8626S1, # Food subsidies
         P8628, P8628S1, # Transport subsidies
         P8630, P8630S1) %>% # Family subsidies
  # Create the primary key that uniquely identifies each observation
  unite("ID", c(DIRECTORIO, SECUENCIA_ENCUESTA, SECUENCIA_P, ORDEN), sep = "", remove = TRUE)

# 3. Identify the target population

Data <- merge(x = Labor, y = Education, by = "ID", all.x = TRUE)
Data <- merge(x = Data, y = Home, by = "ID", all.x = TRUE)

Data <- Data %>%
  filter(P6435 == "1" | P6435 == "2") %>%
  # Adjust the survey response for those who do not report the information
  filter(P6087 != "10" & P6088 != "10") %>%
  filter(!is.na(P6087) & !is.na(P6088)) %>%
  mutate(P8626S1 = case_when(P8626S1 == 98 ~ 0,
                             TRUE ~ P8626S1)) %>%
  mutate(P8628S1 = case_when(P8628S1 == 98 ~ 0,
                             TRUE ~ P8628S1)) %>%
  mutate(P8630S1 = case_when(P8630S1 == 98 ~ 0,
                             TRUE ~ P8630S1))

# 3. Create the response variable

Data$Wage <- rowSums(x = Data[, c("P8624", "P8626S1", "P8628S1", "P8630S1")], na.rm = TRUE)

# Delete the survey response for those who do not report a positive labor income

Data <- subset(Data, !(Wage %in% c(0, 98, 99)))
Data <- Data %>% filter(!is.na(P8587))

Data$Wage <- log(Data$Wage)

# 4. Simple Random Sampling

mas <- sample(x = Data$ID, size = 1500, replace = FALSE)

Data <- Data %>%
  filter(ID %in% mas)

# 5. Create the matrix of explanatory variables

Data <- Data %>%
  mutate(across(c(P6020, P6080, P6087, P6088, P6450, P8587), as.factor))

# Set the reference category and perform dummy encoding

Data$P6020 <- relevel(Data$P6020, ref = "1")
Data$P6080 <- relevel(Data$P6080, ref = "6")
Data$P6087 <- relevel(Data$P6087, ref = "4")
Data$P6088 <- relevel(Data$P6088, ref = "4")
Data$P6450 <- relevel(Data$P6450, ref = "1")
Data$P8587 <- relevel(Data$P8587, ref = "5")

Data <- Data %>%
  select(P6020, P6040, P6080, P6087, P6088, P6450, P8587, Wage) %>%
  mutate(`P6040_2` = P6040^2)

Data <- dummy_cols(.data = Data,
                   select_columns = c("P6020", "P6080", "P6087", "P6088", "P6450", "P8587"),
                   remove_first_dummy = TRUE,
                   remove_selected_columns = TRUE)

# 6. Export database

write_csv(x = Data, file = "~/Downloads/Database - Case study 2.txt")
