#Matthew Faith, 05/11/2025

#Load libraries
library(readr)
library(dplyr)
library(tidyr)

#Load data
setwd()
master_data <- read_csv("Input/GriddedEffortby_FGroup_FishingCountry_Sector.csv")

#Data prep
data_clean <- master_data %>% 
  dplyr::select(-1) %>% 
  mutate(FGroup = as.factor(FGroup),
    Year = as.factor(Year))

#Filter to years 2014-2017
data_clean <- data_clean %>% 
  filter(Year %in% c("2014", "2015", "2016", "2017"))

#All fishing effort
All_effort_summary <- data_clean %>%
  mutate(loc = paste(Lat, Lon)) %>%
  group_by(Lat, Lon, loc) %>%
  summarise(Total_4Yr_Effort = sum(NomActiveHours, na.rm = TRUE), .groups = 'drop') %>%
  mutate(All_Nom_201417_mean_effort = Total_4Yr_Effort / 4) %>% #In fishing database used NaN = 0 fishing
  dplyr::select(Lat, Lon, loc, All_Nom_201417_mean_effort)

#Save
write_csv(All_effort_summary, "All_NomFishing_201417_mean.csv")

#Total Pelagic effort
pelagic_groups <- c("pelagic<30cm", "pelagic30-90cm", "pelagic>=90cm")

Pelagic_effort_summary <- data_clean %>%
  filter(FGroup %in% pelagic_groups) %>%
  mutate(loc = paste(Lat, Lon)) %>%
  group_by(Lat, Lon, loc) %>%
  summarise(Total_4Yr_Effort = sum(NomActiveHours, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Pelagic_Nom_201417_mean_effort = Total_4Yr_Effort / 4) %>%
  dplyr::select(Lat, Lon, loc, Pelagic_Nom_201417_mean_effort)

#Save
write_csv(Pelagic_effort_summary, "Pelagic_NomFishing_201417_mean.csv")

#Effort by SAUP (for vulnerability mapping)
Pelagic_effort_summary_SAUP <- data_clean %>%
  filter(FGroup %in% pelagic_groups) %>%
  mutate(SAUP = as.factor(SAUP)) %>% 
  mutate(loc = paste(Lat, Lon)) %>%
  group_by(Lat, Lon, loc, SAUP) %>%
  summarise(Total_4Yr_Effort = sum(NomActiveHours, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Pelagic_Nom_201417_mean_effort = Total_4Yr_Effort / 4) %>%
  dplyr::select(Lat, Lon, loc, SAUP, Pelagic_Nom_201417_mean_effort)

#Save
write_csv(Pelagic_effort_summary_SAUP, "Pelagic_Nom_201417_mean_bySAUP.csv")

#Calculate Percentage of global fishing effort is pelagic
Total_Global_Effort <- sum(All_effort_summary$All_Nom_201417_mean_effort)
Total_Pelagic_Effort <- sum(Pelagic_effort_summary$Pelagic_Nom_201417_mean_effort)
Percent_Pelagic <- (Total_Pelagic_Effort / Total_Global_Effort) * 100
print(paste0("Global Pelagic Percentage: ", round(Percent_Pelagic, 2), "%"))

