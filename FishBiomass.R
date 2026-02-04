#Matthew Faith, 05/11/2025
#INPUT: CMIP6 projections of CHLa, depth-integrated phytoplankton biomass and FishMIP model outputs

#Load libraries
library(readr)
library(dplyr)
library(raster)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(sp)

#Set WD
setwd("C:/Users/matthewfaith/OneDrive - University of Plymouth/PhD/Chapter 5 - Fishing effort for assessing fisheries climate risks/Data/")

#Load ESM datasets
CHLa_ssp585 <- read_csv("Input/chl_mgm3_ssp585_decadal.csv")
Phyto_Biomass_ssp585 <- read_csv("Input/phyc-vin2_gCm2_ssp585_decadal.csv")
CHLa_ssp126 <- read_csv("Input/chl_mgm3_ssp126_decadal.csv")
Phyto_Biomass_ssp126 <- read_csv("Input/phyc-vin2_gCm2_ssp126_decadal.csv")

#Select SSP
CHLa <- CHLa_ssp585
Phyto_Biomass <- Phyto_Biomass_ssp585
selected_ssp <- "ssp585"  # used to filter FishMIP file (case-insensitive)

#Ref period
future_decades <- c("2040-2050","2050-2060","2060-2070","2070-2080","2080-2090","2090-2100")
ref_period <- "1990-2000"

#Functions based on Atkinson et al., 2024
calc_percent_change <- function(df, ref_col, decades) {
  out <- data.frame(Lon = df$Lon, Lat = df$Lat)
  for (dec in decades) {
    out[[dec]] <- (df[[dec]] - df[[ref_col]]) / df[[ref_col]]
  }
  return(out)
}

calc_nbss_slope <- function(chla_vals) {
  return((-0.0065) * (log10(chla_vals)^2) + 0.0663 * log10(chla_vals) + (-1.0651))
}

calc_nbss_CIs <- function(chla_vals, coeff) {
  return(coeff[1] * (log10(chla_vals)^2) + coeff[2] * log10(chla_vals) + coeff[3])
}

calc_fish_biomass <- function(nbss_vals, phyto_vals) {
  return( 49999.5 * 10^((12 * nbss_vals) + log10((phyto_vals / 0.0000000499995))) )
}

#NBSS calcs of supportable fish biomass changes
CHLa_prcntg_change <- calc_percent_change(CHLa, ref_period, future_decades)
Phyto_Biomass_prcntg_change <- calc_percent_change(Phyto_Biomass, ref_period, future_decades)

#NBSS Slopes
NBSS_slope <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
#coeff_main <- c(-0.0065, 0.0663, -1.0651) 

NBSS_slope[[ref_period]] <- calc_nbss_slope(CHLa[[ref_period]])
for (dec in future_decades) NBSS_slope[[dec]] <- calc_nbss_slope(CHLa[[dec]])

#NBSS slope change
NBSS_slope_change <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
for (dec in future_decades) NBSS_slope_change[[dec]] <- NBSS_slope[[dec]] - NBSS_slope[[ref_period]]

#Supportable Fish Biomass (absolute)
Supportable_fish_biomass <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
Supportable_fish_biomass[[ref_period]] <- calc_fish_biomass(NBSS_slope[[ref_period]], Phyto_Biomass[[ref_period]])
for (dec in future_decades) Supportable_fish_biomass[[dec]] <- calc_fish_biomass(NBSS_slope[[dec]], Phyto_Biomass[[dec]])

#Supportable Fish Biomass (percentage change, fractional)
Supportable_fish_biomass_change <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
for (dec in future_decades) {
  Supportable_fish_biomass_change[[dec]] <- (Supportable_fish_biomass[[dec]] - Supportable_fish_biomass[[ref_period]]) / Supportable_fish_biomass[[ref_period]]
}

#NBSS confidence intervals - coefficients from Atksinson et al., 2024
#LCI
coeff_lci <- c(-0.03658449, 0.04260289, -1.082543)
LCI_NBSS_slope <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
LCI_NBSS_slope[[ref_period]] <- calc_nbss_CIs(CHLa[[ref_period]], coeff_lci)
for (dec in future_decades) LCI_NBSS_slope[[dec]] <- calc_nbss_CIs(CHLa[[dec]], coeff_lci)

LCI_Supportable_fish_biomass <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
LCI_Supportable_fish_biomass[[ref_period]] <- calc_fish_biomass(LCI_NBSS_slope[[ref_period]], Phyto_Biomass[[ref_period]])
for (dec in future_decades) LCI_Supportable_fish_biomass[[dec]] <- calc_fish_biomass(LCI_NBSS_slope[[dec]], Phyto_Biomass[[dec]])

LCI_Supportable_fish_biomass_change <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
for (dec in future_decades) {
  LCI_Supportable_fish_biomass_change[[dec]] <- (LCI_Supportable_fish_biomass[[dec]] - LCI_Supportable_fish_biomass[[ref_period]]) / LCI_Supportable_fish_biomass[[ref_period]]
}

#UCI
coeff_uci <- c(0.0235849, 0.09017938, -1.05387)
UCI_NBSS_slope <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
UCI_NBSS_slope[[ref_period]] <- calc_nbss_CIs(CHLa[[ref_period]], coeff_uci)
for (dec in future_decades) UCI_NBSS_slope[[dec]] <- calc_nbss_CIs(CHLa[[dec]], coeff_uci)

UCI_Supportable_fish_biomass <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
UCI_Supportable_fish_biomass[[ref_period]] <- calc_fish_biomass(UCI_NBSS_slope[[ref_period]], Phyto_Biomass[[ref_period]])
for (dec in future_decades) UCI_Supportable_fish_biomass[[dec]] <- calc_fish_biomass(UCI_NBSS_slope[[dec]], Phyto_Biomass[[dec]])

UCI_Supportable_fish_biomass_change <- data.frame(Lon = CHLa$Lon, Lat = CHLa$Lat)
for (dec in future_decades) {
  UCI_Supportable_fish_biomass_change[[dec]] <- (UCI_Supportable_fish_biomass[[dec]] - UCI_Supportable_fish_biomass[[ref_period]]) / UCI_Supportable_fish_biomass[[ref_period]]
}

#Sense check Rasters
CHLa_Change_rast <- rasterFromXYZ(CHLa_prcntg_change)
Phyto_biomass_change_rast <- rasterFromXYZ(Phyto_Biomass_prcntg_change)
NBSS_slope_rast <- rasterFromXYZ(NBSS_slope)
NBSS_slope_change_rast <- rasterFromXYZ(NBSS_slope_change)
Supportable_fish_change_rast <- rasterFromXYZ(Supportable_fish_biomass_change)

#Save CSVs
write.csv(NBSS_slope,"NBSS_slope.csv", row.names = FALSE)
write.csv(CHLa_prcntg_change,"Prcntg_chng_CHLa_refperiod_1990-2000.csv", row.names = FALSE)
write.csv(Phyto_Biomass_prcntg_change,"Prcntg_chng_PhytoBiomass_refperiod_1990-2000.csv", row.names = FALSE)
write.csv(NBSS_slope_change,"Abslt_chng_NBSS_refperiod_1990-2000.csv", row.names = FALSE)
write.csv(Supportable_fish_biomass,"Supportable_FishBiomass.csv", row.names = FALSE)
write.csv(Supportable_fish_biomass_change,"Prcntg_chng_FishBiomass_refperiod_1990-2000.csv", row.names = FALSE)
write.csv(LCI_Supportable_fish_biomass_change,"Prcntg_chng_FishBiomass_refperiod_1990-2000_LCI.csv", row.names = FALSE)
write.csv(UCI_Supportable_fish_biomass_change,"Prcntg_chng_FishBiomass_refperiod_1990-2000_UCI.csv", row.names = FALSE)

#Add in FishMIP MEMs for % change in supportable fish biomass (TCB)
fishmip <- read_csv("Input/fishmip_tcb_delta_all_decadal.csv", col_types = cols())

#Validate decades exist in FishMIP
available_decades <- intersect(future_decades, names(fishmip))
if(length(available_decades) != length(future_decades)) {
  stop("FishMIP file is missing one or more expected decade columns. Expected: ", paste(future_decades, collapse=", "))
}

#Filter to selected SSP
fishmip_ssp <- fishmip %>%
  filter(tolower(SSP) == tolower(selected_ssp))

#Pivot longer, average across ESMs for each Lon, Lat, MEM, decade
fishmip_long_avg <- fishmip_ssp %>%
  dplyr::select(Lon, Lat, MEM, ESM, all_of(future_decades)) %>%
  pivot_longer(cols = all_of(future_decades), names_to = "Decade", values_to = "PctValue") %>%
  group_by(Lon, Lat, MEM, Decade) %>%
  summarise(mean_pct = mean(PctValue, na.rm = TRUE), .groups = "drop")

#Convert percent
fishmip_long_avg <- fishmip_long_avg %>%
  mutate(frac_change = mean_pct)

#Pivot wider to have decades as columns and MEM as Model
fishmip_wide <- fishmip_long_avg %>%
  dplyr::select(-mean_pct) %>%
  pivot_wider(names_from = Decade, values_from = frac_change)

#Rename MEM -> Model (keep names exactly as in MEM)
fishmip_wide <- fishmip_wide %>%
  rename(Model = MEM)

FishMIP_mean <- fishmip_wide %>%
  group_by(Lon, Lat) %>%
  summarise(across(all_of(future_decades), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Model = "FishMIP-Mean") %>%
  dplyr::select(Lon, Lat, Model, all_of(future_decades))

NBSS_wide <- Supportable_fish_biomass_change %>%
  mutate(Model = "NBSS") %>%
  dplyr::select(Lon, Lat, Model, all_of(future_decades))

fishmip_wide <- fishmip_wide %>%
  dplyr::select(Lon, Lat, Model, all_of(future_decades))

combined_pct_change <- bind_rows(NBSS_wide, fishmip_wide, FishMIP_mean)

#Ensure consistentcy
combined_pct_change <- combined_pct_change %>%
  mutate(across(all_of(future_decades), ~ ifelse(is.nan(.x), NA_real_, .x)))

#Overwrite
write.csv(combined_pct_change, "Prcntg_chng_FishBiomass_refperiod_1990-2000.csv", row.names = FALSE)

#Plotting
#Choose model and decade to plot
model_to_plot <- "NBSS"   #e.g., "NBSS", "FishMIP-Mean", "apecosm", "boats", "dbem", "dbpm"           
#"ecoocean", "ecotroph", "feisty", "macroecological", "zoomss"
decade_to_plot <- "2050-2060"

#Choose colour limits
color_limits <- c(-0.6, 0.6)

#Extract the rows for the chosen model
plot_df <- combined_pct_change %>%
  filter(Model == model_to_plot) %>%
  dplyr::select(Lon, Lat, all_of(decade_to_plot)) %>%
  dplyr::rename(Change = !!sym(decade_to_plot))

if(nrow(plot_df) == 0) stop("No rows found for model_to_plot = '", model_to_plot, "'. Check spelling / MEM names in fishmip_tcb_delta_all_decadal.csv")

#Create raster from the Lon/Lat/Change grid
plot_df_clean <- plot_df %>% filter(!is.na(Change))

#If plot_df_clean is empty, warn and stop
if(nrow(plot_df_clean) == 0) stop("All values are NA for model '", model_to_plot, "' and decade '", decade_to_plot, "'")

layer <- rasterFromXYZ(plot_df_clean)
crs(layer) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#Robinson projection object
robinson <- CRS('+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

#Project the raster to Robinson
layer_proj <- projectRaster(layer, crs = robinson, method = "bilinear")

#Convert projected raster to points
pts <- rasterToPoints(layer_proj)
df_pts <- as.data.frame(pts)

names(df_pts)[3] <- "val"

#Clip extreme values to the fixed color limits
df_pts$val[df_pts$val < color_limits[1]] <- color_limits[1]
df_pts$val[df_pts$val > color_limits[2]] <- color_limits[2]

bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -179.999, xmax = 179.999, ymax = 90, ymin = -90), crs = st_crs(4326)),
  n = 100))
bb_robinson <- st_transform(bb, as.character(robinson))

countries <- ne_countries(scale = 110, returnclass = c("sf"))
sf_use_s2(FALSE)
countries_robinson <- countries %>%
  st_buffer(0) %>%
  st_intersection(st_union(bb)) %>%
  st_transform(robinson)

plot_title <- "zoomss (SSP 5-8.5)"
plot_sub <- decade_to_plot

plot_fixed <- ggplot() +
  theme_void() +
  geom_tile(data = df_pts, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(colors = c("red","white","blue"),
                       limits = color_limits,
                       name = '% Change') +
  geom_sf(data = countries_robinson,
          colour = 'white',
          linetype = 'solid',
          fill = 'black',
          size = 0.2) +
  geom_sf(data = bb_robinson,
          colour = 'black',
          linetype = 'solid',
          fill = NA,
          linewidth = 1) +
  theme(plot.title = element_text(family = "Prata", size = 30, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 25),
        plot.caption = element_text(color = "black", size = 16, hjust = 0.9),
        legend.title = element_text(size = 20),
        legend.position = 'bottom',
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 15)) +
  labs(title = plot_title,
       subtitle = plot_sub,
       caption = "")

print(plot_fixed)

#Save the plot to file
#ggsave(filename = paste0("Map_", model_to_plot, "_", gsub("-", "", decade_to_plot), ".png"),
#      plot = plot_fixed, width = 12, height = 7, dpi = 300)
