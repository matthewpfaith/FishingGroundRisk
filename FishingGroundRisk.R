#Matthew Faith, 05/11/2025

#Load libraries
library(readr) 
library(dplyr)
library(tidyr)
library(raster)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(stringr)
library(ggpattern) 
library(stars)     
library(viridis)   
library(scales)    

#Setwd
setwd()

biomass_csv      <- "Prcntg_chng_FishBiomass_refperiod_1990-2000.csv"
effort_csv       <- "Pelagic_NomFishing_201417_mean.csv"
vuln_csv         <- "Fishingground_Totalvulnerability_FINAL.csv"
missing_csv      <- "Missing_Data.csv"

decade_cols      <- c("2040-2050","2050-2060","2060-2070","2070-2080","2080-2090","2090-2100")
nbss_name        <- "NBSS" 

#Functions
norm_minmax <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  mn <- min(x, na.rm = TRUE); mx <- max(x, na.rm = TRUE)
  if (mn == mx) return(rep(0, length(x)))
  (x - mn) / (mx - mn)
}

find_value_col <- function(df) {
  numcols <- names(df)[sapply(df, is.numeric)]
  valcol <- setdiff(numcols, c("Lon","Lat"))[1]
  if (is.na(valcol)) stop("Could not find numeric value column in a driver CSV")
  return(valcol)
}

#Generic Power Transform Function
calc_power_transform <- function(r, power_val) {
  calc(r, fun = function(x){
    if(all(is.na(x))) return(x)
    x_norm <- norm_minmax(x) 
    min_val <- 1e-18
    if(power_val == 1) return(x_norm)
    ifelse(is.na(x_norm), NA_real_, ifelse(x_norm==0,0,(x_norm+min_val)^power_val - min_val))
  })
}

#Load data
bio_all <- read_csv(biomass_csv, show_col_types = FALSE)
effort  <- read_csv(effort_csv, show_col_types = FALSE)
vuln    <- read_csv(vuln_csv, show_col_types = FALSE)
missing <- read_csv(missing_csv, show_col_types = FALSE)

effort_valcol   <- find_value_col(effort)
vuln_valcol     <- "Vulnerability"
missing_valcol  <- find_value_col(missing)

bio_long <- bio_all %>%
  pivot_longer(cols = all_of(decade_cols), names_to = "Decade", values_to = "Change") %>%
  mutate(
    Year_mid = as.numeric(str_extract(Decade, "\\d{4}(?=-)")) + 5,
    Year_scaled = Year_mid - 1995
  )

#Calculate trend - SLOW!
bio_slopes <- bio_long %>%
  group_by(Model, Lon, Lat) %>%
  summarise(
    Slope = {
      df <- cur_data()
      df <- df[!is.na(df$Change), ]
      if (nrow(df) < 2) NA_real_ else {
        m <- try(lm(Change ~ 0 + Year_scaled, data = df)$coefficients[1], silent = TRUE)
        if (inherits(m, "try-error")) NA_real_ else m
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    biomass_decline = ifelse(!is.na(Slope) & Slope < 0, -Slope, 0),
    biomass_norm = norm_minmax(biomass_decline)
  )

#Rasterise
bio_template <- bio_slopes %>% filter(Model == nbss_name)
template_raster <- rasterFromXYZ(bio_template %>% dplyr::select(Lon, Lat, biomass_norm))
crs(template_raster) <- "+proj=longlat +datum=WGS84 +no_defs"

#Exposure
r_effort_raw <- rasterFromXYZ(effort %>% dplyr::select(Lon, Lat, !!effort_valcol) %>% rename(value = !!effort_valcol))
crs(r_effort_raw) <- crs(template_raster)
r_effort_raw <- resample(r_effort_raw, template_raster, method = "bilinear")
r_effort_raw[is.na(r_effort_raw[])] <- 0
#Create Normalised Effort (0-1) used for Calc
r_effort_norm <- calc(r_effort_raw, fun = norm_minmax)

#Vulnerability
r_vuln <- rasterFromXYZ(vuln %>% dplyr::select(Lon, Lat, !!vuln_valcol) %>% rename(value = !!vuln_valcol))
crs(r_vuln) <- crs(template_raster)
r_vuln <- resample(r_vuln, template_raster, method = "bilinear")
r_vuln[is.na(r_vuln[])] <- 0

#Missing mask
r_missing <- rasterFromXYZ(missing %>% dplyr::select(Lon, Lat, !!missing_valcol) %>% rename(value = !!missing_valcol))
crs(r_missing) <- crs(template_raster)
r_missing <- resample(r_missing, template_raster, method = "ngb")

#Compute risk
models <- unique(bio_slopes$Model)
risk_list_raw <- list()

for(mod in models) {
  bio_mod <- bio_slopes %>% filter(Model == mod)
  bio_r <- rasterFromXYZ(bio_mod %>% dplyr::select(Lon, Lat, biomass_norm))
  crs(bio_r) <- crs(template_raster)
  risk_r <- bio_r * r_effort_norm * r_vuln
  risk_r[r_missing == 1] <- NA
  risk_r <- resample(risk_r, template_raster, method = "bilinear")
  risk_list_raw[[mod]] <- risk_r
}

#NBSS-FishMIP consensus
all_names <- names(risk_list_raw)
fishmip_models <- setdiff(all_names, nbss_name)

r_nbss <- risk_list_raw[[nbss_name]]
r_fishmip_stack <- stack(risk_list_raw[fishmip_models])
r_fishmip_mean <- calc(r_fishmip_stack, mean, na.rm = TRUE)

thresh_nbss <- quantile(values(r_nbss), probs = 0.90, na.rm = TRUE)
thresh_fishmip <- quantile(values(r_fishmip_mean), probs = 0.90, na.rm = TRUE)

bin_nbss <- r_nbss >= thresh_nbss
bin_nbss[is.na(bin_nbss)] <- 0
bin_fishmip <- r_fishmip_mean >= thresh_fishmip
bin_fishmip[is.na(bin_fishmip)] <- 0

consensus_raster <- bin_nbss * bin_fishmip 

#Map prep
robinson <- CRS('+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

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

#Transform
prep_for_plot <- function(r, power_val=1, apply_mask=FALSE) {
  if(apply_mask) { r[r_missing == 1] <- NA }
  r <- calc_power_transform(r, power_val)
  r_proj <- projectRaster(r, crs = robinson, method = "bilinear")
  df <- rasterToPoints(r_proj) %>% as.data.frame()
  names(df) <- c("x", "y", "val")
  df$val <- pmax(pmin(df$val, 1), 0)
  return(df)
}

v3_theme <- theme(
  plot.title = element_text(family = "Prata", size = 25, hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5, size = 18),
  plot.caption = element_text(color = "grey50", size = 12, hjust = 0.9),
  legend.position = "bottom",
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank()
)

#GENERATE COMPONENT MAPS
#GENERATE HAZARD MAP
bio_nbss_r <- rasterFromXYZ(bio_slopes %>% filter(Model == nbss_name) %>% dplyr::select(Lon, Lat, biomass_norm))
crs(bio_nbss_r) <- crs(template_raster)
df_hazard <- prep_for_plot(bio_nbss_r, power_val = 1.0, apply_mask = FALSE) 

p_hazard <- ggplot() +
  theme_void() +
  geom_tile(data = df_hazard, aes(x = x, y = y, fill = val)) +
  scale_fill_viridis_c(
    option = 'rocket', 
    name = "Hazard", 
    na.value = "white",
    limits = c(0, 1)
  ) +
  geom_sf(data=countries_robinson, colour="black", linetype='solid', fill = "lightgrey", size=0.5) +
  geom_sf(data=bb_robinson, colour="black", linetype='solid', fill = NA, linewidth= 1) +
  labs(title = "a) Projected Decline (Hazard)") +
  v3_theme

#GENERATE EXPOSURE MAP
df_exposure <- prep_for_plot(r_effort_norm, power_val = 0.15, apply_mask = FALSE)

p_exposure <- ggplot() +
  theme_void() +
  geom_tile(data = df_exposure, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(
    colors = c("#000033", "aquamarine4", "yellowgreen", "yellow", "gold", "darkorange"),
    name = "Exposure",
    limits = c(0, 1)
  ) +
  geom_sf(data=countries_robinson, colour="black", linetype='solid', fill = "lightgrey", size=0.5) +
  geom_sf(data=bb_robinson, colour="black", linetype='solid', fill = NA, linewidth= 1) +
  labs(title = "b) Fishing Effort (Exposure)") +
  v3_theme

#GENERATE VULNERABILITY MAP
df_vuln <- prep_for_plot(r_vuln, power_val = 0.5, apply_mask = TRUE)

p_vuln <- ggplot() +
  theme_void() +
  geom_tile(data = df_vuln, aes(x = x, y = y, fill = val)) +
  scale_fill_viridis_c(
    option = 'mako', 
    name = "Vulnerability",
    na.value = 'white', 
    limits = c(0, 1)
  ) +
  geom_sf(data=countries_robinson, colour="black", linetype='solid', fill = "lightgrey", size=0.5) +
  geom_sf(data=bb_robinson, colour="black", linetype='solid', fill = NA, linewidth= 1) +
  labs(title = "c) Socioeconomic Dependence (Vulnerability)") +
  v3_theme

#GENERATE RISK MAP
df_risk <- prep_for_plot(r_nbss, power_val = 0.15, apply_mask = TRUE)

st_consensus <- st_as_stars(consensus_raster)
sf_consensus <- st_as_sf(st_consensus, as_points = FALSE, merge = TRUE)
names(sf_consensus)[1] <- "value"
sf_consensus <- sf_consensus %>% filter(value == 1)
sf_consensus_proj <- st_transform(sf_consensus, crs = robinson)

p_risk <- ggplot() +
  theme_void() +
  geom_tile(data = df_risk, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(
    colors = c("seagreen", "yellow","orange","orangered", "red", "darkred"), 
    na.value = 'white',
    name = "Risk",
    limits = c(0, 1)
  ) +
  geom_sf_pattern(
    data = sf_consensus_proj,
    fill = NA,                
    colour = NA,              
    pattern = 'stripe',       
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.02,   
    pattern_spacing = 0.02,   
    pattern_alpha = 1,      
    pattern_size = 0.4) +
  geom_sf(data=countries_robinson, colour="black", linetype='solid', fill = "lightgrey", size=0.5) +
  geom_sf(data=bb_robinson, colour="black", linetype='solid', fill = NA, linewidth= 1) +
  labs(title = "d) Fisheries Climate Risk") +
  v3_theme

#plot
plot(p_hazard)
plot(p_exposure)
plot(p_vuln)
plot(p_risk)

#Save Outputs
#ggsave("Figure3a_Hazard_Final.png", p_hazard, width = 8, height = 5, dpi = 600)
#ggsave("Figure3b_Exposure_Final.png", p_exposure, width = 8, height = 5, dpi = 600)
#ggsave("Figure3c_Vulnerability_Final.png", p_vuln, width = 8, height = 5, dpi = 600)

#ggsave("Figure3d_Risk_Final.png", p_risk, width = 8, height = 5, dpi = 600)
