# FishingGroundRisk

Scripts to generate the risk figures found in: **Faith, M.P., Atkinson, A., Heneghan, R.F., Ostle, C., Fernandes-Salvador, J.A., Thompson, M.S., Serra-Pompei, C., Artioli, Y., Schmidt, K., Rees, S., Holland, M. and McQuatters-Gollop, A. 2025. Mapping global fisheries climate risks to guide sustainable marine management. bioRxiv, https://doi.org/10.1101/2025.05.07.652593**

To generate fish biomass projections, use CMIP6 outputs (the ones used in our study can be found at https://www.isimip.org/) of depth-integrated phytoplankton biomass and Chlorophyll-a, with a reference period (1990-1999 was used in our study) and/or FishMIP model outputs. FishMIP model outputs described in Tittensor et al. (2021; https://doi.org/10.1038/s41558-021-01173-9) can also be found here: https://www.isimip.org/ To generate risk maps, a list of input socioeconomic metrics (and source datasets) can be found in the methods section of: https://doi.org/10.1101/2025.05.07.652593 and fishing effort data can be found at https://doi.org/10.25959/MNGY-0Q43 from the paper Rousseau, Y., Blanchard, J.L., Novaglio, C. et al. A database of mapped global fishing activity 1950–2017. Sci Data 11, 48 (2024). https://doi.org/10.1038/s41597-023-02824-6

# To reproduce the figures from our paper, follow this order: 

**1** - Generate projections of fish biomass (hazard) from ESM output using the scripts found in **FishBiomass**.

![Figure 2](Images/Figure2(fish).png)

**2** - Calculate country-specific and global average fishing effort data (exposure) from https://doi.org/10.25959/MNGY-0Q43 using the scripts found in **Fishingdata_prep**. This script also uses meta data found in Rousseau et al. (2024).

**3** - Generate fishing ground vulnerability data using socioeconomic indicators from https://doi.org/10.1101/2025.05.07.652593 and (country-specific) global average fishing efforts (outputs from **2**) using the scripts found in **VulnerabilityMapping**.

**4** - Generate fishing ground risk maps using projected fish biomass from **1**, fishing effort from **2**, and vulerability outputs from **3** using the scripts found in **FishingGroundRisk**.

![Figure 3](Images/Figure3(ClimateRisk).png)
