# Trait-Filtered Responses of Mammal Communities to Land-Use Change in a Neotropical Dry Forest

## Instructions for Data Usage and Reproducibility

This repository contains the data and code to reproduce the analysis of the associated publication.  
The files are structured as follows:

├── Data/
│   ├── PanTHERIA.txt                  # Species trait data from PanTHERIA
│   ├── bounding_boxes_before_consensus.csv  # Bounding boxes from raw images
│   ├── camop_problem.csv              # Camera operation table
│   ├── camtraps_clean.csv             # Camera trap metadata
│   ├── forest_covariates.csv          # Forest metrics calculated from remote sensing
│   └── species_data.csv               # Post-consensus species data used for modelling
│
├── Remote Sensing/
│   ├── 10S_070W.tif                   # Reference raster for forest metric classification
│   ├── NDVI (Annual means)/           # Stacked annual NDVI means
│   └── site_covs_rs_prep.R            # R script to process remote sensing data
│
└── trait_filtering_analysis.R          # Main script for the analysis

### Notes

- **NDVI products**: We only provide the annual NDVI means, not the raw Sentinel products. The raw products can be retrieved directly from Sentinel using the product names listed in Supporting Material S-3 or requested from the corresponding author.  
- **Camera trap locations**: Exact lat/lon coordinates are not made public to protect ongoing monitoring efforts.  
- **Raw images**: Access to the raw camera trap images can be provided upon reasonable request to the corresponding author.