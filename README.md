# Pacific–Antarctic Ridge Meander Analysis  

This repository provides the datasets, MATLAB functions, and figures used in the analysis of the **characteristics and long-term trends of the Pacific–Antarctic Ridge (PAR) meander** in the Southern Ocean over the **1993–2023** period.

The analyses combine satellite altimetry, RG Argo climatology, and CMEMS surface geostrophic velocities to document meander structure, variability, and associated eddy kinetic energy (EKE).  

---

## Repository structure  

### Core analysis repositories  
- **meander_loc_mon_1993_2023_mean_199303_202303**  
  Contains monthly mean meander positions for 1993–2023 derived from altimetry with probability-based frontal detection.

- **meander_monthly_locations_north_south_boundaries**  
  Provides the north/south latitude boundaries of the meander, enabling calculation of meander width and variability across longitude.  

- **pacific_antarctic_ridge_meander_analysis (this repository)**  
  Includes diagnostic MATLAB scripts and functions for processing, trend estimation, and figure generation.

---

## Key MATLAB functions  

### Gyre and Argo
- **`01-PAR-Argo-Temperature-2004-2023-Build-4D-Cube.m`**  
  Constructs 4-D absolute temperature fields (time × lon × lat × pressure) from RG Argo mean + anomalies (2004–2023), concatenates monthly NetCDF inputs, and saves MATLAB datasets.

- **`02-CMEMS-Gyre-ADT-TimeSeries-Trend-1993-2023.m`**  
  Builds an area-mean monthly ADT anomaly time series for the subtropical gyre (42°–38°S, 150°E–70°W), deseasons it, and computes linear and non-parametric trends (Sen, MK, modified MK, OLS).  

### EKE and geostrophic currents  
- **`04-PAR-CMEMS-EKE-Trend-Analysis-1993-2023.m`**  
  Concatenates CMEMS ugos/vgos files (1993–2023), produces monthly running-mean velocity fields, derives eddy kinetic energy (EKE), maps its 1993–2023 mean, estimates decadal EKE trends, and extracts section-based time series with Sen/MK diagnostics.  

### Meander dynamics  
- **`05-Meander-Width-Speed-Trend-Decomposition.m`**  
  Tests whether meander widening is explained by along-jet speed strengthening. Performs width–speed regression, decomposes width trends into explained and residual parts, and outputs publication-quality figures of raw, explained, and residual trends.  

### Statistical tools  
- **`03-Mann-Kendall-Trend-Test.m`**  
  Classical Mann–Kendall test for monotonic trends. Returns hypothesis test decision, p-value, and S statistic.  

- **`06-Modified-Mann-Kendall-Trend-Test.m`**  
  Modified Mann–Kendall test (Hamed & Rao, 1998) with autocorrelation adjustment. Provides robust trend significance for autocorrelated series.  

- **`08-Theil-Sen-Slope.m`**  
  Theil–Sen slope estimator (median of pairwise slopes). Robust trend slope estimate used throughout the repository.  

### Plotting helpers  
- **`07-Plot-TimeSeries-With-Trend.m`**  
  Generates time series plots with Sen slope lines, Mann–Kendall significance annotations, and publication-ready formatting.  

---

## Figures  

Figures produced with these functions illustrate:  
- **Meander position and trajectory** from altimetry,  
- **Meander width and intensity variability** (north/south boundaries),  
- **Trends in along-jet geostrophic current speed**,  
- **Eddy kinetic energy (EKE) patterns and changes**,  
- **Temperature structure** from RG Argo climatology.

- ---

## License  

MIT License — free to use, modify, and share with attribution.  
