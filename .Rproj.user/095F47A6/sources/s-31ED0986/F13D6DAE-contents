# Analysis for 
# _Estimation of biodiversity with short-range multispectral imaging in a temperate calcareous grassland_
#### Jackson, J, Lawson, C. S., Adelmant, C., Huhtala, E., Fernandes, P., Hodgson, R., King, H., Williamson, L., Maseyk, K., Hawes, N., Hector, A., & Salguero-Gómez, R

#### 2022-09-02
#### Repository created by John Jackson

---

[![DOI](https://zenodo.org/badge/531925353.svg)](https://zenodo.org/badge/latestdoi/531925353)

Full analysis for our study on the link between spectral diversity and biodiversity in a calcareous grassland in Wytham Woods, Oxfordshire. Here we use multispectral image data from a commercially available UAV monitoring system at short range and fine spatial resolution. For the manuscript please see [the bioRxiv entry](https://www.biorxiv.org/content/10.1101/2022.03.08.483493v3). Package version info for this analysis is given below.

Analysis scripts can be found in the `code/` sub-repository, manuscript figures and output in the `output/` sub-repository, and analysis data in the `data/` sub-repository. Raw cropped image data (which is unsuitable for storage online) is available on request from John Jackson.

Scripts are labeled A-E in order of the analysis, and are as follows:

1. `A_biodiversity_data_processing.R` - data cleaning and wrangling of raw biodiversity data for June 2021.
2. `B_raster_extraction_uncalibrated.R` - pipeline for image analysis to calculate spectral distribution characteristics from cropped image data.
3. `C_raster_band_plot.R` - plotting of example raster data for manuscript figure 1. Script C2 is repeated but with a black background.
4. `D_biodiversity_spectral_uncalibrated_models.R` - Bayesian linear modeling framework linking spectral diversity and biodiversity.
5. `E_biodiversity_spectral_uncalibrated_supplementary.R` - supplementary figures and stats for the manuscript.

---

## System Information and Package Versions

<details>
  <summary>Click here to expand</summary>

```
R version 4.0.5 (2021-03-31)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] viridis_0.5.1       viridisLite_0.3.0   rnaturalearth_0.1.0 brms_2.15.0        
 [5] Rcpp_1.0.7          cowplot_1.1.1       ggdist_3.0.0        ggridges_0.5.3     
 [9] flextable_0.6.7     psych_2.1.6         rasterVis_0.50.3    latticeExtra_0.6-29
[13] lattice_0.20-41     terra_1.3-4         raster_3.4-13       sp_1.4-5           
[17] patchwork_1.1.1     forcats_0.5.1       stringr_1.4.0       dplyr_1.0.5        
[21] purrr_0.3.4         readr_1.4.0         tidyr_1.1.3         tibble_3.1.0       
[25] ggplot2_3.3.5       tidyverse_1.3.1    

loaded via a namespace (and not attached):
  [1] uuid_0.1-4           readxl_1.3.1         backports_1.2.1      systemfonts_1.0.2   
  [5] plyr_1.8.6           igraph_1.2.6         splines_4.0.5        crosstalk_1.1.1     
  [9] rstantools_2.1.1     inline_0.3.19        digest_0.6.27        htmltools_0.5.1.1   
 [13] rsconnect_0.8.24     fansi_0.4.2          magrittr_2.0.1       modelr_0.1.8        
 [17] RcppParallel_5.1.4   matrixStats_0.60.0   officer_0.3.19       xts_0.12.1          
 [21] prettyunits_1.1.1    jpeg_0.1-9           colorspace_2.0-0     rvest_1.0.1         
 [25] xfun_0.29            haven_2.3.1          callr_3.7.0          crayon_1.4.1        
 [29] jsonlite_1.7.2       hexbin_1.28.2        lme4_1.1-27.1        zoo_1.8-9           
 [33] glue_1.4.2           gtable_0.3.0         emmeans_1.6.2-1      V8_3.4.2            
 [37] distributional_0.2.2 ggdark_0.2.1         pkgbuild_1.2.0       rstan_2.21.2        
 [41] abind_1.4-5          scales_1.1.1         mvtnorm_1.1-1        DBI_1.1.1           
 [45] miniUI_0.1.1.1       xtable_1.8-4         units_0.7-2          tmvnsim_1.0-2       
 [49] proxy_0.4-26         stats4_4.0.5         StanHeaders_2.21.0-7 DT_0.18             
 [53] htmlwidgets_1.5.3    httr_1.4.2           threejs_0.3.3        RColorBrewer_1.1-2  
 [57] ellipsis_0.3.2       farver_2.1.0         pkgconfig_2.0.3      loo_2.4.1           
 [61] dbplyr_2.1.1         utf8_1.2.1           tidyselect_1.1.1     rlang_0.4.11        
 [65] reshape2_1.4.4       later_1.2.0          munsell_0.5.0        cellranger_1.1.0    
 [69] tools_4.0.5          cli_3.0.1            generics_0.1.0       broom_0.7.9         
 [73] evaluate_0.14        fastmap_1.1.0        knitr_1.33           processx_3.5.2      
 [77] fs_1.5.0             zip_2.2.0            nlme_3.1-152         mime_0.11           
 [81] projpred_2.0.2       xml2_1.3.2           compiler_4.0.5       bayesplot_1.8.1     
 [85] shinythemes_1.2.0    rstudioapi_0.13      curl_4.3.2           gamm4_0.2-6         
 [89] png_0.1-7            e1071_1.7-8          reprex_2.0.1         stringi_1.5.3       
 [93] ps_1.6.0             Brobdingnag_1.2-6    gdtools_0.2.3        Matrix_1.3-2        
 [97] classInt_0.4-3       nloptr_1.2.2.2       markdown_1.1         shinyjs_2.0.0       
[101] vctrs_0.3.8          pillar_1.6.2         lifecycle_1.0.0      bridgesampling_1.1-2
[105] estimability_1.3     data.table_1.14.0    httpuv_1.6.1         R6_2.5.0            
[109] promises_1.2.0.1     KernSmooth_2.23-18   gridExtra_2.3        codetools_0.2-18    
[113] boot_1.3-27          colourpicker_1.1.0   MASS_7.3-53.1        gtools_3.8.2        
[117] assertthat_0.2.1     withr_2.4.2          shinystan_2.5.0      mnormt_2.0.2        
[121] mgcv_1.8-34          parallel_4.0.5       hms_1.1.0            grid_4.0.5          
[125] class_7.3-18         coda_0.19-4          minqa_1.2.4          rmarkdown_2.10      
[129] sf_1.0-2             shiny_1.6.0          lubridate_1.7.10     base64enc_0.1-3     
[133] dygraphs_1.1.1.6           
```
</details>
