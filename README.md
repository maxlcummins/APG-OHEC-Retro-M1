[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14048330.svg)](https://doi.org/10.5281/zenodo.14048330)

# APG-OHEC-1

A repository containing supplementary scripts for data processing and visualisation in Watt and Cummins et al, 2024

## Usage

First, clone the repository and install the required packages using renv:

```         
git clone https://github.com/maxlcummins/APG-OHEC-Retro-M1.git
```

Then, open the R project and run the following:

```         
#Install our packages
renv::restore()
```

This will request a prompt of yes or no to install. Click the output tag below to see what it will look like.

<details>

<summary>Output</summary>

```         
The following package(s) will be updated:

# Bioconductor ---------------------------------------------------------------

-   BiocGenerics [\* -\> 0.48.1]
-   ComplexHeatmap [\* -\> 2.18.0]
-   IRanges [\* -\> 2.36.0]
-   treeio [\* -\> 1.26.0]

# Bioconductor 3.18 ----------------------------------------------------------

-   BiocVersion [\* -\> 3.18.1]
-   ggtree [\* -\> 3.10.1]
-   S4Vectors [\* -\> 0.40.2]

# CRAN -----------------------------------------------------------------------

-   BiocManager [1.30.25 -\> 1.30.22]
-   abind [\* -\> 1.4-5]
-   admisc [\* -\> 0.35]
-   ape [\* -\> 5.8]
-   aplot [\* -\> 0.2.2]
-   askpass [\* -\> 1.2.0]
-   assertthat [\* -\> 0.2.1]
-   backports [\* -\> 1.4.1]
-   base64enc [\* -\> 0.1-3]
-   bit [\* -\> 4.0.5]
-   bit64 [\* -\> 4.0.5]
-   blob [\* -\> 1.2.4]
-   brew [\* -\> 1.0-10]
-   brio [\* -\> 1.1.5]
-   broom [\* -\> 1.0.5]
-   bslib [\* -\> 0.7.0]
-   cachem [\* -\> 1.1.0]
-   callr [\* -\> 3.7.6]
-   car [\* -\> 3.1-2]
-   carData [\* -\> 3.0-5]
-   caret [\* -\> 6.0-94]
-   cellranger [\* -\> 1.1.0]
-   circlize [\* -\> 0.4.16]
-   cli [\* -\> 3.6.2]
-   clipr [\* -\> 0.8.0]
-   clock [\* -\> 0.7.0]
-   clue [\* -\> 0.3-65]
-   colorspace [\* -\> 2.1-0]
-   commonmark [\* -\> 1.9.1]
-   ComplexUpset [\* -\> 1.3.3]
-   conflicted [\* -\> 1.2.0]
-   corrplot [\* -\> 0.92]
-   cowplot [\* -\> 1.1.3]
-   cpp11 [\* -\> 0.4.7]
-   crayon [\* -\> 1.5.2]
-   credentials [\* -\> 2.0.1]
-   curl [\* -\> 5.2.1]
-   data.table [\* -\> 1.15.4]
-   DBI [\* -\> 1.2.2]
-   dbplyr [\* -\> 2.5.0]
-   desc [\* -\> 1.4.3]
-   devtools [\* -\> 2.4.5]
-   diagram [\* -\> 1.6.5]
-   diffobj [\* -\> 0.3.5]
-   digest [\* -\> 0.6.35]
-   doParallel [\* -\> 1.0.17]
-   downlit [\* -\> 0.4.3]
-   dplyr [\* -\> 1.1.4]
-   dtplyr [\* -\> 1.3.1]
-   e1071 [\* -\> 1.7-14]
-   ellipsis [\* -\> 0.3.2]
-   eulerr [\* -\> 7.0.2]
-   evaluate [\* -\> 0.23]
-   fansi [\* -\> 1.0.6]
-   farver [\* -\> 2.1.2]
-   fastmap [\* -\> 1.2.0]
-   fontawesome [\* -\> 0.5.2]
-   forcats [\* -\> 1.0.0]
-   foreach [\* -\> 1.5.2]
-   fs [\* -\> 1.6.4]
-   future [\* -\> 1.33.2]
-   future.apply [\* -\> 1.11.2]
-   gargle [\* -\> 1.5.2]
-   generics [\* -\> 0.1.3]
-   GenSA [\* -\> 1.1.14]
-   gert [\* -\> 2.0.1]
-   GetoptLong [\* -\> 1.0.5]
-   ggalluvial [\* -\> 0.12.5]
-   ggfun [\* -\> 0.1.4]
-   ggplot2 [\* -\> 3.5.1]
-   ggplotify [\* -\> 0.1.2]
-   ggpubr [\* -\> 0.6.0]
-   ggrepel [\* -\> 0.9.5]
-   ggsci [\* -\> 3.0.3]
-   ggsignif [\* -\> 0.6.4]
-   gh [\* -\> 1.4.1]
-   gitcreds [\* -\> 0.1.2]
-   GlobalOptions [\* -\> 0.1.2]
-   globals [\* -\> 0.16.3]
-   glue [\* -\> 1.7.0]
-   googledrive [\* -\> 2.1.1]
-   googlesheets4 [\* -\> 1.1.1]
-   gower [\* -\> 1.0.1]
-   gridExtra [\* -\> 2.3]
-   gridGraphics [\* -\> 0.5-1]
-   gtable [\* -\> 0.3.5]
-   hardhat [\* -\> 1.3.1]
-   harrietr [\* -\> 0.2.3]
-   haven [\* -\> 2.5.4]
-   here [\* -\> 1.0.1]
-   highcharter [\* -\> 0.9.4]
-   highr [\* -\> 0.10]
-   hms [\* -\> 1.1.3]
-   htmltools [\* -\> 0.5.8.1]
-   htmlwidgets [\* -\> 1.6.4]
-   httpuv [\* -\> 1.6.15]
-   httr [\* -\> 1.4.7]
-   httr2 [\* -\> 1.0.1]
-   hues [\* -\> 0.2.0]
-   ids [\* -\> 1.0.1]
-   igraph [\* -\> 2.0.3]
-   ini [\* -\> 0.3.1]
-   ipred [\* -\> 0.9-14]
-   isoband [\* -\> 0.2.7]
-   iterators [\* -\> 1.0.14]
-   jquerylib [\* -\> 0.1.4]
-   jsonlite [\* -\> 1.8.8]
-   knitr [\* -\> 1.46]
-   labeling [\* -\> 0.4.3]
-   later [\* -\> 1.3.2]
-   lava [\* -\> 1.8.0]
-   lazyeval [\* -\> 0.2.2]
-   lifecycle [\* -\> 1.0.4]
-   listenv [\* -\> 0.9.1]
-   lme4 [\* -\> 1.1-35.3]
-   lubridate [\* -\> 1.9.3]
-   magrittr [\* -\> 2.0.3]
-   MatrixModels [\* -\> 0.5-3]
-   matrixStats [\* -\> 1.3.0]
-   memoise [\* -\> 2.0.1]
-   mime [\* -\> 0.12]
-   miniUI [\* -\> 0.1.1.1]
-   minqa [\* -\> 1.2.6]
-   ModelMetrics [\* -\> 1.2.2.2]
-   modelr [\* -\> 0.1.11]
-   munsell [\* -\> 0.5.1]
-   nloptr [\* -\> 2.0.3]
-   numDeriv [\* -\> 2016.8-1.1]
-   openssl [\* -\> 2.1.2]
-   paletteer [\* -\> 1.6.0]
-   parallelly [\* -\> 1.37.1]
-   patchwork [\* -\> 1.2.0]
-   pbkrtest [\* -\> 0.5.2]
-   pheatmap [\* -\> 1.0.12]
-   pillar [\* -\> 1.9.0]
-   pkgbuild [\* -\> 1.4.4]
-   pkgconfig [\* -\> 2.0.3]
-   pkgdown [\* -\> 2.0.9]
-   pkgload [\* -\> 1.3.4]
-   plyr [\* -\> 1.8.9]
-   png [\* -\> 0.1-8]
-   polyclip [\* -\> 1.10-6]
-   polylabelr [\* -\> 0.2.0]
-   polynom [\* -\> 1.4-1]
-   praise [\* -\> 1.0.0]
-   prettyunits [\* -\> 1.2.0]
-   prismatic [\* -\> 1.1.2]
-   pROC [\* -\> 1.18.5]
-   processx [\* -\> 3.8.4]
-   prodlim [\* -\> 2023.08.28]
-   profvis [\* -\> 0.3.8]
-   progress [\* -\> 1.2.3]
-   progressr [\* -\> 0.14.0]
-   promises [\* -\> 1.3.0]
-   proxy [\* -\> 0.4-27]
-   ps [\* -\> 1.7.6]
-   purrr [\* -\> 1.0.2]
-   quantmod [\* -\> 0.4.26]
-   quantreg [\* -\> 5.97]
-   R6 [\* -\> 2.5.1]
-   ragg [\* -\> 1.3.0]
-   rappdirs [\* -\> 0.3.3]
-   rcmdcheck [\* -\> 1.4.0]
-   RColorBrewer [\* -\> 1.1-3]
-   Rcpp [\* -\> 1.0.12]
-   RcppArmadillo [\* -\> 0.12.8.2.1]
-   RcppEigen [\* -\> 0.3.4.0.0]
-   readr [\* -\> 2.1.5]
-   readxl [\* -\> 1.4.3]
-   recipes [\* -\> 1.0.10]
-   rematch [\* -\> 2.0.0]
-   rematch2 [\* -\> 2.1.2]
-   remotes [\* -\> 2.5.0]
-   reprex [\* -\> 2.1.0]
-   reshape2 [\* -\> 1.4.4]
-   rjson [\* -\> 0.2.21]
-   rlang [\* -\> 1.1.3]
-   rlist [\* -\> 0.4.6.2]
-   rmarkdown [\* -\> 2.26]
-   roxygen2 [\* -\> 7.3.1]
-   rprojroot [\* -\> 2.0.4]
-   rstatix [\* -\> 0.7.2]
-   rstudioapi [\* -\> 0.16.0]
-   rversions [\* -\> 2.1.2]
-   rvest [\* -\> 1.0.4]
-   sass [\* -\> 0.4.9]
-   scales [\* -\> 1.3.0]
-   selectr [\* -\> 0.4-2]
-   sessioninfo [\* -\> 1.2.2]
-   shape [\* -\> 1.4.6.1]
-   shiny [\* -\> 1.8.1.1]
-   sourcetools [\* -\> 0.1.7-1]
-   SparseM [\* -\> 1.81]
-   splitstackshape [\* -\> 1.4.8]
-   SQUAREM [\* -\> 2021.1]
-   stringi [\* -\> 1.8.3]
-   stringr [\* -\> 1.5.1]
-   sys [\* -\> 3.4.2]
-   systemfonts [\* -\> 1.0.6]
-   testthat [\* -\> 3.2.1.1]
-   textshaping [\* -\> 0.3.7]
-   tibble [\* -\> 3.2.1]
-   tidyr [\* -\> 1.3.1]
-   tidyselect [\* -\> 1.2.1]
-   tidytree [\* -\> 0.4.6]
-   tidyverse [\* -\> 2.0.0]
-   timechange [\* -\> 0.3.0]
-   timeDate [\* -\> 4032.109]
-   tinytex [\* -\> 0.50]
-   TTR [\* -\> 0.24.4]
-   tzdb [\* -\> 0.4.0]
-   urlchecker [\* -\> 1.0.1]
-   usethis [\* -\> 2.2.3]
-   utf8 [\* -\> 1.2.4]
-   uuid [\* -\> 1.2-0]
-   vctrs [\* -\> 0.6.5]
-   venn [\* -\> 1.12]
-   viridisLite [\* -\> 0.4.2]
-   vroom [\* -\> 1.6.5]
-   waldo [\* -\> 0.5.2]
-   whisker [\* -\> 0.4.1]
-   withr [\* -\> 3.0.0]
-   xfun [\* -\> 0.43]
-   XML [\* -\> 3.99-0.16.1]
-   xml2 [\* -\> 1.3.6]
-   xopen [\* -\> 1.0.1]
-   xtable [\* -\> 1.8-4]
-   xts [\* -\> 0.13.2]
-   yaml [\* -\> 2.3.8]
-   yulab.utils [\* -\> 0.1.4]
-   zip [\* -\> 2.3.1]
-   zoo [\* -\> 1.8-12]

# GitHub ---------------------------------------------------------------------

-   abricateR [\* -\> [maxlcummins/abricateR\@HEAD](mailto:maxlcummins/abricateR@HEAD){.email}]
-   ggVennDiagram [\* -\> [gaospecial/ggVennDiagram\@HEAD](mailto:gaospecial/ggVennDiagram@HEAD){.email}]

Do you want to proceed? [Y/n]:
```

</details>

## Quality Control

Four genomes failed QC as they typed by MLST as being non-*E. coli*. Run the `markdown/AusTrakka_TEA_Quality_Control.Rmd` script to generate the QC figures and see which genomes failed in the resulting files in the delims folder

## Data Aggregation

Data generated by [pipelord](https://github.com/maxlcummins/pipelord) is aggregated using the `markdown/AusTrakka_TEA_Data_Aggregation.Rmd` script. This script generates the data used in the manuscript. Database versions utilised are available in `databases`.

## Data Visualisation

Data is visualised using the `markdown/General_analysis.Rmd` script. This script generates the figures used in the manuscript.
