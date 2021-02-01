Analysis of COVID-19 patients with Neurological Conditions
================

## Install

You can install the development version of **Phase2.1NeuroRPackage**
from GitHub with remotes:

``` r
# install.packages('remotes') # uncomment to install devtools
remotes::install_github('covidclinical/Phase2.1NeuroRPackage',
                        subdir = 'FourCePhase2.1Neuro',
                        upgrade = FALSE)
```

## How to run

The main function `runAnalysis()` has 4 required arguments. See your
site obfuscation parameters (`mask_thres` and `blur_abs`)
[here](https://docs.google.com/spreadsheets/d/1Xl9juDBXt86P3xQtsoTaBl2zPl1BIiAG9DI3Rotyqp8/edit#gid=212461777).
If the site’s ICD codes are primarily of version 9, set
`icd_version = 9`. If your site does not have race information, set
`include_race = FALSE`.

``` r
library(Phase2.1NeuroRPackage)

runAnalysis(mask_thres = 10, blur_abs = 0, icd_version = 10, include_race = TRUE)
```

Finally, please submit the results to the central repository:

``` r
submitAnalysis()
```

If you run into any problem adapting this code to your data, let us
(@meghutch and @trang1618) know via Slack or [submit an
issue](https://github.com/covidclinical/Phase2.1NeuroRPackage/issues/new).

Thank you very much for your support!
