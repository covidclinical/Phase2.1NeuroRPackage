# manual set up for running/testing the package

library(tidyr)
library(readr)
library(forcats)
library(tibble)
library(dplyr)
library(purrr)
library(remotes)
library(stringr)
library(lubridate)
library(broom)
library("devtools")
#install_github("andland/logisticPCA")
library(logisticPCA)
library(survival)
library(survminer)

setwd("/4ceData/Phase2.1NeuroRPackage/FourCePhase2.1Neuro")
# load site params
load("data/site_params.rda")

# load icd10 codes
load("data/neuro_icds_10.rda")

# load in necessary R functions
source("R/run.R")
source("R/summary-stats.R")
source("R/comorbidity-map-funcs.R")
source("R/utils.R")

#####################################

## When testing obfuscation

# site_params <- site_params %>%
#   mutate(mask_thres = if_else(siteid == "NWU", 3, mask_thres),
#          blur_abs = if_else(siteid == "NWU", 2, blur_abs))

##

# Run code in the runAnalysis up to the last part where we save everything into the results() object

#####################################

# Then run:
######################################

obs_raw = obs_first_hosp
demo_processed = demo_processed_first

#######################################

### Part 2: PNS vs CNS
message("Start CNS vs PNS Analysis")

neuro_types <- c("None", "Peripheral", "Central")


###
df = scores_unique_adults
include_race = TRUE
tcut=30
is_pediatric = FALSE
elix = "LPCA"
survtable_interval = 10
analysis = "Adults_30_days_lpca"
