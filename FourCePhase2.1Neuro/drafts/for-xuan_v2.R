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
#library("devtools")
#install_github("andland/logisticPCA")
library(logisticPCA)
library(survival)

setwd("/4ceData/Phase2.1NeuroRPackage/FourCePhase2.1Neuro")

set.seed(446) # for obfuscation posterity

select <- dplyr::select
filter <- dplyr::filter

# load site params
load("data/site_params.rda")

# load icd10 codes
load("data/neuro_icds_10.rda")

# load in necessary R functions
source("R/run.R")
source("R/summary-stats.R")
source("R/comorbidity-map-funcs.R")
source("R/utils.R")

## Review the below parameters to make sure they are set appropriately for your site ##

# Set to FALSE if not running VA data
VA_site <- FALSE

# set to 4CE Data directory
data_dir = "/4ceData/Input"

# change siteID
currSiteId = "NWU"

#### Start running analysis ####

CurrSiteId <- toupper(currSiteId)
site_specs <- site_params %>%
  dplyr::filter(tolower(siteid) == tolower(currSiteId))

mask_thres <- site_specs$mask_thres
blur_abs <- site_specs$blur_abs
icd_version <- site_specs$icd_version
include_race <- site_specs$include_race

if (VA_site) {
  clin_raw <- clin_raw %>%
    mutate(patient_num = as.character(patient_num))
  obs_raw <- obs_raw %>%
    mutate(patient_num = as.character(patient_num))
  demo_raw <- demo_raw %>%
    mutate(
      patient_num = as.character(patient_num),
      across(
        ends_with("_date"),
        ~ na_if(., "1900-01-01")
      )
    )
} else {
  clin_raw <-
    readr::read_csv(
      file.path(data_dir, "LocalPatientClinicalCourse.csv"),
      col_types = list(patient_num = readr::col_character())
    )
  obs_raw <-
    readr::read_csv(
      file.path(data_dir, "LocalPatientObservations.csv"),
      col_types = list(patient_num = readr::col_character())
    )
  demo_raw <-
    readr::read_csv(
      file.path(data_dir, "LocalPatientSummary.csv"),
      col_types = list(patient_num = readr::col_character()),
      na = c("1900-01-01", "1/1/1900")
    )
}

if (CurrSiteId == "MGB") {
  clin_raw <- clin_raw %>%
    group_by(patient_num) %>%
    arrange(patient_num, days_since_admission) %>%
    mutate(in_hospital = case_when(
      in_hospital == 1 & lead(in_hospital) == 0 & lag(in_hospital) == 0 ~ 0,
      TRUE ~ in_hospital
    )) %>%
    ungroup()
}

demo_raw <- demo_raw %>%
  mutate(
    across(
      ends_with("_date") & tidywhere(is.character),
      ~ lubridate::parse_date_time(.x, orders = c("mdy", "ymd"))
    ),
    #death_date = if_else(death_date > Sys.Date(), NA, as.Date(death_date)),
    last_discharge_date = pmin(death_date, last_discharge_date, na.rm = TRUE),
    time_to_last_discharge = subtract_days(admission_date, last_discharge_date)
  )

obs_raw <- obs_raw %>%
  filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9")) %>%
  mutate(concept_code = stringr::str_sub(concept_code, 1, 3))

if (icd_version == 9) {
  neuro_icds <- neuro_icds_9
} else {
  neuro_icds <- neuro_icds_10
}

#####
# start analysis
comp_readmissions <- clin_raw %>%
  group_by(patient_num) %>%
  arrange(days_since_admission, .by_group = TRUE) %>%
  mutate(delta_hospitalized = diff(c(in_hospital[1], in_hospital))) %>%
  mutate(
    first_out =
      delta_hospitalized == -1 & !duplicated(delta_hospitalized == -1),
    first_change =
      first_out |
      (delta_hospitalized == 1 &
         !duplicated(delta_hospitalized == 1)),
    first_discharge_date = case_when(first_out==TRUE ~ lag(calendar_date)),
    first_discharge_date = min(first_discharge_date, na.rm = TRUE)
  ) %>%
  ungroup()

n_readms <- comp_readmissions %>%
  filter(
    delta_hospitalized != 0,
    in_hospital == 1
  ) %>%
  add_count(patient_num, name = "n_readmissions") %>%
  arrange(desc(n_readmissions)) %>%
  select(patient_num, n_readmissions) %>%
  distinct()

readmissions <- comp_readmissions %>%
  filter(patient_num %in% n_readms$patient_num, first_change) %>%
  select(patient_num, delta_hospitalized, days_since_admission) %>%
  pivot_wider(
    names_from = delta_hospitalized,
    values_from = days_since_admission
  ) %>%
  mutate(time_to_first_readmission = `1` - `-1`) %>%
  select(patient_num, time_to_first_readmission) %>%
  left_join(n_readms, by = "patient_num")

nstay_df <- comp_readmissions %>%
  filter(first_out) %>%
  transmute(patient_num, time_to_first_discharge = days_since_admission - 1)

# patients with neuro conditions prior to admission
pre_neuro <- obs_raw %>%
  filter(days_since_admission >= -365 & days_since_admission <= -15) %>%
  right_join(neuro_icds, by = c("concept_code" = "icd")) %>%
  filter(!is.na(patient_num)) %>%
  count(patient_num, pns_cns) %>%
  mutate(value = log(n + 1)) %>%
  pivot_wider(id_cols = patient_num, names_from = pns_cns, values_from = value, values_fill = 0) %>%
  rename(pre_admission_cns = Central,
         pre_admission_pns = Peripheral)

demo_processed_first <- demo_raw %>%
  mutate(
    time_to_severe = subtract_days(admission_date, severe_date),
    time_to_severe = if_else(time_to_severe < 0, NA_real_, time_to_severe),
    time_to_death = subtract_days(admission_date, death_date),
    time_to_death = if_else(time_to_death < 0, NA_real_, time_to_death),
    readmitted = patient_num %in% readmissions$patient_num,
    sex = as.factor(sex),
    race = as.factor(race),
    age_group = as.factor(age_group),
    Severity = as.factor(severe) %>%
      fct_recode(Severe = "1", `Non-severe` = "0"),
    Survival = as.factor(deceased) %>%
      fct_recode(Alive = "0", Deceased = "1")
  ) %>%
  left_join(nstay_df, by = "patient_num") %>%
  left_join(comp_readmissions %>% distinct(patient_num, first_discharge_date), by = "patient_num") %>%
  left_join(readmissions, by = "patient_num") %>%
  left_join(pre_neuro, by = "patient_num") %>%
  replace_na(list(n_readmissions = 0,
                  pre_admission_cns = 0,
                  pre_admission_pns = 0))

# handle cases when a site may not have any cases of pre_admission codes
if ("pre_admission_cns" %in% colnames(demo_processed_first)) {
  print("pre_admission_cns available")
} else {
  demo_processed_first$pre_admission_cns <- 0
}

if ("pre_admission_pns" %in% colnames(demo_processed_first)) {
  print("pre_admission_pns available")
} else {
  demo_processed_first$pre_admission_pns <- 0
}


# identify patient's who are still in the hospital
# this information will help us select the icd codes during first hospitalization
in_hospital <- demo_processed_first %>%
  filter(is.na(time_to_first_discharge)) %>%
  distinct(patient_num, still_in_hospital)

temp_neuro <-
  temporal_neuro(comp_readmissions, obs_raw, neuro_icds, readmissions, in_hospital)

obs_first_hosp <- temp_neuro$obs_first_hosp

# 8.31.2021 - will remove this function for now. I don't think this is excluding our "Both" patients
# propagated_codes <- temp_neuro$propagated_codes %>%
#   blur_it(c("n_early_codes", "n_new_codes"), blur_abs, mask_thres) %>%
#   mutate(
#     prop_new_codes = if_else(n_early_codes == 0, 0, prop_new_codes),
#     prob_repeated = if_else(n_early_codes == 0, 0, prob_repeated),
#     prob_at_least_one_new =
#       if_else(n_early_codes == 0, 0, prob_at_least_one_new)
#   )

# need to identify patients with "Both" neuro cns and pns codes during neuro admission
neuro_patients <- obs_first_hosp %>%
  filter(days_since_admission >= 0) %>%
  right_join(neuro_icds, by = c("concept_code" = "icd")) %>%
  filter(!is.na(patient_num)) %>%
  # for each patient, calculate time to first cns and time to first pns
  group_by(patient_num, pns_cns) %>%
  mutate(time_to_cns = ifelse(pns_cns == "Central", min(days_since_admission), NA),
         time_to_pns = ifelse(pns_cns == "Peripheral", min(days_since_admission), NA)) %>%
  distinct(patient_num, concept_code, pns_cns, time_to_cns, time_to_pns) %>%
  ungroup() %>%
  group_by(patient_num) %>%
  mutate(nerv_sys_count = length(unique(pns_cns))) %>%
  ungroup() %>%
  mutate(neuro_type = case_when(
    nerv_sys_count == 2 ~ "Both",
    TRUE ~ as.character(pns_cns)
  ))

# count number of Both patients
both_pts <- neuro_patients %>%
  filter(neuro_type == "Both") %>%
  distinct(patient_num)

both <- both_pts %>%
  count() %>%
  as.integer() %>%
  blur_mask_int_num(., blur_abs, mask_thres)

# remove both from neuro_patients df
neuro_patients <- neuro_patients %>%
  filter(!neuro_type == "Both") %>%
  # consolidate the 'time_to_cns' and 'time_to_pns' as 'time_to_neuro'
  mutate(time_to_neuro = if_else(pns_cns == "Central", time_to_cns, time_to_pns)) %>%
  select(-time_to_cns, -time_to_pns)

neuro_pt_post <- unique(neuro_patients$patient_num)

# calculate mean and median number of codes per patient and neuro_type status
n_codes_per_patient <- neuro_patients %>%
  group_by(patient_num, neuro_type) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  distinct(patient_num, neuro_type, count) %>%
  select(neuro_type, count) %>%
  group_by(neuro_type) %>%
  mutate(mean_count = mean(count),
         median_count = median(count),
         iqr25 = quantile(count, probs = 0.25),
         iqr75 = quantile(count, probs = 0.75)) %>%
  ungroup() %>%
  distinct(neuro_type, mean_count, iqr25, median_count, iqr75)

non_neuro_patients <-
  data.frame(patient_num = setdiff(demo_processed_first$patient_num, neuro_pt_post)) %>%
  mutate(concept_code = "NN") %>%
  # remove both patients
  filter(!patient_num %in% both_pts$patient_num)

# remove both patients from obs_first_hosp and demo_processed_first dataframe
obs_first_hosp <- obs_first_hosp %>% filter(!patient_num %in% both_pts$patient_num)
demo_processed_first <- demo_processed_first %>% filter(!patient_num %in% both_pts$patient_num)
demo_raw <- demo_raw %>% filter(!patient_num %in% both_pts$patient_num)

comorb_list <- get_elix_mat(obs_first_hosp, icd_version)

index_scores_elix <- comorb_list$index_scores_elix

# ensure data is formatted correctly
index_scores_elix$patient_num <- as.character(index_scores_elix$patient_num)
demo_raw$patient_num <- as.character(demo_raw$patient_num)
obs_raw$patient_num <- as.character(obs_raw$patient_num)
demo_processed_first$patient_num <- as.character(demo_processed_first$patient_num)
index_scores_elix$patient_num <- as.character(index_scores_elix$patient_num)
nstay_df$patient_num <- as.character(nstay_df$patient_num)

index_scores_elix <- index_scores_elix %>%
  right_join0(select(demo_raw, patient_num), by = "patient_num")

mapped_codes_table <- comorb_list$mapped_codes_table

elix_mat <- cor(select(
  index_scores_elix,
  -c(patient_num, elixhauser_score)
))

# individual comorbidity matrixes
cns_pts <- neuro_patients %>%
  filter(pns_cns == "Central") %>%
  distinct(patient_num)

pns_pts <- neuro_patients %>%
  filter(pns_cns == "Peripheral") %>%
  distinct(patient_num)

index_scores_elix_cns <- index_scores_elix %>% filter(patient_num %in% cns_pts$patient_num)
index_scores_elix_pns <- index_scores_elix %>% filter(patient_num %in% pns_pts$patient_num)

elix_mat_cns <-
  cor(select(
    index_scores_elix_cns,
    -c(patient_num, elixhauser_score)
  ))

elix_mat_pns <-
  cor(select(
    index_scores_elix_pns,
    -c(patient_num, elixhauser_score)
  ))

elix_pca <- index_scores_elix %>%
  select(-elixhauser_score) %>%
  tibble::column_to_rownames("patient_num") %>%
  as.matrix()

lpca_fit <- logisticPCA::logisticPCA(elix_pca, k = min(nrow(elix_pca), 10), m = 0)
# k = 10 principal components, m is solved for

deviance_expl <- lpca_fit$prop_deviance_expl

pca_covariates <- lpca_fit$PCs %>%
  data.frame() %>%
  `colnames<-`(paste0(".fittedPC", 1:10)) %>%
  tibble::rownames_to_column("patient_num")

obs_raw = obs_first_hosp
demo_processed = demo_processed_first

### Part 2: PNS vs CNS
message("Start CNS vs PNS Analysis")
binary = FALSE
print(binary == FALSE)

nstay_df <- neuro_patients %>%
  data.frame() %>%
  select(patient_num, concept_code) %>%
  rbind(non_neuro_patients) %>%
  left_join(demo_processed, by = "patient_num") %>%
  mutate(concept_code = fct_reorder(concept_code, time_to_first_discharge)) %>%
  left_join(neuro_icds, by = c("concept_code" = "icd"))

comorb_names_elix <- get_quan_elix_names()

icd_tables <- get_tables(
  c("no_neuro_cond", "neuro_cond"),
  nstay_df,
  right_join0(index_scores_elix, nstay_df, by = "patient_num"),
  comorb_names_elix,
  blur_abs,
  mask_thres,
  "concept_code"
)[-3] # last element is not useful, remove

neuro_types <- c("None", "Peripheral", "Central")

demo_df <- demo_processed %>%
  filter(!patient_num %in% both_pts$patient_num) %>%
  left_join(distinct(select(neuro_patients, patient_num, neuro_type, time_to_neuro)),
            by = "patient_num"
  ) %>%
  replace_na(list(neuro_type = "None")) %>%
  mutate(neuro_post = forcats::fct_relevel(neuro_type, neuro_types))

scores_unique <- index_scores_elix %>%
  right_join0(demo_df, by = "patient_num") %>%
  left_join(pca_covariates, by = "patient_num")

obfus_tables <- get_tables(
  neuro_types,
  demo_df,
  scores_unique,
  comorb_names_elix,
  blur_abs,
  mask_thres
) %>%
  lapply(function(x) mutate(x, site = currSiteId))

## -------------------------------------------------------------------------
# Create individual tables for those who met outcome on admission and are
# exclude from main survival analysis

surv_exclude_pts <- function(df, time_to_outcome) {

  if (time_to_outcome == "time_to_severe") {
    demo_subset_df <- df %>%
      filter(time_to_severe == 0 | time_to_death == 0)
  } else if (time_to_outcome == "time_to_death") {
    demo_subset_df <- df %>%
      filter(time_to_death == 0)
  } else if (time_to_outcome == "time_to_first_discharge") {
    demo_subset_df <- df %>%
      mutate(death_before_outcome = case_when(time_to_death <= 0 ~ 1,
                                              TRUE ~ 0)) %>%
      filter(time_to_first_discharge == 0 & death_before_outcome == 0)
  } else if (time_to_outcome == "time_to_last_discharge") {
    demo_subset_df <- df %>%
      mutate(death_before_outcome = case_when(time_to_death <= 0 ~ 1,
                                              TRUE ~ 0)) %>%
      filter(time_to_last_discharge == 0 & death_before_outcome == 0)
  }

  # dataframe for the survival analysis
  scores_unique <- index_scores_elix %>%
    right_join0(demo_subset_df, by = "patient_num") %>%
    left_join(pca_covariates, by = "patient_num")

  obfus_tables <- tryCatch(
    {
      get_tables(
        neuro_types,
        demo_subset_df,
        scores_unique,
        comorb_names_elix,
        blur_abs,
        mask_thres
      ) %>%
        lapply(function(x) mutate(x, site = currSiteId))
    },
    error = function(cond) {
      message("Original error message:")
      message(cond)
      message("No data to subset. Skipping for now...")
      return(NULL) # return NA in case of error
    }
  )
  return(obfus_tables)
}

severe_adm <- surv_exclude_pts(demo_df, "time_to_severe")
death_adm <- surv_exclude_pts(demo_df, "time_to_death")
first_adm <- surv_exclude_pts(demo_df, "time_to_first_discharge")
last_adm <- surv_exclude_pts(demo_df, "time_to_last_discharge")

# Run Xuan's new code in from-xuan-runcox0725022.r

df <- scores_unique
ind_vars <- get_ind_vars(df, include_race)
tcut=60
source('from-xuan-runcox07252022.r')

xuan_results <- list()
deceased_results <- run_coxregression(df, 'deceased', ind_vars, tcut=tcut)
severe_results <- run_coxregression(df, 'severe', ind_vars, tcut=tcut)
los_results <- run_coxregression(df, 'time_to_first_discharge', ind_vars, tcut=tcut)

xuan_results <- c("deceased_results" = deceased_results,
                  "severe_results" = severe_results,
                  "los_results" = los_results)


## Old results ------------------------------------------------------------

## -------------------------------------------------------------------------
#reg_results <- run_regressions(scores_unique, include_race)

## Run Cox-Regression models
sub_reg_results <- run_coxregressions(scores_unique, include_race, blur_abs, mask_thres)

## ----save-results---------------------------------------------------------

run_coxregression <- function(df, depend_var, ind_vars, blur_abs, mask_thres) {
  if (length(unique(df[, depend_var, drop = T])) <= 1) {
    return(NULL)
  }

  independ_vars <- paste(ind_vars, collapse = " + ")

  df <- df %>%
    mutate(across(starts_with("time_to"), as.numeric),
           days_since_admission = as.numeric(days_since_admission))

  if (depend_var == "deceased") {
    surv_df <- df %>%
      mutate(time = if_else(deceased == 1, time_to_death, time_to_last_discharge),
             time = if_else((is.na(time_to_last_discharge) & deceased == 0), days_since_admission, time),
             delta = deceased) %>%
      select(patient_num, delta, time, all_of(ind_vars)) %>%
      filter(!(time == 0 & delta == 1))
  } else if (depend_var == "severe") {
    surv_df <- df %>%
      mutate(
        time = if_else(severe == 1 | deceased == 1, pmin(time_to_severe, time_to_death, na.rm = TRUE), time_to_last_discharge),
        time = if_else(is.na(time), days_since_admission, time),
        delta = if_else(severe == 1 | deceased == 1, 1, 0)) %>%
      select(patient_num, delta, time, all_of(ind_vars)) %>%
      filter(!(time == 0 & delta == 1))
  } else if (depend_var == "time_to_last_discharge") {
    surv_df <- df %>%
      mutate(
        time = case_when(
          time_to_death <= time_to_last_discharge ~ time_to_death,
          is.na(time_to_last_discharge) & is.na(time_to_death) ~ days_since_admission,
          TRUE ~ time_to_last_discharge),
        delta = case_when(
          is.na(time_to_last_discharge) & is.na(time_to_death) ~ 3, #patients still in hospital
          time_to_death <= time_to_last_discharge ~ 2, #censored patients (died)
          TRUE ~ 1 #patients who are discharged
        )
      ) %>%
      select(patient_num, delta, time, all_of(ind_vars)) %>%
      filter(!(time == 0 & delta == 1))
  } else if (depend_var == "time_to_first_discharge") {
    surv_df <- df %>%
      mutate(
        time = case_when(
          time_to_death <= time_to_first_discharge ~ time_to_death,
          is.na(time_to_first_discharge) & is.na(time_to_death) ~ days_since_admission,
          TRUE ~ time_to_first_discharge),
        delta = case_when(
          is.na(time_to_first_discharge) & is.na(time_to_death) ~ 3, #patients still in hospital
          time_to_death <= time_to_first_discharge ~ 2, #censored patients (died)
          TRUE ~ 1 #patients who are discharged
        )
      ) %>%
      select(patient_num, delta, time, all_of(ind_vars)) %>%
      filter(!(time == 0 & delta == 1))

  }

  output <- tryCatch(
    {
      list(
        cox = "survival::Surv(time, delta==1)" %>%
          paste("~", independ_vars) %>%
          as.formula() %>%
          survival::coxph(data = surv_df) %>% #blocking out inclusion of #id = patient_num) to see if it fixes BCH problem
          summary(),
        life = "survival::Surv(time, delta==1) ~ neuro_post" %>%
          as.formula() %>%
          survival::survfit(data = surv_df) %>%
          summary()
      )
    },
    error = function(cond) {
      message(paste("Error when regressing", depend_var))
      message("Original error message:")
      message(cond)
      message("Skipping for now...")
      return(NULL) # return NA in case of error
    }
  )

  message('Print length(output). This should be two')
  print(length(output))
  message('Print names(output). This should be cox and life')
  print(names(output))

  event_table_obfs <- tryCatch(

    {

      if (is.null(dim(output$life$n.censor)) == TRUE)  {
        print("n.censor table is correct")
        print(names(output$life))
      } else {
        output$life$n.censor = data.frame(x=unlist(as.data.frame(output$life$n.censor)))
        colnames(output$life$n.censor) <- "n.censor"
        print(names(output$life))
      }


      message("generating event_tables for cox model")
      event_table = data.frame(status = output$life$strata, time = output$life$time, n.risk = output$life$n.risk, n.event = output$life$n.event, n.censor = output$life$n.censor)

      # mask and blur for obfuscation
      message("blurring event_tables for adjusted survival")
      event_table_obfs <- blur_it(event_table, vars = c("n.risk", "n.event", "n.censor"), blur_abs, mask_thres)
    },
    error = function(cond) {
      message(paste("Error when regressing", depend_var))
      message("Original error message:")
      message(cond)
      message("Skipping for now...")
      return(NULL) # return NA in case of error
    }
  )

  output$event_table_obfs <- event_table_obfs

  if (!is.null(output)) {
    output$cox$deviance.resid <- NULL
    output$cox$na.action <- NULL
    output$cox$terms <- NULL
    output$cox$residuals <- NULL
    output$cox$n <- NULL
    output$cox$nevent <- NULL
    output$cox$y <- NULL
    output$cox$linear.predictors <- NULL
    output$life$n.censor <- NULL
    output$life$n <- NULL
    output$life$n.event <- NULL
    output$life$n.risk <- NULL
    output$life$table <- NULL
  }

  if (length(unique(surv_df$neuro_post)) == 3) {

    # average survival functions
    average_survival <- tryCatch(

      {

        # [,-1] removes the intercept term
        covariate=model.matrix(as.formula(paste("survival::Surv(time,delta==1)", '~',
                                                independ_vars)),data=surv_df )[,-1]

        data=data.frame( cbind('time'=surv_df$time,'delta'=surv_df$delta,covariate) )

        # (-1:2) removes time and delta columns
        cox=survival::coxph(as.formula(paste("survival::Surv(time,delta==1)", '~',
                                             paste(colnames(data[,-(1:2)]),collapse='+'))),data=data)

        # calculate mean per covariate
        meancovariate=apply(covariate,2,mean)

        # create a new data frame to create indicators for neuro status
        # each neuro status will have the same mean for each covariate
        newdata=data.frame(rbind(c(0,0,meancovariate[-(1:2)]),
                                 c(1,0,meancovariate[-(1:2)]),
                                 c(0,1,meancovariate[-(1:2)])))
        colnames(newdata)=names(meancovariate)
        survout=survival::survfit(cox,newdata)
        survout = survout %>% summary()

        message("generating event_tables for adjusted survival curves")
        event_table1 = data.frame(time = survout$time, n.risk = survout$n.risk, n.event = survout$n.event, n.censor = survout$n.censor, neuro_status = "neuro_postNone")
        event_table2 = data.frame(time = survout$time, n.risk = survout$n.risk, n.event = survout$n.event, n.censor = survout$n.censor, neuro_status = "neuro_postPeripheral")
        event_table3 = data.frame(time = survout$time, n.risk = survout$n.risk, n.event = survout$n.event, n.censor = survout$n.censor, neuro_status = "neuro_postCentral")

        event_table_surv_adjust <- rbind(event_table1, event_table2, event_table3)

        # mask and blur for obfuscation
        message("blurring event_tables for adjusted survival curves")
        event_table_obfs <- blur_it(event_table_surv_adjust, vars = c("n.risk", "n.event", "n.censor"), blur_abs, mask_thres)

        cox <- cox %>% summary()

        average_survival = list("cox"=cox,"survf"=survout, "event_table_obfs" = event_table_obfs)

      },
      error = function(cond) {
        message(paste("Error when regressing", depend_var))
        message("Original error message:")
        message(cond)
        message('Skipping for now...')
        return(NULL) # return NA in case of error
      }
    )
    if (!is.null(average_survival)){
      average_survival$cox$deviance.resid <- NULL
      average_survival$cox$na.action <- NULL
      average_survival$cox$terms <- NULL
      average_survival$cox$residuals <- NULL
      average_survival$cox$n <- NULL
      average_survival$cox$nevent <- NULL
      average_survival$cox$y <- NULL
      average_survival$cox$linear.predictors <- NULL
      average_survival$survf$n <- NULL
      average_survival$survf$n.censor <- NULL
      average_survival$life$n.censor <- NULL
      average_survival$survf$n.event <- NULL
      average_survival$survf$n.risk <- NULL
    }

    output$average_survival <- average_survival
  }

  return(output)

}


run_coxregressions <- function(df, include_race = TRUE, blur_abs, mask_thres) {
  ind_vars <- get_ind_vars(df, include_race)

  time_severe_reg_elix <-
    run_coxregression(df, "severe", ind_vars, blur_abs, mask_thres)

  time_deceased_reg_elix <-
    run_coxregression(df, "deceased", ind_vars, blur_abs, mask_thres)

  time_last_discharge_reg_elix <-
    run_coxregression(df, "time_to_last_discharge", ind_vars, blur_abs, mask_thres)

  time_first_discharge_reg_elix <-
    run_coxregression(df, "time_to_first_discharge", ind_vars, blur_abs, mask_thres)

  list(
    time_severe_reg_elix = time_severe_reg_elix,
    time_deceased_reg_elix = time_deceased_reg_elix,
    time_last_discharge_reg_elix = time_last_discharge_reg_elix,
    time_first_discharge_reg_elix = time_first_discharge_reg_elix
  )
}


