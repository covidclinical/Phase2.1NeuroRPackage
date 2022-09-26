#' Run analysis
#'
#' @return NULL. Result files are written out to the `getProjectOutputDirectory()` directory.
#'
#' @keywords 4CE
#' @export
#' @import dplyr
#' @importFrom stats cor glm lm median sd setNames as.formula
#' @importFrom forcats fct_recode fct_reorder
#' @importFrom tidyr pivot_longer pivot_wider replace_na
#'
runAnalysis <- function(is_docker = TRUE, currSiteId=NULL, data_dir= "/4ceData/Input", output_dir=NULL) {

  #sink("analysis_output.txt")

  # set timer for analysis
  print('Set timer - this analysis will take some time to run')
  time.start.analysis=Sys.time()

  # set seed to maintain obfuscation prosterity
  set.seed(446)

  # specify use of dplyr for `select()` and `filter()` functions
  select <- dplyr::select
  filter <- dplyr::filter

  ## make sure this instance has the latest version of the quality control and data wrangling code available
  # remotes::install_github("covidclinical/Phase2.1DataRPackage@77d32fe", subdir="FourCePhase2.1Data", upgrade=FALSE)

  ## get the site identifier assocaited with the files stored in the /4ceData/Input directory that
  ## is mounted to the container

  if(is_docker == TRUE) {
    data_dir <- FourCePhase2.1Data::getInputDataDirectoryName()
    currSiteId <- FourCePhase2.1Data::getSiteId()
  } else {
    data_dir = data_dir
    currSiteId = currSiteId
  }

  CurrSiteId <- toupper(currSiteId)
  site_specs <- site_params %>%
    dplyr::filter(tolower(siteid) == tolower(currSiteId))


  # specify obfuscation, blurring, and site specific params
  mask_thres <- site_specs$mask_thres
  blur_abs <- site_specs$blur_abs
  icd_version <- site_specs$icd_version
  include_race <- site_specs$include_race

  ## run the quality control
  print('run QC')
  remotes::install_github(
    "https://github.com/covidclinical/Phase2.1DataRPackage",
    subdir = "FourCePhase2.1Data", upgrade = FALSE
  )

  FourCePhase2.1Data::runQC(currSiteId)

  ## load and apply pre-processing of 4CE Phase2.1 files
  # apply special parameters for VA sites
  VA_site <- FALSE

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

  # remove patients who were admitted after October 31, 2021 in order to avoid including patients with omicron variant.
  demo_raw <- demo_raw %>%
    filter(!admission_date >= '2021-11-01')


  # apply special params for MGB sites
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

  # load in and process all 4CE data from other sites
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

  # add time_period to demo_raw
  #demo_raw <- calc_time_period(demo_raw)

  obs_raw <- obs_raw %>%
    filter(concept_type %in% c(paste0("DIAG-ICD",icd_version))) %>%
    mutate(concept_code = stringr::str_sub(concept_code, 1, 3))

  if (icd_version == 9) {
    neuro_icds <- neuro_icds_9
  } else {
    neuro_icds <- neuro_icds_10
  }

  ## estimate last admission date
  site_last_admission_date <- demo_raw %>%
    arrange(admission_date) %>%
    tail(., 1) %>%
    select(admission_date)

  site_last_admission_date <- format(site_last_admission_date$admission_date,"%m-%Y")

  # eval whether a site only has clin_raw `in_hospital`==1
  if(all(clin_raw$in_hospital == 1)) {
    clin_raw <- clin_raw %>%
      group_by(patient_num) %>%
      group_modify(count_sequences_hospitalisation) %>%
      select(siteid, patient_num, days_since_admission, calendar_date, in_hospital,
             severe, deceased)
  }


  #####
  ## start analysis

  # compute admission timelines
  print('compute readmission timelines')
  tryCatch({
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
      first_discharge_date = case_when(first_out==TRUE ~ lag(calendar_date))) %>%
    ungroup()

  # identify first_discharge_date for each patient
  print('determine first discharge date')
  df_first_discharge_date = comp_readmissions %>%
         group_by(patient_num) %>%
         arrange(days_since_admission, .by_group = TRUE) %>%
         filter(!is.na(first_discharge_date)) %>%
         slice(1L)

  # add back to main comp_readmissions df
  print('add first discharge date back to readmission df')
  comp_readmissions <- comp_readmissions %>%
         select(-first_discharge_date) %>%
         left_join(., df_first_discharge_date %>%
                     select(patient_num, first_discharge_date),
                                   by = "patient_num")
  },
  error = function(cond) {
    message("Original error message:")
    message(cond)
    return(NULL) # return NA in case of error
  }
  )

  # compute number of readmissions
  n_readms <- comp_readmissions %>%
    filter(
      delta_hospitalized != 0,
      in_hospital == 1
    ) %>%
    add_count(patient_num, name = "n_readmissions") %>%
    arrange(desc(n_readmissions)) %>%
    select(patient_num, n_readmissions) %>%
    distinct()

  # compute time_to_first_readmission
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

  # determine time_to_first_discharge
  nstay_df <- comp_readmissions %>%
    filter(first_out) %>%
    transmute(patient_num, time_to_first_discharge = days_since_admission - 1)

  # identify patients with neuro conditions prior to admission
  pre_neuro <- obs_raw %>%
    filter(days_since_admission >= -365 & days_since_admission <= -15) %>%
    right_join(neuro_icds, by = c("concept_code" = "icd")) %>%
    filter(!is.na(patient_num)) %>%
    count(patient_num, pns_cns) %>%
    mutate(value = log(n + 1)) %>%
    pivot_wider(id_cols = patient_num,
                names_from = pns_cns,
                values_from = value,
                values_fill = 0) %>%
    rename(pre_admission_cns = Central,
           pre_admission_pns = Peripheral)

  # prepare demographics df and determine time to events for our outcomes
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
                    pre_admission_pns = 0)) %>%
    # determine patients discharged from covid admission
    mutate(covid_discharged = if_else(still_in_hospital == 1 & is.na(as.factor(first_discharge_date)), "Not Discharged", "Discharged"))

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

  temp_neuro <- temporal_neuro(comp_readmissions, obs_raw, neuro_icds, readmissions, in_hospital)

  obs_first_hosp <- temp_neuro$obs_first_hosp

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

  # count both patients by adults & pediatrics
  both_adult_pts <- demo_processed_first %>%
    filter(patient_num %in% both_pts$patient_num,
           age_group %in% c("18to25", "26to49", "50to69", "70to79", "80plus")) %>%
    distinct(patient_num) %>%
    count() %>%
    as.integer() %>%
    blur_mask_int_num(., blur_abs, mask_thres)

  both_ped_pts <- demo_processed_first %>%
    filter(patient_num %in% both_pts$patient_num,
           age_group %in% c("00to02", "06to11", "12to17")) %>%
    distinct(patient_num) %>%
    count() %>%
    as.integer() %>%
    blur_mask_int_num(., blur_abs, mask_thres)

  # save all 'both' counts - these have all been obfuscated via `blur_mask_int_num`
  both_counts <- list(both = both,
                      both_adult_pts = both_adult_pts,
                      both_ped_pts = both_ped_pts)

  # remove all of the 'both' patients from neuro_patients df
  neuro_patients <- neuro_patients %>%
    filter(!neuro_type == "Both") %>%
    # consolidate the 'time_to_cns' and 'time_to_pns' as 'time_to_neuro'
    mutate(time_to_neuro = if_else(pns_cns == "Central", time_to_cns, time_to_pns)) %>%
    select(-time_to_cns, -time_to_pns)

  neuro_pt_post <- unique(neuro_patients$patient_num)

  # indicate adult vs pediatric patients
  demo_processed_first <- demo_processed_first %>%
    mutate(adult_ped = ifelse(age_group %in% c("18to25", "26to49", "50to69", "70to79", "80plus"), "adult", "pediatric"))


  # calculate mean and median number of codes per patient and neuro_type status
  n_codes_per_patient <- neuro_patients %>%
    left_join(., demo_processed_first %>%
                select(patient_num, adult_ped),
              by = "patient_num") %>%
    group_by(patient_num, neuro_type) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    distinct(patient_num, neuro_type, adult_ped, count) %>%
    select(neuro_type, adult_ped, count) %>%
    group_by(neuro_type, adult_ped) %>%
    mutate(mean_count = mean(count),
           median_count = median(count),
           iqr25 = quantile(count, probs = 0.25),
           iqr75 = quantile(count, probs = 0.75)) %>%
    ungroup() %>%
    distinct(neuro_type, adult_ped, mean_count, iqr25, median_count, iqr75)

  non_neuro_patients <-
    data.frame(patient_num = setdiff(demo_processed_first$patient_num, neuro_pt_post)) %>%
    mutate(concept_code = "NN") %>%
    # remove both patients
    filter(!patient_num %in% both_pts$patient_num)

  # remove both patients from obs_first_hosp and demo_processed_first dataframe
  obs_first_hosp <- obs_first_hosp %>% filter(!patient_num %in% both_pts$patient_num)
  demo_processed_first <- demo_processed_first %>% filter(!patient_num %in% both_pts$patient_num)
  demo_raw <- demo_raw %>% filter(!patient_num %in% both_pts$patient_num)
  pre_neuro <- pre_neuro %>% filter(!patient_num %in% both_pts$patient_num)

  # ensure data is formatted correctly
  demo_raw$patient_num <- as.character(demo_raw$patient_num)
  obs_raw$patient_num <- as.character(obs_raw$patient_num)
  demo_processed_first$patient_num <- as.character(demo_processed_first$patient_num)
  nstay_df$patient_num <- as.character(nstay_df$patient_num)

  ## additional formatting of neuro patient data for modeling
  nstay_df <- neuro_patients %>%
    data.frame() %>%
    select(patient_num, concept_code) %>%
    rbind(non_neuro_patients) %>%
    left_join(demo_processed_first, by = "patient_num") %>%
    mutate(concept_code = fct_reorder(concept_code, time_to_first_discharge)) %>%
    left_join(neuro_icds, by = c("concept_code" = "icd"))

  # process comorbidity data by adults & pediatrics
  # first, create separate dataframes for adult & pediatric patients
  adult_obs <- demo_processed_first %>%
    filter(adult_ped == "adult") %>%
    distinct(patient_num) %>%
    left_join(., obs_first_hosp)

  ped_obs <- demo_processed_first %>%
    filter(adult_ped == "pediatric") %>%
    distinct(patient_num) %>%
    left_join(., obs_first_hosp)


  ## debugging
  # check that numbers align with what we expect
  # nrow(unique(data.frame(demo_raw$patient_num)))
  # nrow(unique(data.frame(obs_raw$patient_num)))
  # nrow(unique(data.frame(demo_processed_first$patient_num)))
  # nrow(unique(data.frame(obs_first_hosp$patient_num)))
  # nrow(unique(data.frame(adult_obs$patient_num)))
  # nrow(unique(data.frame(ped_obs$patient_num)))
  # table(demo_processed_first$adult_ped)


  # process comorbidity data for survival analysis covariates and comorobidity/neuro risk analysis
  print("processing adult comorbidities")
  tryCatch({
    comorb_adults <- process_comorb_data(adult_obs, demo_raw, nstay_df, neuro_patients, icd_version, is_pediatric = FALSE, blur_abs, mask_thres)
    },
  error = function(cond) {
    message("Original error message:")
    message(cond)
    return(NULL) # return NA in case of error
  }
    )
  print("processing pediatric comorbidities")
  tryCatch({
    comorb_pediatrics <- process_comorb_data(ped_obs, demo_raw, nstay_df, neuro_patients, icd_version, is_pediatric = TRUE, blur_abs, mask_thres)
  },
  error = function(cond) {
    message("Original error message:")
    message(cond)
    return(NULL) # return NA in case of error
  }
  )

  # manually add empty dataframes if they do not exist due to lack of pediatric/adult patients
  # this will prevent future downstream errors with consolidating the list of results
    if(!exists('comorb_adults')) {
      comorb_adults <- data.frame()
    }

    if(!exists('comorb_pediatrics')) {
      comorb_pediatrics <- data.frame()
    }

  # perform the run_hosps to compute our primary analysis and save all results to a list object
  results <- list(
    site = CurrSiteId,
    site_last_admission_date = site_last_admission_date,
    comorbidities = list(comorb_adults = comorb_adults,
                         comorb_pediatrics = comorb_pediatrics),
    first_hosp_results = run_hosps(
      both_pts,
      both_counts,
      mask_thres,
      blur_abs,
      include_race,
      currSiteId,
      readmissions,
      demo_processed_first,
      obs_first_hosp,
      neuro_patients,
      neuro_icds,
      comorb_adults,
      comorb_pediatrics
    )
  )

  # # remove additional icd_tables from the comorb_* lists - these are saved under first_hosp_results
  results$comorbidities$comorb_adults$icd_tables <- NULL
  results$comorbidities$comorb_pediatrics$icd_tables <- NULL

  # remove list variables with PHI
  results$comorbidities$comorb_adults$pca_covariates <- NULL
  results$comorbidities$comorb_adults$index_scores_elix <- NULL

  results$comorbidities$comorb_pediatrics$pca_covariates <- NULL
  results$comorbidities$comorb_pediatrics$index_scores_elix <- NULL

  # process the pre-CNS/PNS code statistics
  # stratify by age group
  pre_neuro_adults <- demo_processed_first %>%
    filter(adult_ped == "adult") %>%
    distinct(patient_num, pre_admission_cns, pre_admission_pns)

  pre_neuro_pediatrics <- demo_processed_first %>%
    filter(adult_ped == "pediatric") %>%
    distinct(patient_num, pre_admission_cns, pre_admission_pns)

  pre_cns_adults_summary <- summary(pre_neuro_adults$pre_admission_cns)
  pre_pns_adults_summary <- summary(pre_neuro_adults$pre_admission_pns)

  pre_cns_pediatrics_summary <- summary(pre_neuro_pediatrics$pre_admission_cns)
  pre_pns_pediatrics_summary <- summary(pre_neuro_pediatrics$pre_admission_pns)

  # create a list of the pre-admission summary results
  pre_admission_summary <- list(pre_cns_adults_summary = pre_cns_adults_summary,
                                pre_pns_adults_summary = pre_pns_adults_summary,
                                pre_cns_pediatrics_summary = pre_cns_pediatrics_summary,
                                pre_pns_pediatrics_summary = pre_pns_pediatrics_summary,
                                n_codes_per_patient = n_codes_per_patient)


  # add the pre cns and pns summaries
  results$pre_admission_summary <- pre_admission_summary

  # end and print total run time
  time.end.analysis=Sys.time()
  time.analysis=time.end.analysis-time.start.analysis
  print(paste('The analysis completed in', time.analysis))

  # remove unneccessary R objects
  rm(list = setdiff(ls(), c("CurrSiteId", "results")))

  # append siteId to results file
  site_results <- paste0(CurrSiteId, "_results")
  assign(site_results, results)

  if(is_docker==TRUE) {
    # save results
    save(
      list = site_results,
      file = file.path(
        getProjectOutputDirectory(),
        paste0(CurrSiteId, "_results.rda")
      )
    )
    save(
      list = site_results,
      file = file.path(
        FourCePhase2.1Data::get4ceRootDirectoryName(),
        paste0(CurrSiteId, "_results.rda")
      )
    )
    cat(
      "Result is saved in",
      file.path(
        getProjectOutputDirectory(),
        paste0(CurrSiteId, "_results.rda")
      ),
      "\nPlease submit the result file by running submitAnalysis()\n"
    )
  } else {
    cat("Saving results...")
    save(
      list = site_results,
      file = file.path(
        output_dir,
        paste0(CurrSiteId, "_results.rda")
      )
    )
    cat('analysis complete')
  }

  sink()

}
