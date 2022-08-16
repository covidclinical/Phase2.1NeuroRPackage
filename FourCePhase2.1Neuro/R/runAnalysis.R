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
runAnalysis <- function() {

  print('Set timer - this analysis will take sometime to run')
  time.start.analysis=Sys.time()

  set.seed(446) # for obfuscation posterity

  select <- dplyr::select
  filter <- dplyr::filter

  ## make sure this instance has the latest version of the quality control and data wrangling code available
  # remotes::install_github("covidclinical/Phase2.1DataRPackage@77d32fe", subdir="FourCePhase2.1Data", upgrade=FALSE)

  ## get the site identifier assocaited with the files stored in the /4ceData/Input directory that
  ## is mounted to the container
  data_dir <- FourCePhase2.1Data::getInputDataDirectoryName()
  currSiteId <- FourCePhase2.1Data::getSiteId()
  CurrSiteId <- toupper(currSiteId)
  site_specs <- site_params %>%
    dplyr::filter(tolower(siteid) == tolower(currSiteId))

  mask_thres <- site_specs$mask_thres
  blur_abs <- site_specs$blur_abs
  icd_version <- site_specs$icd_version
  include_race <- site_specs$include_race

  ## run the quality control
  remotes::install_github(
    "https://github.com/covidclinical/Phase2.1DataRPackage",
    subdir = "FourCePhase2.1Data", upgrade = FALSE
  )

  FourCePhase2.1Data::runQC(currSiteId)

  # count mask threshold
  # absolute max of blurring range
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

  ## additional formatting of data for modeling
  nstay_df <- neuro_patients %>%
    data.frame() %>%
    select(patient_num, concept_code) %>%
    rbind(non_neuro_patients) %>%
    left_join(demo_processed_first, by = "patient_num") %>%
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

  # perform the run_hosps to compute our primary analysis and save all results to a list object
  results <- list(
    site = CurrSiteId,
    elix_mat = elix_mat,
    elix_mat_cns = elix_mat_cns,
    elix_mat_pns = elix_mat_cns,
    deviance_expl = deviance_expl,
    first_hosp_results = run_hosps(
      neuro_patients,
      neuro_pt_post,
      non_neuro_patients,
      both_pts,
      both,
      mask_thres,
      blur_abs,
      include_race,
      currSiteId,
      readmissions,
      demo_processed_first,
      obs_first_hosp,
      neuro_icds,
      index_scores_elix,
      pca_covariates
    )
  )

  ## obfuscate comorbidity table
  mapped_codes_table_obfus <- blur_it(mapped_codes_table, vars = 'n_patients', blur_abs, mask_thres)
  mapped_codes_table_obfus <- mask_it(mapped_codes_table_obfus, var = 'n_patients', blur_abs, mask_thres)

  # remove categories with 0 patients
  results$mapped_codes_table_obfus <- mapped_codes_table_obfus %>%
    filter(!n_patients == 0)

  # remove "Both" patients from pre_neuro
  pre_neuro <- pre_neuro %>%
    filter(!patient_num %in% both_pts$patient_num)

  pre_cns_summary <- summary(pre_neuro$pre_admission_cns)
  pre_pns_summary <- summary(pre_neuro$pre_admission_pns)

  # add the pre cns and pns summaries
  results$pre_cns_summary <- pre_cns_summary
  results$pre_pns_summary <- pre_pns_summary

  # add n_codes_per_patients dataframe
  results$n_codes_per_patient <- n_codes_per_patient

 # timing for testing when we are not submitting results
 # time.end.analysis=Sys.time()
 # time.analysis=time.end.analysis-time.start.analysis
 # print(paste('The analysis completed in', time.analysis))


  rm(list = setdiff(ls(), c("CurrSiteId", "results")))

  site_results <- paste0(CurrSiteId, "_results")
  assign(site_results, results)
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

  time.end.analysis=Sys.time()
  time.analysis=time.end.analysis-time.start.analysis
  print(paste('The analysis completed in', time.analysis))

}
