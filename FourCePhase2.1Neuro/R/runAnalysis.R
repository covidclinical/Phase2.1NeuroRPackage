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
      death_date = if_else(death_date > lubridate::today("EST"), lubridate::NA_Date_, as.Date(death_date)),
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
            !duplicated(delta_hospitalized == 1))
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

  temp_neuro <-
    temporal_neuro(comp_readmissions, obs_raw, neuro_icds, readmissions)
  obs_first_hosp <- temp_neuro$obs_first_hosp
  propagated_codes <- temp_neuro$propagated_codes %>%
    blur_it(c("n_early_codes", "n_new_codes"), blur_abs, mask_thres) %>%
    mutate(
      prop_new_codes = if_else(n_early_codes == 0, 0, prop_new_codes),
      prob_repeated = if_else(n_early_codes == 0, 0, prob_repeated),
      prob_at_least_one_new =
        if_else(n_early_codes == 0, 0, prob_at_least_one_new)
    )

  nstay_df <- comp_readmissions %>%
    filter(first_out) %>%
    transmute(patient_num, time_to_first_discharge = days_since_admission - 1)
  pre_neuro <- obs_raw %>%
    filter(days_since_admission >= -365 & days_since_admission <= -15) %>%
    right_join(neuro_icds, by = c("concept_code" = "icd")) %>%
    filter(!is.na(patient_num)) %>%
    count(patient_num, pns_cns) %>%
    mutate(value = log(n + 1)) %>%
    pivot_wider(id_cols = patient_num, names_from = pns_cns, values_from = value, values_fill = 0) %>%
    rename(pre_admission_cns = Central,
           pre_admission_pns = Peripheral)

  pre_cns_summary <- summary(pre_neuro$pre_admission_cns)
  pre_pns_summary <- summary(pre_neuro$pre_admission_pns)

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
    left_join(readmissions, by = "patient_num") %>%
    left_join(pre_neuro, by = "patient_num") %>%
    replace_na(list(n_readmissions = 0,
                    pre_admission_cns = 0,
                    pre_admission_pns = 0)) %>%
    # refactor age and race variables
    mutate(age_group_rf = if_else(age_group == "00to02" |
                                    age_group == "03to05" |
                                    age_group == "06to11" |
                                    age_group == "12to17" |
                                  age_group == "18to25" |
                                  age_group == "26to49",
                                  "00to49", "50to80plus"),
           # age_group_rf = if_else(age_group == "18to25" |
           #                          age_group == "26to49",
           #                          "18to49", age_group_rf),
           # age_group_rf = if_else(age_group == "50to69" |
           #                          age_group == "70to79" |
           #                          age_group == "80plus",
           #                        "50to80plus", age_group_rf),
           race_rf = if_else(race == "white", "white", "NA"),
           race_rf = if_else(race == "black", "black", race_rf),
           race_rf = if_else(race == "asian" | race == "hawaiian_pacific_islander", "asian_hawaiian_pacific_islander", race_rf),
           race_rf = if_else(race == "american_indian" | race == "other", "other", race_rf))

  comorb_list <- get_elix_mat(obs_first_hosp, icd_version)

  index_scores_elix <- comorb_list$index_scores_elix
  index_scores_elix <- index_scores_elix %>%
    right_join0(select(demo_raw, patient_num), by = "patient_num")

  mapped_codes_table <- comorb_list$mapped_codes_table

  elix_mat <- cor(select(
    index_scores_elix,
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

  output_log <- file(paste0(getProjectOutputDirectory(), "output_log.txt"), open = "wt") # File name of output log
  #sink(file = file.path(getProjectOutputDirectory(), paste0(currSiteId, "_log.txt")), split = TRUE, append = FALSE)
  sink(file = output_log, type="message")

  results <- list(
    site = CurrSiteId,
    elix_mat = elix_mat,
    deviance_expl = deviance_expl,
    propagated_codes = propagated_codes,
    first_hosp_results = run_hosps(
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

  # ## reset message sink and close the file connection
  sink(type="message")

  closeAllConnections() # Close connection to log file

  ## obfuscate comorbidity table
  mapped_codes_table_obfus <- blur_it(mapped_codes_table, vars = 'n_patients', blur_abs, mask_thres)
  mapped_codes_table_obfus <- mask_it(mapped_codes_table_obfus, var = 'n_patients', blur_abs, mask_thres)

  # remove categories with 0 patients
  results$mapped_codes_table_obfus <- mapped_codes_table_obfus %>%
    filter(!n_patients == 0)

  # add the pre cns and pns summaries
  results$pre_cns_summary <- pre_cns_summary
  results$pre_pns_summary <- pre_pns_summary


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
}
