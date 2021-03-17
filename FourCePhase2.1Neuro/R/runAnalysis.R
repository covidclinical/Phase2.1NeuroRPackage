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
  ## make sure this instance has the latest version of the quality control and data wrangling code available
  # remotes::install_github("covidclinical/Phase2.1DataRPackage@77d32fe", subdir="FourCePhase2.1Data", upgrade=FALSE)

  ## get the site identifier assocaited with the files stored in the /4ceData/Input directory that
  ## is mounted to the container
  data_dir <- FourCePhase2.1Data::getInputDataDirectoryName()
  currSiteId <- FourCePhase2.1Data::getSiteId()
  site_specs <- site_params %>%
    dplyr::filter(tolower(siteid) == tolower(currSiteId))

  mask_thres <- site_specs$mask_thres
  blur_abs <- site_specs$blur_abs
  icd_version <- site_specs$icd_version
  include_race <- site_specs$include_race

  ## run the quality control
  FourCePhase2.1Data::runQC(currSiteId)

  # count mask threshold
  # absolute max of blurring range
  clin_raw <-
    readr::read_csv(
      file.path(data_dir, 'LocalPatientClinicalCourse.csv'),
      col_types = list(patient_num = readr::col_character())
    )
  demo_raw <-
    readr::read_csv(
      file.path(data_dir, 'LocalPatientSummary.csv'),
      col_types = list(patient_num = readr::col_character()),
      na = c('1900-01-01', '1/1/1900')
    ) %>%
    # mutate_at(vars(which(sapply(., is.character) & names(contains('_date')))),
    # lubridate::mdy) %>%
    mutate(across(
      ends_with('_date') &
        tidywhere(is.character),
      lubridate::mdy
    )) %>%
    mutate(
      last_discharge_date = if_else(
        !is.na(death_date) & death_date < last_discharge_date,
        death_date,
        last_discharge_date
      )
    )

  obs_raw <-
    readr::read_csv(
      file.path(data_dir, 'LocalPatientObservations.csv'),
      col_types = list(patient_num = readr::col_character())
    ) %>%
    filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"))

  if (icd_version == 9) {
    neuro_icds <- neuro_icds_9
  } else {
    neuro_icds <- neuro_icds_10
  }

  #####
  # start analysis
  comp_readmissions <- clin_raw %>%
    group_by(patient_num) %>%
    mutate(delta_hospitalized = diff(c(in_hospital[1], in_hospital))) %>%
    mutate(
      first_out =
        (
          delta_hospitalized == -1 & !duplicated(delta_hospitalized == -1)
        ),
      first_change =
        first_out |
        (delta_hospitalized == 1 &
           !duplicated(delta_hospitalized == 1))
    ) %>%
    ungroup()

  n_readms <- comp_readmissions %>%
    filter(delta_hospitalized != 0,
           in_hospital == 1) %>%
    add_count(patient_num, name = 'n_readmissions') %>%
    arrange(desc(n_readmissions)) %>%
    select(patient_num, n_readmissions) %>%
    distinct()

  readmissions <- comp_readmissions %>%
    filter(patient_num %in% n_readms$patient_num, first_change) %>%
    select(patient_num, delta_hospitalized, days_since_admission) %>%
    pivot_wider(names_from = delta_hospitalized,
                values_from = days_since_admission) %>%
    mutate(time_to_first_readmission = `1` - `-1`) %>%
    select(patient_num, time_to_first_readmission) %>%
    left_join(n_readms, by = 'patient_num')

  temp_neuro <-
    temporal_neuro(comp_readmissions, obs_raw, neuro_icds, readmissions)
  obs_first_hosp <- temp_neuro$obs_first_hosp
  propagated_codes <- temp_neuro$propagated_codes %>%
    blur_it(c('n_early_codes', 'n_new_codes'), blur_abs, mask_thres) %>%
    mutate(
      prop_new_codes = if_else(n_early_codes == 0, 0, prop_new_codes),
      prob_repeated = if_else(n_early_codes == 0, 0, prob_repeated),
      prob_at_least_one_new =
        if_else(n_early_codes == 0, 0, prob_at_least_one_new)
    )

  nstay_df <- comp_readmissions %>%
    filter(first_out) %>%
    select(patient_num, n_stay = days_since_admission)

  demo_processed_first <- demo_raw %>%
    left_join(nstay_df, by = 'patient_num') %>%
    mutate(
      time_to_severe = as.numeric(severe_date - admission_date, 'days'),
      time_to_death = as.numeric(death_date - admission_date, 'days'),
      time_to_severe = if_else(
        time_to_severe < 0 |
          time_to_severe > n_stay,
        NA_real_,
        time_to_severe
      ),
      time_to_death = if_else(
        time_to_death < 0 |
          time_to_death > n_stay,
        NA_real_,
        time_to_death
      ),
      severe =  if_else(is.na(time_to_severe), 0, 1),
      deceased =  if_else(is.na(time_to_death), 0, 1),
      readmitted = patient_num %in% readmissions$patient_num,
      sex = as.factor(sex),
      race = as.factor(race),
      age_group = as.factor(age_group),
      Severity = as.factor(severe) %>%
        fct_recode(Severe = "1", `Non-severe` = "0"),
      Survival = as.factor(deceased) %>%
        fct_recode(Alive = "0", Deceased = "1")
    ) %>%
    # left_join(days_count_min_max, by = 'patient_num') %>%
    left_join(readmissions, by = 'patient_num') %>%
    replace_na(list(n_readmissions = 0))

  demo_processed_all <- demo_raw %>%
    mutate(
      time_to_severe = severe_date - admission_date,
      time_to_severe = ifelse(time_to_severe < 0, NA, time_to_severe),
      time_to_death = death_date - admission_date,
      time_to_death = ifelse(time_to_death < 0, NA, time_to_death),
      readmitted = patient_num %in% readmissions$patient_num,
      sex = as.factor(sex),
      race = as.factor(race),
      age_group = as.factor(age_group),
      Severity = as.factor(severe) %>%
        fct_recode(Severe = "1", `Non-severe` = "0"),
      Survival = as.factor(deceased) %>%
        fct_recode(Alive = "0", Deceased = "1"),
      n_stay = as.numeric(last_discharge_date - admission_date,
                          units = "days")
    ) %>%
    # left_join(days_count_min_max, by = 'patient_num') %>%
    left_join(readmissions, by = 'patient_num') %>%
    replace_na(list(n_readmissions = 0))

  index_scores_elix <- get_elix_mat(obs_first_hosp, icd_version)
  elix_mat <-
    cor(select(index_scores_elix,-c(patient_num, elixhauser_score)))

  results <- list(
    site = currSiteId,
    elix_mat = elix_mat,
    propagated_codes = propagated_codes,
    all_hosp_results = run_hosps(
      mask_thres,
      blur_abs,
      include_race,
      currSiteId,
      readmissions,
      demo_processed_all,
      obs_first_hosp,
      neuro_icds,
      index_scores_elix
    ),
    first_hosp_results = run_hosps(
      mask_thres,
      blur_abs,
      include_race,
      currSiteId,
      readmissions,
      demo_processed_first,
      obs_first_hosp,
      neuro_icds,
      index_scores_elix
    )
  )

  rm(list = setdiff(ls(), c('currSiteId', 'results')))

  site_results <- paste0(currSiteId, '_results')
  assign(site_results, results)
  save(list = site_results,
       file = file.path(
         getProjectOutputDirectory(),
         paste0(currSiteId, '_results.rda')
       ))
  save(list = site_results,
       file = paste0(currSiteId, '_results.rda'))
  cat(
    'Result is saved in',
    file.path(
      getProjectOutputDirectory(),
      paste0(currSiteId, '_results.rda')
    ),
    '\nPlease submit the result file by running submitAnalysis()\n'
  )
}
