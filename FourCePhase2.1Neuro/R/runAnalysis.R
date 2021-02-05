#' Run analysis
#'
#' @param mask_thres Obfuscation small count mask threshold (e.g. 10)
#' @param blur_abs Absolute max of obfuscation blurring range
#' @param icd_version Version of ICD code the site uses. Must be EITHER 9 or 10.
#' @param include_race Boolean. Whether race should be included
#' in the regression model.
#' @param data_dir Optional. Directory where datasets are located.
#' Defaults to 'Input'.
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
runAnalysis <-
  function(mask_thres,
           blur_abs,
           icd_version = 10,
           include_race = TRUE,
           data_dir = 'Input') {

    ## make sure this instance has the latest version of the quality control and data wrangling code available
    # remotes::install_github("covidclinical/Phase2.1DataRPackage@77d32fe", subdir="FourCePhase2.1Data", upgrade=FALSE)

    ## get the site identifier assocaited with the files stored in the /4ceData/Input directory that
    ## is mounted to the container
    currSiteId <- FourCePhase2.1Data::getSiteId()

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
        na = '1900-01-01'
      ) %>%
      mutate(across(ends_with('_date') &
                      where(is.character), lubridate::mdy)) %>%
      # mutate(
      #   last_discharge_date = if_else(
      #     !is.na(death_date) & death_date < last_discharge_date,
      #     death_date,
      #     last_discharge_date
      #   )
      # ) %>%
      {.}

    obs_raw <-
      readr::read_csv(
        file.path(data_dir, 'LocalPatientObservations.csv'),
        col_types = list(patient_num = readr::col_character())
      ) %>%
      filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"))

    if (icd_version == 9){
      neuro_icds <- neuro_icds_9
    } else {
      neuro_icds <- neuro_icds_10
    }

    comp_readmissions <- clin_raw %>%
      group_by(patient_num) %>%
      mutate(delta_hospitalized = diff(c(in_hospital[1], in_hospital))) %>%
      mutate(
        first_out =
          (delta_hospitalized == -1 & !duplicated(delta_hospitalized == -1)),
        first_change =
          first_out|
          (delta_hospitalized == 1 & !duplicated(delta_hospitalized == 1))) %>%
      ungroup() %>%
      {.}

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

    obs_first_hosp <- comp_readmissions %>%
      filter(first_out) %>%
      # days since admission the patient is out of hospital
      transmute(patient_num, dsa = days_since_admission) %>%
      right_join(obs_raw, by = 'patient_num') %>%
      filter(days_since_admission < dsa) %>%
      select(-dsa)

    results <- list(
      all_hosp_results = run_hosps(
        mask_thres,
        blur_abs,
        include_race,
        currSiteId,
        readmissions,
        demo_raw,
        obs_raw,
        neuro_icds
      ),
      first_hosp_results = run_hosps(
        mask_thres,
        blur_abs,
        include_race,
        currSiteId,
        readmissions,
        demo_raw,
        obs_first_hosp,
        neuro_icds
      )
    )

    site_results <- paste0(currSiteId, '_results')
    assign(site_results, results)
    save(list = site_results,
         file = file.path(getProjectOutputDirectory(),
                          paste0(currSiteId, '_results.rda')))
    cat('Result is saved in',
        file.path(getProjectOutputDirectory(),
                  paste0(currSiteId, '_results.rda')),
        '\n')
  }
