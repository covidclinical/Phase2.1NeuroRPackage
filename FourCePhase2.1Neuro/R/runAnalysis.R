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
      )
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

    ## -------------------------------------------------------------------------

    neuro_patients <- obs_raw %>%
      filter(days_since_admission >= 0) %>%
      right_join(neuro_icds, by = c('concept_code' = 'icd')) %>%
      filter(!is.na(patient_num)) %>%
      distinct(patient_num, concept_code, pns_cns) %>%
      group_by(patient_num) %>%
      mutate(nerv_sys_count = length(unique(pns_cns))) %>%
      ungroup() %>%
      mutate(neuro_type = case_when(nerv_sys_count == 2 ~ 'Both',
                                    TRUE ~ as.character(pns_cns)))

    neuro_pt_post <- unique(neuro_patients$patient_num)

    non_neuro_patients <-
      data.frame(patient_num = setdiff(demo_raw$patient_num, neuro_pt_post)) %>%
      mutate(concept_code = 'NN')

    comp_readmissions <- clin_raw %>%
      group_by(patient_num) %>%
      mutate(delta_hospitalized = diff(c(in_hospital[1], in_hospital))) %>%
      mutate(
        first_adm =
          (delta_hospitalized == -1 & !duplicated(delta_hospitalized == -1))|
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
      filter(patient_num %in% n_readms$patient_num, first_adm) %>%
      select(patient_num, delta_hospitalized, days_since_admission) %>%
      pivot_wider(names_from = delta_hospitalized,
                  values_from = days_since_admission) %>%
      mutate(time_to_first_readmission = `1` - `-1`) %>%
      select(patient_num, time_to_first_readmission) %>%
      left_join(n_readms, by = 'patient_num')

    ## -------------------------------------------------------------------------
    days_count_min_max <- obs_raw %>%
      group_by(patient_num) %>%
      summarise(
        distinct_days = n_distinct(days_since_admission),
        min_hos = min(days_since_admission),
        .groups = 'drop'
      )

    demo_processed <- demo_raw %>%
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
      left_join(days_count_min_max, by = 'patient_num') %>%
      left_join(readmissions, by = 'patient_num') %>%
      replace_na(list(n_readmissions = 0))


    ## -------------------------------------------------------------------------
    # for elixhauser
    comorb_names_elix <- get_quan_elix_names()

    # t1: earliest time point to consider comorbidities
    # t2: latest time point to consider comorbidities
    # example <- t1 = -365, and t2 = -1 will map all all codes up to
    # a year prior but before admission (admission = day 0)

    comorb_elix <- map_char_elix_codes(
      df = obs_raw,
      comorb_names = comorb_names_elix,
      t1 = -365,
      t2 = -15,
      map_type = 'elixhauser'
    )

    index_scores_elix <- comorb_elix$index_scores %>%
      rename('elixhauser_score' = van_walraven_score)
    # van Walraven is a modification of Elixhauser comorbidity measure
    # doi.org/10.1097/MLR.0b013e31819432e5
    mapped_codes_table_elix <- comorb_elix$mapped_codes_table
    elix_mat <- cor(select(index_scores_elix, - c(patient_num, elixhauser_score)))

    ## -------------------------------------------------------------------------
    nstay_df <- neuro_patients %>%
      select(patient_num, concept_code) %>%
      bind_rows(non_neuro_patients) %>%
      left_join(demo_processed, by = 'patient_num') %>%
      mutate(concept_code = fct_reorder(concept_code, n_stay)) %>%
      left_join(neuro_icds, by = c('concept_code' = 'icd'))

    icd_tables <- get_tables(
      c('no_neuro_cond', 'neuro_cond'),
      nstay_df,
      right_join0(index_scores_elix, nstay_df, by = 'patient_num'),
      comorb_names_elix,
      blur_abs,
      mask_thres,
      'concept_code'
    )[-3] # last element is not useful, remove

    ## -------------------------------------------------------------------------
    # Part 1: Binary outcome: neuro vs. non_neuro
    demo_df <- demo_processed %>%
      mutate(neuro_post = patient_num %in% neuro_pt_post %>%
               as.factor() %>%
               fct_recode(neuro_cond = "TRUE",
                          no_neuro_cond = "FALSE"))

    scores_unique <- right_join0(index_scores_elix, demo_df, by = 'patient_num')

    obfus_tables <- get_tables(
      c('no_neuro_cond', 'neuro_cond'),
      demo_df,
      scores_unique,
      comorb_names_elix,
      blur_abs,
      mask_thres
    ) %>%
      lapply(function(x) mutate(x, site = currSiteId))

    ## -------------------------------------------------------------------------
    reg_results <- run_regressions(scores_unique, include_race)
    sub_reg_results <- run_subgroup_regs(scores_unique, include_race)

    ## ----save-results---------------------------------------------------------
    binary_results <- c(obfus_tables, reg_results, sub_reg_results)

    ## -------------------------------------------------------------------------
    ### Part 2: PNS vs CNS
    neuro_types <- c('None', 'Peripheral', 'Central', 'Both')

    demo_df <- demo_processed %>%
      left_join(distinct(select(neuro_patients, patient_num, neuro_type)),
                by = 'patient_num') %>%
      replace_na(list(neuro_type = 'None')) %>%
      mutate(neuro_post = forcats::fct_relevel(neuro_type, neuro_types))

    scores_unique <- right_join0(index_scores_elix, demo_df, by = 'patient_num')

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
    reg_results <- run_regressions(scores_unique, include_race)
    sub_reg_results <- run_subgroup_regs(scores_unique, include_race)

    ## ----save-results---------------------------------------------------------
    cpns_results <- c(obfus_tables, reg_results, sub_reg_results)

    results <- list(site = currSiteId,
                    icd_tables = icd_tables,
                    elix_mat = elix_mat,
                    binary_results = binary_results,
                    cpns_results = cpns_results)

    site_results <- paste0(currSiteId, '_results')
    assign(site_results, results)
    save(list = site_results,
         file = file.path(getProjectOutputDirectory(),
                          paste0(currSiteId, '_results.rda')))
  }
