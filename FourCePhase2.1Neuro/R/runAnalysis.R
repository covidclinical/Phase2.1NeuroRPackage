#' Run analysis
#'
#' @param currSiteId Your site acronym
#' @param mask_thres Obfuscation small count mask threshold (e.g. 10)
#' @param blur_abs Absolute max of obfuscation blurring range
#' @param include_race Boolean. Whether race should be included
#' in the regression model.
#' @param getProjectOutputDirectory() Optional. Directory where resulting outputs are located.
#' Defaults to 'output'.
#' @param data_dir Optional. Directory where datasets are located.
#' Defaults to 'Input'.
#'
#' @return NULL. Result files are written out to the `getProjectOutputDirectory()` directory.
#'
#' @keywords 4CE
#' @export
#' @import dplyr
#' @importFrom stats cor glm lm median sd
#' @importFrom forcats fct_recode fct_reorder
#' @importFrom tidyr pivot_longer pivot_wider replace_na
#'
runAnalysis <-
  function(mask_thres,
           blur_abs,
           include_race = TRUE,
           data_dir = 'Input') {

    ## make sure this instance has the latest version of the quality control and data wrangling code available
    devtools::install_github("covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)

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
      )

    ## -------------------------------------------------------------------------
    neuro_patients <- obs_raw %>%
      filter(days_since_admission >= 0,
             concept_code %in% neuro_icds_10$icd) %>%
      distinct(patient_num, concept_code)

    neuro_pt_post <- unique(neuro_patients$patient_num)

    non_neuro_patients <-
      data.frame(patient_num = setdiff(demo_raw$patient_num, neuro_pt_post)) %>%
      mutate(concept_code = 'NN')

    readmissions <- clin_raw %>%
      group_by(patient_num) %>%
      mutate(delta_hospitalized = diff(c(in_hospital[1], in_hospital))) %>%
      ungroup() %>%
      filter(delta_hospitalized != 0,
             in_hospital == 1) %>%
      add_count(patient_num, name = 'n_readmissions') %>%
      arrange(desc(n_readmissions)) %>%
      select(patient_num, n_readmissions) %>%
      distinct()

    ## -------------------------------------------------------------------------
    days_count_min_max <- obs_raw %>%
      group_by(patient_num) %>%
      summarise(
        distinct_days = n_distinct(days_since_admission),
        min_hos = min(days_since_admission),
        .groups = 'drop'
      )

    demo_df <- demo_raw %>%
      mutate(
        time_to_severe = severe_date - admission_date,
        time_to_severe = ifelse(time_to_severe < 0, NA, time_to_severe),
        time_to_death = death_date - admission_date,
        time_to_death = ifelse(time_to_death < 0, NA, time_to_death),
        readmitted = patient_num %in% readmissions$patient_num,
        neuro_post = patient_num %in% neuro_pt_post %>%
          as.factor() %>%
          fct_recode(neuro_cond = "TRUE",
                     no_neuro_cond = "FALSE"),
        Survival = as.factor(deceased) %>%
          fct_recode(Alive = "0", Deceased = "1"),
        sex = as.factor(sex),
        race = as.factor(race),
        age_group = as.factor(age_group),
        Severity = as.factor(severe) %>%
          fct_recode(Severe = "1", `Non-severe` = "0"),
        n_stay = as.numeric(last_discharge_date - admission_date,
                            units = "days")
      ) %>%
      left_join(days_count_min_max, by = 'patient_num')

    ## -------------------------------------------------------------------------
    vars_to_obfs <- c("sex",
                      'age_group',
                      "race",
                      "Severity",
                      "Survival",
                      "readmitted")
    get_stats <-
      function(x)
        demo_stats(demo_df, x, blur_abs, mask_thres)
    demo_obfus_table <- lapply(vars_to_obfs, get_stats) %>%
      do.call(rbind, .)

    ## -------------------------------------------------------------------------
    nstay_obfus_table <- nstay_stats(demo_df, blur_abs, mask_thres)
    severity_obfus_table <-
      severity_stats(demo_df, blur_abs, mask_thres)
    survival_obfus_table <-
      survival_stats(demo_df, blur_abs, mask_thres)
    readmission_obfus_table <- demo_df %>%
      left_join(readmissions, by = 'patient_num') %>%
      replace_na(list(n_readmissions = 0)) %>%
      readmission_stats(blur_abs, mask_thres)

    ## -------------------------------------------------------------------------
    n_patients <- nrow(demo_raw)
    obs_raw <- obs_raw %>%
      filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"))

    ## -------------------------------------------------------------------------
    # for elixhauser
    comorb_names_elix <- get_quan_elix_names()
    comorbs_elix <- as.vector(comorb_names_elix$Abbreviation)

    # t1: earliest time point to consider comorbidities
    # t2: latest time point to consider comorbidities
    # example <- t1 = -365, and t2 = -1 will map all all codes up to a year prior but before admission (admission = day 0)

    comorb_elix <- map_char_elix_codes(
      df = obs_raw,
      comorb_names = comorb_names_elix,
      t1 = -365,
      t2 = -15,
      map_type = 'elixhauser'
    )

    ## -------------------------------------------------------------------------
    # elixhauser
    index_scores_elix <- comorb_elix$index_scores %>%
      rename('elixhauser_score' = van_walraven_score)
    # van Walraven is a modification of Elixhauser comorbidity measure
    # doi.org/10.1097/MLR.0b013e31819432e5
    mapped_codes_table_elix <- comorb_elix$mapped_codes_table
    comorb_names_elix$Abbreviation <-
      as.character(comorb_names_elix$Abbreviation)

    ## -------------------------------------------------------------------------
    n_comorbs <- colSums(index_scores_elix[, comorbs_elix])
    pos_comorbs <- names(n_comorbs[n_comorbs > 0])
    elix_mat <- cor(index_scores_elix[, pos_comorbs])
    comorb_unique <- index_scores_elix %>%
      select(patient_num, elixhauser_score) %>%
      left_join(demo_df, by = 'patient_num')

    elix_obfus_table <-
      elix_stats(comorb_unique, blur_abs, mask_thres)

    other_obfus_table <-
      bind_rows(
        nstay_obfus_table,
        readmission_obfus_table,
        severity_obfus_table,
        survival_obfus_table,
        elix_obfus_table
      )

    ## -------------------------------------------------------------------------
    scores_unique <-
      left_join(index_scores_elix, demo_df, by = 'patient_num')

    scores_neuro <- obs_raw %>%
      # 1 patient can have different code but each only counted once
      distinct(patient_num, concept_code) %>%
      left_join(neuro_icds_10, by = c('concept_code' = 'icd')) %>%
      left_join(scores_unique, by = 'patient_num') %>%
      filter(!is.na(elixhauser_score)) %>%
      mutate(
        concept_code = case_when(
          is.na(`Neurological Disease Category`) ~ 'NN',
          TRUE ~ concept_code
        ) %>%
          as.factor() %>%
          fct_reorder(-elixhauser_score),
        `Neurological Disease Category` =
          as.factor(`Neurological Disease Category`) %>%
          fct_reorder(elixhauser_score)
      ) %>%
      {
        .
      }

    elix_obfus_table1 <-
      Reduce(
        function(...)
          left_join(..., by = c("Comorbidity", "Abbreviation")),
        lapply(
          c('no_neuro_cond', 'neuro_cond'),
          list_table1,
          df = scores_unique,
          num_pats = nrow(demo_df),
          comorb_names = comorb_names_elix,
          blur_abs = blur_abs,
          mask_thres = mask_thres
        )
      ) %>%
      mutate(n_Total = n_no_neuro_cond + n_neuro_cond,
             prop_Total = n_Total / nrow(demo_raw)) %>%
      arrange(desc(n_Total))


    ## -------------------------------------------------------------------------
    severe_reg_elix <-
      glm(
        severe ~ neuro_post + elixhauser_score + sex + age_group + race,
        data = scores_unique,
        family = 'binomial'
      )

    deceased_reg_elix <-
      glm(
        deceased ~ neuro_post + elixhauser_score + sex + age_group + race,
        data = scores_unique,
        family = 'binomial'
      )

    n_stay_reg_elix <-
      lm(n_stay ~ neuro_post + elixhauser_score + race + sex + age_group,
         data = scores_unique)


    ## -------------------------------------------------------------------------
    nstay_df <- neuro_patients %>%
      bind_rows(non_neuro_patients) %>%
      left_join(demo_df, by = 'patient_num') %>%
      left_join(index_scores_elix, by = 'patient_num') %>%
      # mutate(concept_code = fct_reorder(concept_code, n_stay)) %>%
      left_join(neuro_icds_10, by = c('concept_code' = 'icd')) %>%
      mutate(
        full_icd = case_when(
          concept_code == 'NN' ~ 'No neurological condition',
          TRUE ~ paste0(`ICD-10 Description`, ' (', concept_code, ')')
        ) %>%
          as.factor() %>% fct_reorder(n_stay)
      )

    summarised_obfus_icd <- nstay_df %>%
      group_by(concept_code) %>%
      summarise(
        mean_stay = mean(n_stay),
        median_stay = median(n_stay),
        sd_stay = sd(n_stay),
        mean_elix = mean(elixhauser_score, na.rm = TRUE),
        median_elix = median(elixhauser_score, na.rm = TRUE),
        sd_elix = sd(elixhauser_score, na.rm = TRUE),
        n_patients = n(),
        prop_deceased = mean(deceased),
        prop_severe = mean(severe),
        .groups = 'drop'
      ) %>%
      blur_it('n_patients', blur_abs, mask_thres)


    ## ----save-results-----------------------------------------------------------------------------------------------------------------------
    list_results = list(
      demo_table = demo_obfus_table,
      other_obfus_table = other_obfus_table,
      elix_obfus_table1 = elix_obfus_table1,
      summarised_obfus_icd = summarised_obfus_icd
    )

    list_results <-
      lapply(list_results, function(x)
        mutate(x, site = currSiteId))
    list_results <- c(
      list_results,
      list(
        site = currSiteId,
        elix_mat = elix_mat,
        n_stay_reg_elix = n_stay_reg_elix,
        severe_reg_elix = severe_reg_elix,
        deceased_reg_elix = deceased_reg_elix
      )
    )


    ## -------------------------------------------------------------------------
    ### Part 2: PNS vs CNS

    ## -------------------------------------------------------------------------
    neuro_types <- c('None', 'Peripheral', 'Central', 'Both')

    neuro_patients <- obs_raw %>%
      filter(days_since_admission >= 0, ) %>%
      right_join(neuro_icds_10, by = c('concept_code' = 'icd')) %>%
      filter(!is.na(patient_num)) %>%
      distinct(patient_num, concept_code, pns_cns) %>%
      group_by(patient_num) %>%
      mutate(nerv_sys_count = length(unique(pns_cns))) %>%
      ungroup() %>%
      mutate(neuro_type = case_when(nerv_sys_count == 2 ~ 'Both',
                                    TRUE ~ as.character(pns_cns)))

    patients_pns_cns <- neuro_patients %>%
      distinct(patient_num, neuro_type)

    neuro_pt_post <- unique(neuro_patients$patient_num)

    non_neuro_patients <-
      data.frame(patient_num = setdiff(demo_raw$patient_num, neuro_pt_post)) %>%
      mutate(concept_code = 'NN')

    readmissions <- clin_raw %>%
      group_by(patient_num) %>%
      mutate(delta_hospitalized = diff(c(in_hospital[1], in_hospital))) %>%
      ungroup() %>%
      filter(delta_hospitalized != 0,
             in_hospital == 1) %>%
      add_count(patient_num, name = 'n_readmissions') %>%
      arrange(desc(n_readmissions)) %>%
      select(patient_num, n_readmissions) %>%
      distinct()


    ## -------------------------------------------------------------------------
    days_count_min_max <- obs_raw %>%
      group_by(patient_num) %>%
      summarise(
        distinct_days = n_distinct(days_since_admission),
        min_hos = min(days_since_admission),
        .groups = 'drop'
      )

    demo_df <- demo_raw %>%
      mutate(
        time_to_severe = severe_date - admission_date,
        time_to_severe = ifelse(time_to_severe < 0, NA, time_to_severe),
        time_to_death = death_date - admission_date,
        time_to_death = ifelse(time_to_death < 0, NA, time_to_death),
        readmitted = patient_num %in% readmissions$patient_num,
        Survival = as.factor(deceased) %>%
          fct_recode(Alive = "0", Deceased = "1"),
        sex = as.factor(sex),
        race = as.factor(race),
        age_group = as.factor(age_group),
        Severity = as.factor(severe) %>%
          fct_recode(Severe = "1", `Non-severe` = "0"),
        n_stay = as.numeric(last_discharge_date - admission_date,
                            units = "days")
      ) %>%
      left_join(select(neuro_patients, patient_num, neuro_type),
                by = 'patient_num') %>%
      replace_na(list(neuro_type = 'None')) %>%
      mutate(neuro_post = forcats::fct_relevel(neuro_type, neuro_types)) %>%
      left_join(days_count_min_max, by = 'patient_num')


    ## -------------------------------------------------------------------------
    get_stats <-
      function(x)
        demo_stats(demo_df, x, blur_abs, mask_thres)

    demo_obfus_table <- lapply(vars_to_obfs, get_stats) %>%
      do.call(rbind, .)


    ## -------------------------------------------------------------------------
    nstay_obfus_table <- nstay_stats(demo_df, blur_abs, mask_thres)
    severity_obfus_table <-
      severity_stats(demo_df, blur_abs, mask_thres)
    survival_obfus_table <-
      survival_stats(demo_df, blur_abs, mask_thres)
    readmission_obfus_table <- demo_df %>%
      left_join(readmissions, by = 'patient_num') %>%
      replace_na(list(n_readmissions = 0)) %>%
      readmission_stats(blur_abs, mask_thres)

    other_obfus_table <-
      bind_rows(
        nstay_obfus_table,
        readmission_obfus_table,
        severity_obfus_table,
        survival_obfus_table
      )

    ## -------------------------------------------------------------------------
    scores_unique <-
      left_join(index_scores_elix, demo_df, by = 'patient_num')

    scores_neuro <- obs_raw %>%
      # 1 patient can have different code but each only counted once
      distinct(patient_num, concept_code) %>%
      left_join(neuro_icds_10, by = c('concept_code' = 'icd')) %>%
      left_join(scores_unique, by = 'patient_num') %>%
      filter(!is.na(elixhauser_score)) %>%
      mutate(
        concept_code = case_when(
          is.na(`Neurological Disease Category`) ~ 'NN',
          TRUE ~ concept_code
        ) %>%
          as.factor() %>%
          fct_reorder(-elixhauser_score),
        `Neurological Disease Category` =
          as.factor(`Neurological Disease Category`) %>%
          fct_reorder(elixhauser_score)
      ) %>%
      {
        .
      }


    elix_obfus_table1 <-
      Reduce(
        function(...)
          left_join(..., by = c("Comorbidity", "Abbreviation")),
        lapply(
          neuro_types,
          list_table1,
          df = scores_unique,
          num_pats = nrow(demo_df),
          comorb_names = comorb_names_elix,
          blur_abs = blur_abs,
          mask_thres = mask_thres
        )
      ) %>%
      mutate(
        n_Total = n_None + n_Central + n_Peripheral + n_Both,
        prop_Total = n_Total / nrow(demo_raw)
      ) %>%
      arrange(desc(n_Total))

    ## -------------------------------------------------------------------------
    if (include_race){
      severe_reg_elix <-
        glm(
          severe ~ neuro_post + elixhauser_score + sex + age_group + race,
          data = scores_unique,
          family = 'binomial'
        )
      deceased_reg_elix <-
        glm(
          deceased ~ neuro_post + elixhauser_score + sex + age_group + race,
          data = scores_unique,
          family = 'binomial'
        )

      n_stay_reg_elix <-
        lm(n_stay ~ neuro_post + elixhauser_score + race + sex + age_group,
           data = scores_unique)
    } else {
      severe_reg_elix <-
        glm(
          severe ~ neuro_post + elixhauser_score + sex + age_group,
          data = scores_unique,
          family = 'binomial'
        )
      deceased_reg_elix <-
        glm(
          deceased ~ neuro_post + elixhauser_score + sex + age_group,
          data = scores_unique,
          family = 'binomial'
        )

      n_stay_reg_elix <-
        lm(n_stay ~ neuro_post + elixhauser_score + sex + age_group,
           data = scores_unique)
    }


    ## ----save-results---------------------------------------------------------
    list_results_cpns = list(
      demo_table = demo_obfus_table,
      other_obfus_table = other_obfus_table,
      elix_obfus_table1 = elix_obfus_table1,
      summarised_obfus_icd = summarised_obfus_icd
    )

    list_results_cpns <-
      lapply(list_results_cpns, function(x)
        mutate(x, site = currSiteId))
    list_results_cpns <- c(
      list_results_cpns,
      list(
        site = currSiteId,
        elix_mat = elix_mat,
        n_stay_reg_elix = n_stay_reg_elix,
        severe_reg_elix = severe_reg_elix,
        deceased_reg_elix = deceased_reg_elix
      )
    )

    results <- list(list_results = list_results,
                    list_results_cpns = list_results_cpns)
    site_results <- paste0(currSiteId, '_results')
    assign(site_results, results)
    save(list = site_results,
         file = file.path(getProjectOutputDirectory(),
                          paste0(currSiteId, '_results.rda')))


  }
