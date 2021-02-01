
run_regression <-
  function(df, depend_var, binary = TRUE, include_race = TRUE) {
    independ_vars <- '~ neuro_post + elixhauser_score + sex + age_group'
    if (include_race)
      independ_vars <- paste(independ_vars, '+ race')

    if (binary) {
      # is dependent variable binary?
      glm(as.formula(paste(depend_var, independ_vars)),
          family = 'binomial', data = df) %>%
        summary()
    } else {
      lm(as.formula(paste(depend_var, independ_vars)), data = df) %>%
        summary()
    }
  }

run_regressions <- function(df, include_race = TRUE) {
  severe_reg_elix <-
    run_regression(df, 'severe', TRUE, include_race)

  deceased_reg_elix <-
    run_regression(df, 'deceased', TRUE, include_race)

  n_stay_reg_elix <-
    run_regression(df, 'n_stay', FALSE, include_race)

  n_readmit_reg_elix <-
    run_regression(df, 'n_readmissions', FALSE, include_race)

  readmit_reg_elix <-
    run_regression(df, 'readmitted', TRUE, include_race)

  list(
    n_stay_reg_elix = n_stay_reg_elix,
    severe_reg_elix = severe_reg_elix,
    deceased_reg_elix = deceased_reg_elix,
    n_readmit_reg_elix = n_readmit_reg_elix,
    readmit_reg_elix = readmit_reg_elix
  )
}

run_subgroup_regs <- function(df, include_race = TRUE) {
  time_severe_reg_elix <-
    run_regression(df, 'time_to_severe', FALSE, include_race)

  time_deceased_reg_elix <-
    run_regression(df, 'time_to_death', FALSE, include_race)

  time_reg_elix <-
    run_regression(df, 'time_to_first_readmission', FALSE, include_race)

  list(
    time_severe_reg_elix = time_severe_reg_elix,
    time_deceased_reg_elix = time_deceased_reg_elix,
    time_reg_elix = time_reg_elix
  )
}

run_hosps <- function(mask_thres,
                      blur_abs,
                      include_race,
                      currSiteId,
                      readmissions,
                      demo_raw,
                      obs_raw,
                      neuro_icds) {
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
  # optional
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
}
