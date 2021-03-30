
run_regression <-
  function(df, depend_var, ind_vars, binary = TRUE) {
    if (length(unique(df[, depend_var, drop = T])) <= 1)
      return(NULL)

    independ_vars <- paste(ind_vars, collapse = ' + ')

    if (binary) {
      # is dependent variable binary?
      output <- tryCatch(
        {glm(as.formula(paste(depend_var, '~', independ_vars)),
              family = 'binomial', data = df) %>%
            summary()},
        error = function(cond) {
          message(paste("Error when regressing", depend_var))
          message("Original error message:")
          message(cond)
          message('Skipping for now...')
          return(NULL) # return NA in case of error
        }
      )
      if (!is.null(output)){
        output$deviance.resid <- NULL
        output$na.action <- NULL
        output$terms <- NULL
      }

    } else {

      output <- tryCatch(
        {lm(as.formula(paste(depend_var, '~', independ_vars)), data = df) %>%
            summary()},
        error = function(cond) {
          message(paste("Error when regressing", depend_var))
          message("Original error message:")
          message(cond)
          message('Skipping for now...')
          return(NULL) # return NA in case of error
        }
      )
      if (!is.null(output)){
        output$residuals <- NULL
        output$na.action <- NULL
        output$terms <- NULL
      }
    }

    output
  }

get_ind_vars <- function(df, include_race){
  unique_cols <- apply(df, 2, function(x) length(unique(x)))

  ind_vars <- setdiff(
    c('neuro_post', 'sex', 'age_group',
      paste0('.fittedPC', 1:10)),
    names(unique_cols)[unique_cols == 1])

  if (include_race)
    ind_vars <- c(ind_vars, 'race')
  ind_vars
}

run_regressions <- function(df, include_race = TRUE) {
  ind_vars <- get_ind_vars(df, include_race)

  severe_reg_elix <-
    run_regression(df, 'severe', ind_vars, TRUE)

  deceased_reg_elix <-
    run_regression(df, 'deceased', ind_vars, TRUE)

  n_stay_reg_elix <-
    run_regression(df, 'n_stay', ind_vars, FALSE)

  n_readmit_reg_elix <-
    run_regression(df, 'n_readmissions', ind_vars, FALSE)

  readmit_reg_elix <-
    run_regression(df, 'readmitted', ind_vars, TRUE)

  list(
    n_stay_reg_elix = n_stay_reg_elix,
    severe_reg_elix = severe_reg_elix,
    deceased_reg_elix = deceased_reg_elix,
    n_readmit_reg_elix = n_readmit_reg_elix,
    readmit_reg_elix = readmit_reg_elix
  )
}

run_subgroup_regs <- function(df, include_race = TRUE) {
  ind_vars <- get_ind_vars(df, include_race)

  time_severe_reg_elix <-
    run_regression(df, 'time_to_severe', ind_vars, FALSE)

  time_deceased_reg_elix <-
    run_regression(df, 'time_to_death', ind_vars, FALSE)

  time_readmit_reg_elix <-
    run_regression(df, 'time_to_first_readmission', ind_vars, FALSE)

  list(
    time_severe_reg_elix = time_severe_reg_elix,
    time_deceased_reg_elix = time_deceased_reg_elix,
    time_readmit_reg_elix = time_readmit_reg_elix
  )
}

run_hosps <- function(mask_thres,
                      blur_abs,
                      include_race,
                      currSiteId,
                      readmissions,
                      demo_processed,
                      obs_raw,
                      neuro_icds,
                      index_scores_elix,
                      pca_covariates) {
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
    data.frame(patient_num = setdiff(demo_processed$patient_num, neuro_pt_post)) %>%
    mutate(concept_code = 'NN')

  ## -------------------------------------------------------------------------
  # days_count_min_max <- obs_raw %>%
  #   group_by(patient_num) %>%
  #   summarise(
  #     distinct_days = n_distinct(days_since_admission),
  #     min_hos = min(days_since_admission),
  #     .groups = 'drop'
  #   )


  ## -------------------------------------------------------------------------
  # optional
  nstay_df <- neuro_patients %>%
    select(patient_num, concept_code) %>%
    bind_rows(non_neuro_patients) %>%
    left_join(demo_processed, by = 'patient_num') %>%
    mutate(concept_code = fct_reorder(concept_code, n_stay)) %>%
    left_join(neuro_icds, by = c('concept_code' = 'icd'))

  comorb_names_elix <- get_quan_elix_names()

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

  scores_unique <- index_scores_elix %>%
    right_join0(demo_df, by = 'patient_num') %>%
    left_join(pca_covariates, by = 'patient_num')

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

  scores_unique <- index_scores_elix %>%
    right_join0(demo_df, by = 'patient_num') %>%
    left_join(pca_covariates, by = 'patient_num')

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

  results <- list(icd_tables = icd_tables,
                  binary_results = binary_results,
                  cpns_results = cpns_results)
}

get_elix_mat <- function(obs_raw, icd_version, t1 = -365, t2 = -15, map_type = 'elixhauser'){
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
    icd_version = icd_version,
    t1 = t1,
    t2 = t2,
    map_type = map_type
  )

  index_scores_elix <- comorb_elix$index_scores %>%
    rename('elixhauser_score' = van_walraven_score)
  # van Walraven is a modification of Elixhauser comorbidity measure
  # doi.org/10.1097/MLR.0b013e31819432e5
  # mapped_codes_table_elix <- comorb_elix$mapped_codes_table

  index_scores_elix
}

temporal_neuro <- function(comp_readmissions, obs_raw, neuro_icds, readmissions){
  obs_first_hosp <- comp_readmissions %>%
    filter(first_out) %>%
    # days since admission the patient is out of hospital
    transmute(patient_num, dsa = days_since_admission) %>%
    right_join(obs_raw, by = 'patient_num') %>%
    filter(days_since_admission < dsa) %>%
    select(-dsa)

  first_neuro_conds <- obs_first_hosp %>%
    filter(days_since_admission >= 0) %>%
    right_join(neuro_icds, by = c('concept_code' = 'icd')) %>%
    distinct(patient_num, early_code = concept_code) %>%
    group_by(patient_num) %>%
    summarise_all(list(~ list(.))) %>%
    {.}

  obs_later_hosp <- comp_readmissions %>%
    filter(first_out) %>%
    # days since admission the patient is out of hospital
    transmute(patient_num,
              dsa = days_since_admission) %>%
    right_join(obs_raw, by = 'patient_num') %>%
    filter(days_since_admission >= dsa) %>%
    right_join(neuro_icds, by = c('concept_code' = 'icd')) %>%
    distinct(patient_num, later_code = concept_code) %>%
    group_by(patient_num) %>%
    summarise_all(list( ~ list(.))) %>%
    ungroup() %>%
    right_join(first_neuro_conds, by = 'patient_num') %>%
    mutate(
      n_new_code = purrr::map2(later_code, early_code, ~ length(setdiff(.x, .y))),
      repeated_code = purrr::map2(later_code, early_code, intersect),
      readmitted = patient_num %in% readmissions$patient_num,
      prop_new_code = purrr::map2(n_new_code, later_code, ~ .x/length(.y))
    ) %>%
    select(- patient_num)

  new_codes <- obs_later_hosp %>%
    filter(readmitted) %>%
    tidyr::unnest(c(early_code, n_new_code, prop_new_code)) %>%
    mutate_at(vars(prop_new_code),
              ~ replace(., is.nan(.), 0)) %>%
    group_by(early_code) %>%
    summarise(
      n_early_codes = n(),
      n_new_codes = sum(n_new_code),
      n_no_new_codes = sum(n_new_code == 0),
      at_least_one_new_code = sum(n_new_code > 0),
      prop_new_codes = sum(prop_new_code)
    )

  propagated_codes <- obs_later_hosp %>%
    filter(readmitted) %>%
    tidyr::unnest(repeated_code) %>%
    pull(repeated_code) %>%
    table() %>%
    data.frame() %>%
    `colnames<-`(c('early_code', 'repeated')) %>%
    right_join0(new_codes, by = 'early_code') %>%
    transmute(
      early_code,
      n_early_codes,
      n_new_codes,
      prop_new_codes,
      prob_repeated = repeated / n_early_codes,
      prob_at_least_one_new = at_least_one_new_code / n_early_codes
    )

  list(obs_first_hosp = obs_first_hosp,
       propagated_codes = propagated_codes)

}
