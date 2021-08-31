
run_regression <-
  function(df, depend_var, ind_vars, binary = TRUE) {
    if (length(unique(df[, depend_var, drop = T])) <= 1) {
      return(NULL)
    }

    independ_vars <- paste(ind_vars, collapse = " + ")

    if (binary) {
      # is dependent variable binary?
      output <- tryCatch(
        {
          glm(as.formula(paste(depend_var, "~", independ_vars)),
            family = "binomial", data = df
          ) %>%
            summary()
        },
        error = function(cond) {
          message(paste("Error when regressing", depend_var))
          message("Original error message:")
          message(cond)
          message("Skipping for now...")
          return(NULL) # return NA in case of error
        }
      )
      if (!is.null(output)) {
        output$deviance.resid <- NULL
        output$na.action <- NULL
        output$terms <- NULL
      }
    } else {
      output <- tryCatch(
        {
          lm(as.formula(paste(depend_var, "~", independ_vars)), data = df) %>%
            summary()
        },
        error = function(cond) {
          message(paste("Error when regressing", depend_var))
          message("Original error message:")
          message(cond)
          message("Skipping for now...")
          return(NULL) # return NA in case of error
        }
      )
      if (!is.null(output)) {
        output$residuals <- NULL
        output$na.action <- NULL
        output$terms <- NULL
      }
    }

    output
  }

get_ind_vars <- function(df, include_race) {
  unique_cols <- apply(df, 2, function(x) length(unique(x)))

  ind_vars <- setdiff(
    c(
      "neuro_post", "sex", "age_group",
      "pre_admission_cns", "pre_admission_pns",
      paste0(".fittedPC", 1:10)
    ),
    names(unique_cols)[unique_cols == 1]
  )

  if (include_race) {
    ind_vars <- c(ind_vars, "race")
  }
  ind_vars
}

run_regressions <- function(df, include_race = TRUE) {
  ind_vars <- get_ind_vars(df, include_race)

  n_readmit_reg_elix <-
    run_regression(df, "n_readmissions", ind_vars, FALSE)

  n_readmit_reg_elix$fstatistic$dendf <- NULL
  n_readmit_reg_elix$df <- NULL

  list(n_readmit_reg_elix = n_readmit_reg_elix)
}

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
          survival::coxph(data = surv_df, id = patient_num) %>%
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

  event_table_obfs <- tryCatch(
    {
      message("generating event_tables for cox model")
      event_table = data.frame(neuro_status = output$life$strata, time = output$life$time, n.risk = output$life$n.risk, n.event = output$life$n.event, n.censor = output$life$n.censor)

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

  if (length(unique(surv_df$neuro_post)) == 4) {

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
      newdata=data.frame(rbind(c(0,0,0,meancovariate[-(1:3)]),
                               c(1,0,0,meancovariate[-(1:3)]),
                               c(0,1,0,meancovariate[-(1:3)]),
                               c(0,0,1,meancovariate[-(1:3)])))
      colnames(newdata)=names(meancovariate)
      survout=survival::survfit(cox,newdata)
      message("generating event_tables for adjusted survival curves")
      event_table1 = data.frame(time = survout[1]$time, n.risk = survout[1]$n.risk, n.event = survout[1]$n.event, n.censor = survout[1]$n.censor, neuro_status = "neuro_postNone")
      event_table2 = data.frame(time = survout[2]$time, n.risk = survout[2]$n.risk, n.event = survout[2]$n.event, n.censor = survout[2]$n.censor, neuro_status = "neuro_postPeripheral")
      event_table3 = data.frame(time = survout[3]$time, n.risk = survout[3]$n.risk, n.event = survout[3]$n.event, n.censor = survout[3]$n.censor, neuro_status = "neuro_postCentral")

      event_table_surv_adjust <- rbind(event_table1, event_table2, event_table3)

      # mask and blur for obfuscation
      message("blurring event_tables for adjusted survival curves")
      event_table_obfs <- blur_it(event_table_surv_adjust, vars = c("n.risk", "n.event", "n.censor"), blur_abs, mask_thres)

      cox <- cox %>% summary()

      average_survival =list("cox"=cox,"survf"=survout, "event_table_obfs" = event_table_obfs)

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
    right_join(neuro_icds, by = c("concept_code" = "icd")) %>%
    filter(!is.na(patient_num)) %>%
    distinct(patient_num, concept_code, pns_cns) %>%
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
    filter(!neuro_type == "Both")

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
    data.frame(patient_num = setdiff(demo_processed$patient_num, neuro_pt_post)) %>%
    mutate(concept_code = "NN") %>%
    # remove both patients
    filter(!patient_num %in% both_pts$patient_num)


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

  ## -------------------------------------------------------------------------
  # Part 1: Binary outcome: neuro vs. non_neuro
  demo_df <- demo_processed %>%
    filter(!patient_num %in% both_pts$patient_num) %>%
    mutate(neuro_post = patient_num %in% neuro_pt_post %>%
      as.factor() %>%
      fct_recode(
        neuro_cond = "TRUE",
        no_neuro_cond = "FALSE"
      ))

  scores_unique <- index_scores_elix %>%
    right_join0(demo_df, by = "patient_num") %>%
    left_join(pca_covariates, by = "patient_num")

  obfus_tables <- get_tables(
    c("no_neuro_cond", "neuro_cond"),
    demo_df,
    scores_unique,
    comorb_names_elix,
    blur_abs,
    mask_thres
  ) %>%
    lapply(function(x) mutate(x, site = currSiteId))


  ## -------------------------------------------------------------------------
  reg_results <- run_regressions(scores_unique, include_race)
  sub_reg_results <- run_coxregressions(scores_unique, include_race, blur_abs, mask_thres)

  ## ----save-results---------------------------------------------------------
  binary_results <- c(obfus_tables, reg_results, sub_reg_results)

  ## -------------------------------------------------------------------------
  ### Part 2: PNS vs CNS
  neuro_types <- c("None", "Peripheral", "Central")

  demo_df <- demo_processed %>%
    filter(!patient_num %in% both_pts$patient_num) %>%
    left_join(distinct(select(neuro_patients, patient_num, neuro_type)),
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
        filter(time_to_first_discharge == 0)
    } else if (time_to_outcome == "time_to_last_discharge") {
      demo_subset_df <- df %>%
        filter(time_to_last_discharge == 0)
    }

    scores_unique <- index_scores_elix %>%
      right_join0(demo_subset_df, by = "patient_num") %>%
      left_join(pca_covariates, by = "patient_num")

    obfus_tables <- get_tables(
      neuro_types,
      demo_subset_df,
      scores_unique,
      comorb_names_elix,
      blur_abs,
      mask_thres
    ) %>%
      lapply(function(x) mutate(x, site = currSiteId))

    return(obfus_tables)
  }

  severe_adm <- surv_exclude_pts(demo_df, "time_to_severe")
  death_adm <- surv_exclude_pts(demo_df, "time_to_death")
  first_adm <- surv_exclude_pts(demo_df, "time_to_first_discharge")
  last_adm <- surv_exclude_pts(demo_df, "time_to_last_discharge")

  ## -------------------------------------------------------------------------
  reg_results <- run_regressions(scores_unique, include_race)
  sub_reg_results <- run_coxregressions(scores_unique, include_race, blur_abs, mask_thres)

  ## ----save-results---------------------------------------------------------
  cpns_results <- c(obfus_tables, reg_results, sub_reg_results)

  results <- list(
    icd_tables = icd_tables,
    binary_results = binary_results,
    cpns_results = cpns_results,
    severe_adm = severe_adm,
    death_adm = death_adm,
    first_adm = first_adm,
    last_adm = last_adm,
    n_codes_per_patient = n_codes_per_patient,
    both
  )

  return(results)
}

get_elix_mat <- function(obs_raw, icd_version, t1 = -365, t2 = -15, map_type = "elixhauser") {
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
    rename("elixhauser_score" = van_walraven_score)
  # van Walraven is a modification of Elixhauser comorbidity measure
  # doi.org/10.1097/MLR.0b013e31819432e5

  mapped_codes_table <- comorb_elix$mapped_codes_table

  comorb_list <- list(index_scores_elix = index_scores_elix,
                      mapped_codes_table = mapped_codes_table)

  return(comorb_list)
}

temporal_neuro <- function(comp_readmissions, obs_raw, neuro_icds, readmissions) {
  obs_first_hosp <- comp_readmissions %>%
    filter(first_out) %>%
    # days since admission the patient is out of hospital
    transmute(patient_num, dsa = days_since_admission) %>%
    right_join(obs_raw, by = "patient_num") %>%
    filter(days_since_admission < dsa) %>%
    select(-dsa)

  first_neuro_conds <- obs_first_hosp %>%
    filter(days_since_admission >= 0) %>%
    right_join(neuro_icds, by = c("concept_code" = "icd")) %>%
    distinct(patient_num, early_code = concept_code) %>%
    group_by(patient_num) %>%
    summarise_all(list(~ list(.)))

  obs_later_hosp <- comp_readmissions %>%
    filter(first_out) %>%
    # days since admission the patient is out of hospital
    transmute(patient_num,
      dsa = days_since_admission
    ) %>%
    right_join(obs_raw, by = "patient_num") %>%
    filter(days_since_admission >= dsa) %>%
    right_join(neuro_icds, by = c("concept_code" = "icd")) %>%
    distinct(patient_num, later_code = concept_code) %>%
    group_by(patient_num) %>%
    summarise_all(list(~ list(.))) %>%
    ungroup() %>%
    right_join(first_neuro_conds, by = "patient_num")

  #8.31.2021 - will remove the functions for propagated codes for now. I don't think this is excluding our "Both" patients
  #%>%
  #   mutate(
  #     n_new_code = purrr::map2(later_code, early_code, ~ length(setdiff(.x, .y))),
  #     repeated_code = purrr::map2(later_code, early_code, intersect),
  #     readmitted = patient_num %in% readmissions$patient_num,
  #     prop_new_code = purrr::map2(n_new_code, later_code, ~ .x / length(.y))
  #   ) %>%
  #   select(-patient_num)
  #
  # new_codes <- obs_later_hosp %>%
  #   filter(readmitted) %>%
  #   tidyr::unnest(c(early_code, n_new_code, prop_new_code)) %>%
  #   mutate_at(
  #     vars(prop_new_code),
  #     ~ replace(., is.nan(.), 0)
  #   ) %>%
  #   group_by(early_code) %>%
  #   summarise(
  #     n_early_codes = n(),
  #     n_new_codes = sum(n_new_code),
  #     n_no_new_codes = sum(n_new_code == 0),
  #     at_least_one_new_code = sum(n_new_code > 0),
  #     prop_new_codes = sum(prop_new_code)
  #   )

  # propagated_codes <- NULL
  #
  # if (sum(obs_later_hosp$readmitted) > 0) {
  #   propagated_codes <- obs_later_hosp %>%
  #     filter(readmitted) %>%
  #     tidyr::unnest(repeated_code) %>%
  #     pull(repeated_code) %>%
  #     table() %>%
  #     data.frame() %>%
  #     `colnames<-`(c("early_code", "repeated")) %>%
  #     right_join0(new_codes, by = "early_code") %>%
  #     transmute(
  #       early_code,
  #       n_early_codes,
  #       n_new_codes,
  #       prop_new_codes,
  #       prob_repeated = repeated / n_early_codes,
  #       prob_at_least_one_new = at_least_one_new_code / n_early_codes
  #     )
  # }


  list(
    obs_first_hosp = obs_first_hosp#,
    #propagated_codes = propagated_codes
  )
}
