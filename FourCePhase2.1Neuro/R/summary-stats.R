concat <- function(x, y){
  paste0(x, ' (', round(y, 3)*100, '%)')
}

concat_median <- function(med, mi, ma){
  paste0(med, ' [', mi, ', ', ma, ']')
}

concat_mean <- function(mea, s, acc = 0){
  paste0(round(mea, acc), ' (', round(s, acc), ')')
}

count_stats <- function(df, count_var, neg_var, group_var, ...) {
  # summary statistics for survival/severity status
  # count values are obfuscated

  count_var <- sym(count_var)
  neg_var <- sym(neg_var)
  Count_var <- sym(stringr::str_to_title(count_var))
  group_var <- sym(group_var)

  df %>%
    select(group_var, count_var) %>%
    group_by(!!group_var) %>%
    summarise(!!neg_var := sum(!!count_var == 0),
              Total = n(),
              .groups = 'drop') %>%
    blur_it(c(neg_var, 'Total'), ...) %>%
    mutate(!!count_var := Total - !!neg_var) %>%
    transmute(
      !!group_var,
      !!neg_var := concat(!!neg_var, !!neg_var/Total),
      !!Count_var := concat(!!count_var, !!count_var/Total)) %>%
    pivot_longer(- !!group_var) %>%
    pivot_wider(names_from = !!group_var, values_from = value)
}

continuous_stats <- function(df, cont_var, name, group_var,...) {
  # summary statistics for length of stay
  # count values are obfuscated
  med <- sym(paste('Median', name, '[Min, Max]'))
  mea <- sym(paste('Mean', name, '(SD)'))
  cont_var <- sym(cont_var)
  group_var <- sym(group_var)

  df %>%
    select(group_var, cont_var) %>%
    group_by(!!group_var) %>%
    summarise(median_var = median(!!cont_var, na.rm = TRUE),
              min_var = min(!!cont_var, na.rm = TRUE),
              max_var = max(!!cont_var, na.rm = TRUE),
              mean_var = mean(!!cont_var, na.rm = TRUE),
              sd_var = sd(!!cont_var, na.rm = TRUE),
              .groups = 'drop') %>%
    transmute(
      !!group_var,
      !!med := concat_median(median_var, min_var, max_var),
      !!mea := concat_mean(mean_var, sd_var)) %>%
    pivot_longer(- !!group_var) %>%
    pivot_wider(names_from = !!group_var, values_from = value)
}

demo_stats <- function(var, df, group_var, ...){
  group_var <- sym(group_var)
  svar <- sym(var)
  df %>%
    group_by(!!svar) %>%
    count(!!group_var, name = 'n_var') %>%
    as.data.frame() %>%
    blur_it('n_var', ...) %>%
    group_by(!!group_var) %>%
    mutate(both_neuro = sum(n_var)) %>%
    ungroup() %>%
    mutate(prop = n_var/both_neuro,
           pres = concat(n_var, n_var/both_neuro)) %>%
    pivot_wider(- c(n_var, both_neuro), names_from = !!group_var,
                values_from = c(n_var, prop, pres)) %>%
    mutate(variable = paste(var, !!svar, sep = '.')) %>%
    replace_na(list(pres_no_neuro_cond = '0 (0%)',
                    pres_neuro_cond = '0 (0%)')) %>%
    select(- !!svar)
}

blur_it <- function(df, vars, blur_abs, mask_thres){
  # Obfuscate count values.
  # If blurring range is +/-3, or blur_abs = 3,
  # the count receive a small addition of a random number from -3 to 3.
  # If a count is less than mask_thres, set that count to 0.

  for (var in vars){
    var <- sym(var)
    blur_vec <- sample(seq(- blur_abs, blur_abs), nrow(df), replace = TRUE)
    df <- df %>%
      mutate(!!var := !!var + blur_vec,
             !!var := ifelse(abs(!!var) < mask_thres, 0, !!var))
  }
  df
}

get_tables <- function(neuro_types,
                       demo_df,
                       scores_unique,
                       comorb_names_elix,
                       blur_abs,
                       mask_thres,
                       group_var = 'neuro_post',
                       vars_to_obfs = c('sex',
                                        'age_group',
                                        'race',
                                        'Severity',
                                        'Survival',
                                        'readmitted')) {

  total_patients <- length(unique(demo_df$patient_num))
  demo_obfus_table <- lapply(
    vars_to_obfs,
    demo_stats,
    df = demo_df,
    blur_abs = blur_abs,
    mask_thres = mask_thres,
    group_var = group_var
  ) %>%
    bind_rows()

  other_obfus_table <-
    bind_rows(
      continuous_stats(demo_df, 'n_stay', 'length of stay', group_var),
      count_stats(demo_df, 'severe', 'Nonsevere', group_var, blur_abs, mask_thres),
      continuous_stats(demo_df, 'time_to_severe', 'time to severe', group_var),
      count_stats(demo_df, 'deceased', 'Alive', group_var, blur_abs, mask_thres),
      continuous_stats(demo_df, 'time_to_death', 'time to death', group_var),
      continuous_stats(demo_df, 'n_readmissions', 'number of readmissions', group_var),
      continuous_stats(
        demo_df,
        'time_to_first_readmission',
        'time to first readmission',
        group_var
      ),
      continuous_stats(
        scores_unique,
        'elixhauser_score',
        'Elixhauser score',
        group_var
      )
    )

  elix_obfus_table1 <-
    Reduce(
      function(...)
        left_join(..., by = c("Comorbidity", "Abbreviation")),
      lapply(
        neuro_types,
        list_table1,
        df = scores_unique,
        num_pats = total_patients,
        comorb_names = comorb_names_elix,
        group_var = group_var,
        blur_abs = blur_abs,
        mask_thres = mask_thres
      )
    ) %>%
    mutate(
      n_Total = rowSums(select(., starts_with('n_'))),
      prop_Total = n_Total / total_patients
    ) %>%
    arrange(desc(n_Total))

  list(
    demo_table = demo_obfus_table,
    other_obfus_table = other_obfus_table,
    elix_obfus_table1 = elix_obfus_table1
  )
}
