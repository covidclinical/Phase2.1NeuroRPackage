
count_stats <- function(df, count_var, neg_var, group_var, ...) {
  # summary statistics for survival/severity status
  # count values are obfuscated

  count_var <- sym(count_var)
  neg_var <- sym(neg_var)
  Count_var <- sym(stringr::str_to_title(count_var))
  group_var <- sym(group_var)

  set.seed(446)
  df %>%
    select(group_var, count_var) %>%
    group_by(!!group_var) %>%
    summarise(!!neg_var := sum(!!count_var == 0),
              Total = n(),
              .groups = 'drop') %>%
    blur_it(c(neg_var, 'Total'), ...) %>%
    mutate(!!count_var := Total - !!neg_var) %>%
    mask_it(count_var, ...) %>%
    # mutate_at(vars(all_of(count_var), all_of(neg_var)),
    #           ~ replace(., is.nan(.)|is.na(.), 0)) %>%
    transmute(
      !!group_var,
      !!neg_var :=
        concat(!!neg_var, if_else(Total == 0, 0, !!neg_var/Total)),
      !!Count_var :=
        concat(!!count_var, if_else(Total == 0, 0, !!count_var/Total))) %>%
    pivot_longer(- !!group_var) %>%
    pivot_wider(names_from = !!group_var,
                values_from = value)
}

continuous_stats <- function(df, cont_var, name, group_var,...) {
  # summary statistics for length of stay
  # count values are obfuscated
  med <- sym(paste('Median', name, '[Min, Max]'))
  mea <- sym(paste('Mean', name, '(SD)'))
  cont_var <- sym(cont_var)
  group_var <- sym(group_var)

  set.seed(446)
  df %>%
    select(group_var, cont_var) %>%
    group_by(!!group_var) %>%
    summarise(median_var = round(median(!!cont_var, na.rm = TRUE),1),
              min_var = round(min(!!cont_var, na.rm = TRUE),1),
              max_var = round(max(!!cont_var, na.rm = TRUE),1),
              mean_var = mean(!!cont_var, na.rm = TRUE),
              sd_var = sd(!!cont_var, na.rm = TRUE),
              .groups = 'drop') %>%
    mutate_at(vars(sd_var, median_var, mean_var),
              ~ replace(., is.nan(.)|is.na(.), 0)) %>%
    transmute(
      !!group_var,
      !!med := concat_median(median_var, min_var, max_var),
      !!mea := concat_mean(mean_var, sd_var)) %>%
    pivot_longer(- !!group_var) %>%
    pivot_wider(names_from = !!group_var, values_from = value)
}

demo_stats <- function(var, df, group_var, blur_abs, mask_thres, ...){

  set.seed(446)

  group_var <- sym(group_var)
  svar <- sym(var)
  df %>%
    group_by(!!svar) %>%
    count(!!group_var, name = 'n_var') %>%
    as.data.frame() %>%
    blur_it('n_var', blur_abs, mask_thres) %>%
    group_by(!!group_var) %>%
    mutate(both_neuro = sum(n_var)) %>%
    ungroup() %>%
    mutate(prop = n_var/both_neuro,
           pres = concat(n_var, n_var/both_neuro)) %>%
    pivot_wider(id_cols = - c(both_neuro),
                names_from = !!group_var,
                values_from = c(n_var, prop, pres),
                values_fill = list(n_var = 0, prop = 0, pres = '0 (0%)')) %>%
    mutate(variable = paste(var, !!svar, sep = '.')) %>%
    select(- !!svar) %>%
    select(variable, everything())
}

blur_it <- function(df, vars, blur_abs, mask_thres){
  # Obfuscate count values.
  # If blurring range is +/-3, or blur_abs = 3,
  # the count receive a small addition of a random number from -3 to 3.
  # If a count is less than mask_thres, set that count to 0.

  set.seed(446)

  for (var in vars){
    var <- sym(var)
    blur_vec <- sample(seq(- blur_abs, blur_abs), nrow(df), replace = TRUE)
    df <- df %>%
      mutate(!!var := !!var + blur_vec,
             !!var := ifelse(abs(!!var) < mask_thres, 0, !!var))
  }
  df
}

blur_mask_int_num <- function(int_or_num_to_blur, blur_abs, mask_thres){
  # Obfuscate single numeric or interger value.
  # If blurring range is +/-3, or blur_abs = 3,
  # the count receive a small addition of a random number from -3 to 3.
  # If a count is less than mask_thres, set that count to 0.

    set.seed(446)

    blur_vec <- sample(seq(- blur_abs, blur_abs), 1, replace = TRUE)

    int_or_num_to_blur = int_or_num_to_blur + blur_vec

    int_or_num_to_blur = ifelse(int_or_num_to_blur < mask_thres, 0, int_or_num_to_blur)

  return(int_or_num_to_blur)
}

mask_it <- function(df, var, blur_abs, mask_thres){
  # Obfuscate count values.
  # If a count is less than mask_thres, set that count to 0.

  set.seed(446)
  var <- sym(var)
  df %>%
    mutate(!!var := ifelse(abs(!!var) < mask_thres, 0, !!var))
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
                                        'readmitted',
                                        'covid_discharged')) {

  total_patients <- length(unique(demo_df$patient_num))

    demo_obfus_table <- lapply(
    X = vars_to_obfs,
    FUN = demo_stats,
    df = demo_df,
    blur_abs = blur_abs,
    mask_thres = mask_thres,
    group_var = group_var
  ) %>%
    bind_rows()

    demo_df <- demo_df %>%
      mutate('covid_discharged' = if_else(covid_discharged == 'Discharged', 1, 0))

  suppressWarnings(
  other_obfus_table <-
    bind_rows(
      continuous_stats(demo_df, 'time_to_first_discharge', 'time to first discharge', group_var),
      continuous_stats(demo_df, 'time_to_last_discharge', 'time to last discharge', group_var),
      count_stats(demo_df, 'severe', 'Nonsevere', group_var, blur_abs, mask_thres),
      continuous_stats(demo_df, 'time_to_severe', 'time to severe', group_var),
      count_stats(demo_df, 'deceased', 'Alive', group_var, blur_abs, mask_thres),
      continuous_stats(demo_df, 'time_to_death', 'time to death', group_var),
      count_stats(demo_df, 'covid_discharged', 'not discharged', group_var, blur_abs, mask_thres),
      continuous_stats(demo_df, 'n_readmissions', 'number of readmissions', group_var),
      continuous_stats(demo_df, 'pre_admission_cns', 'pre admission cns', group_var),
      continuous_stats(demo_df, 'pre_admission_pns', 'pre admission pns', group_var),
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
  )
  list(
    demo_table = demo_obfus_table,
    other_obfus_table = other_obfus_table
  )
}

get_comorb_table <- function(neuro_types,
                             demo_df,
                             scores_unique,
                             comorb_names_elix,
                             blur_abs,
                             mask_thres,
                             group_var = 'neuro_post') {

  total_patients <- length(unique(demo_df$patient_num))

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

  return(elix_obfus_table1)

}

get_table1 <- function(
  df, num_pats, comorbidities,
  pat_col = 'patients', ...
)
{
  col_char <- paste('n', pat_col, sep = '_')
  npat_col <- sym(col_char)
  proppat_col = sym(paste('prop', pat_col, sep = '_'))
  comorbidities_map = comorbidities$Abbreviation

  df <- df %>%
    select(all_of(comorbidities_map)) %>%
    colSums() %>%
    data.frame(n_patients = .)

  df <- df %>% tibble::rownames_to_column("Abbreviation")

  df <- df %>%
    blur_it('n_patients', ...) %>%
    mutate(prop_patients = n_patients/num_pats) %>%
    rename(!!npat_col := n_patients,
           !!proppat_col := prop_patients) %>%
    right_join(comorbidities, ., by = "Abbreviation")

  df
}

list_table1 <- function(x, df, num_pats, comorb_names, group_var, ...){
  group_var <- sym(group_var)
  get_table1(
    df %>% filter(!!group_var == x),
    num_pats,
    comorbidities = comorb_names,
    pat_col = x, ...)
}

count_sequences_hospitalisation <- function(df, ...) {

  date_df = data.frame(calendar_date = seq(min(df$calendar_date), max(df$calendar_date), by = "1 day"))

  seq_hospitalisation_df <- data.frame(days_since_admission = seq(min(df$days_since_admission),
                                                                  max(df$days_since_admission))
  ) %>%
    cbind(., date_df) %>%
    left_join(., df, by = c("days_since_admission", "calendar_date")) %>%
    arrange(days_since_admission) %>%
    replace_na(list(in_hospital = 0)) %>%
    #fill(patient_num, .direction = "downup") %>%
    fill(siteid, .direction = "downup") %>%
    fill(severe, .direction = "down") %>%
    fill(deceased, .direction = "down")

  return(seq_hospitalisation_df)
}


calc_time_period <- function(df) {

  start_date = as.Date('2020-01-01')
  end_date = as.Date('2022-12-31')

  timeline = seq(start_date, end_date, by = "3 months")

 df_date =  tibble( date = seq( start_date, end_date, by = "1 day" ) ) %>%
   mutate(time_period = case_when(
     date >= timeline[1] & date < timeline[2] ~ timeline[1],
     date >= timeline[2] & date < timeline[3] ~ timeline[2],
     date >= timeline[3] & date < timeline[4] ~ timeline[3],
     date >= timeline[4] & date < timeline[5] ~ timeline[4],
     date >= timeline[5] & date < timeline[6] ~ timeline[5],
     date >= timeline[6] & date < timeline[7] ~ timeline[6],
     date >= timeline[7] & date < timeline[8] ~ timeline[7],
     date >= timeline[8] & date < timeline[9] ~ timeline[8],
     date >= timeline[9] & date < timeline[10] ~ timeline[9],
     date >= timeline[10] & date < timeline[11] ~ timeline[10],
     date >= timeline[11] & date < timeline[12] ~ timeline[11],
     date >= timeline[12] ~ timeline[12]
   )) %>%
   rename('admission_date' = date)

 df <- df %>%
   left_join(., df_date, by = 'admission_date')

    return(df)
}

