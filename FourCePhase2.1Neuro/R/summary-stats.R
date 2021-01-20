concat <- function(x, y){
  paste0(x, ' (', round(y, 3)*100, '%)')
}

concat_median <- function(med, mi, ma){
  paste0(med, ' [', mi, ', ', ma, ']')
}

concat_mean <- function(mea, s, acc = 0){
  paste0(round(mea, acc), ' (', round(s, acc), ')')
}

severity_stats <- function(df, ...) {
  # summary statistics for severity status
  # count values are obfuscated

  df %>%
    select(neuro_post, time_severe = time_to_severe) %>%
    group_by(neuro_post) %>%
    summarise(median_time = median(time_severe, na.rm = TRUE),
              min_time = min(time_severe, na.rm = TRUE),
              max_time = max(time_severe, na.rm = TRUE),
              mean_time = mean(time_severe, na.rm = TRUE),
              sd_time = sd(time_severe, na.rm = TRUE),
              non_severe = sum(is.na(time_severe)),
              Total = n(),
              .groups = 'drop') %>%
    blur_it(c('non_severe', 'Total'), ...) %>%
    mutate(severe = Total - non_severe) %>%
    transmute(
      neuro_post,
      Nonsevere = concat(non_severe, non_severe/Total),
      Severe = concat(severe, severe/Total),
      `Median time to severity onset [Min, Max] (days)` = concat_median(median_time, min_time, max_time),
      `Mean time to severity onset (SD) (days)` = concat_mean(mean_time, sd_time)) %>%
    pivot_longer(-neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value) %>%
    tibble::column_to_rownames('name')
}

survival_stats <- function(df, ...) {
  # summary statistics for survival status
  # count values are obfuscated

  df %>%
    select(neuro_post, time_death = time_to_death) %>%
    group_by(neuro_post) %>%
    summarise(median_time = median(time_death, na.rm = TRUE),
              min_time = min(time_death, na.rm = TRUE),
              max_time = max(time_death, na.rm = TRUE),
              mean_time = mean(time_death, na.rm = TRUE),
              sd_time = sd(time_death, na.rm = TRUE),
              alive = sum(is.na(time_death)),
              Total = n(),
              .groups = 'drop') %>%
    blur_it(c('alive', 'Total'), ...) %>%
    mutate(deceased = Total - alive) %>%
    transmute(
      neuro_post,
      Alive = concat(alive, alive/Total),
      Deceased = concat(deceased, deceased/Total),
      `Median time to death [Min, Max] (days)` = concat_median(median_time, min_time, max_time),
      `Mean time to death (SD) (days)` = concat_mean(mean_time, sd_time)) %>%
    pivot_longer(-neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value) %>%
    tibble::column_to_rownames('name')
}

nstay_stats <- function(df, ...) {
  # summary statistics for length of stay
  # count values are obfuscated

  df %>%
    select(neuro_post, n_stay) %>%
    group_by(neuro_post) %>%
    summarise(median_time = median(n_stay, na.rm = TRUE),
              min_time = min(n_stay, na.rm = TRUE),
              max_time = max(n_stay, na.rm = TRUE),
              mean_time = mean(n_stay, na.rm = TRUE),
              sd_time = sd(n_stay, na.rm = TRUE),
              .groups = 'drop') %>%
    transmute(
      neuro_post,
      `Median length of stay [Min, Max] (days)` = concat_median(median_time, min_time, max_time),
      `Mean length of stay (SD) (days)` = concat_mean(mean_time, sd_time)) %>%
    pivot_longer(-neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value) %>%
    tibble::column_to_rownames('name')
}

readmission_stats <- function(df, ...) {
  # summary statistics for length of stay
  # count values are obfuscated

  df %>%
    select(neuro_post, n_readmissions) %>%
    group_by(neuro_post) %>%
    summarise(median_readmis = median(n_readmissions, na.rm = TRUE),
              min_readmis = min(n_readmissions, na.rm = TRUE),
              max_readmis = max(n_readmissions, na.rm = TRUE),
              mean_readmis = mean(n_readmissions, na.rm = TRUE),
              sd_readmis = sd(n_readmissions, na.rm = TRUE),
              .groups = 'drop') %>%
    transmute(
      neuro_post,
      `Median number of readmissions [Min, Max]` =
        concat_median(median_readmis, min_readmis, max_readmis),
      `Mean number of readmissions (SD)` =
        concat_mean(mean_readmis, sd_readmis)) %>%
    pivot_longer(- neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value) %>%
    tibble::column_to_rownames('name')
}

elix_stats <- function(df, ...) {
  # summary statistics for Elixhauser score
  # count values are obfuscated

  df %>%
    select(neuro_post, elixhauser_score) %>%
    group_by(neuro_post) %>%
    summarise(median_elix = median(elixhauser_score, na.rm = TRUE),
              min_elix = min(elixhauser_score, na.rm = TRUE),
              max_elix = max(elixhauser_score, na.rm = TRUE),
              mean_elix = mean(elixhauser_score, na.rm = TRUE),
              sd_elix = sd(elixhauser_score, na.rm = TRUE),
              .groups = 'drop') %>%
    transmute(
      neuro_post,
      `Median Elixhauser score [Min, Max]` =
        concat_median(median_elix, min_elix, max_elix),
      `Mean Elixhauser score (SD)` =
        concat_mean(mean_elix, sd_elix)) %>%
    pivot_longer(- neuro_post) %>%
    pivot_wider(names_from = neuro_post, values_from = value) %>%
    tibble::column_to_rownames('name')
}

demo_stats <- function(df, var, ...){
  svar <- sym(var)
  df %>%
    group_by(!!svar) %>%
    count(neuro_post, name = 'n_var') %>%
    as.data.frame() %>%
    blur_it('n_var', ...) %>%
    group_by(neuro_post) %>%
    mutate(both_neuro = sum(n_var)) %>%
    ungroup() %>%
    mutate(prop = n_var/both_neuro,
           pres = concat(n_var, n_var/both_neuro)) %>%
    pivot_wider(- c(n_var, both_neuro), names_from = neuro_post,
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
