
map_char_elix_codes <- function(df, comorb_names, t1, t2, map_type, truncate = TRUE) {

  df <- df %>%
    filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"),
           days_since_admission >= t1 & days_since_admission <= t2)

  # Create separate df frames for ICD9 and 10 Codes
  # icd package does not support simultaneous processing of both ICD code types
  # we will recombine after the initial processing
  icd10 <- df %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD10") %>%
    distinct()

  icd9 <- df %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD9") %>%
    distinct()

  # select map_type for analysis
  #if (map_type == "charlson") {
  #  icd10_comorb_map = icd::icd10_map_charlson
  #  icd9_comorb_map = icd::icd9_map_charlson
  #}

  # quan prefixes revised charlson & elixhauser mapping
  if (map_type == "charlson") {
    icd10_comorb_map = icd::icd10_map_quan_deyo
    icd9_comorb_map = icd::icd9_map_quan_deyo
  }
  if (map_type == "elixhauser") {
    icd10_comorb_map = icd::icd10_map_quan_elix
    icd9_comorb_map = icd::icd9_map_quan_elix
  }

  ## Because the 4CE has truncated ICD codes, we will also truncate the icd package index maps
  # Function to select first 3 characters of the ICD Code in all lists of the index map
  if (truncate == TRUE) {
    icd10_comorb_map <- lapply(icd10_comorb_map, first_3)
    icd9_comorb_map <- lapply(icd9_comorb_map, first_3)
  }

  if (truncate == FALSE) {
    # convert icd code to short format (without decimals) to facilitate mapping
    # where diagnosis code is you non-truncated icd column
    icd10 <- icd10 %>%
      mutate(diagnosis_code = decimal_to_short(diagnosis_code)) %>%
      select(-concept_code) %>%
      rename(concept_code = diagnosis_code)

    icd9 <- icd9 %>%
      mutate(diagnosis_code = decimal_to_short(diagnosis_code)) %>%
      select(-concept_code) %>%
      rename(concept_code = diagnosis_code)

  }

  # perform the mapping
  icd10_map <-
    icd::icd10_comorbid(
      icd10,
      map = icd10_comorb_map,
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE,
    )

  icd9_map <-
    icd::icd9_comorbid(
      icd9,
      map = icd9_comorb_map,
      icd_name = "concept_code",
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE
    )

  # If multiple rows due to a patient having both ICD 9 and 10 codes, we will take the max of the column
  # This will allow us to capture the 1s indicating that the comorbidity is present
  # the try wrapper is important in cases there is not an instance of a specific comorbidity in the data - try silences errors
  icd_map <- rbind(icd9_map, icd10_map) %>% # results of 9 and 10 mapping
    group_by(patient_num) %>%
    summarise(
      across(everything(), ~ try(max(.x), silent = TRUE)),
      .groups = 'drop'
    )

  ## Calculate Index Scores
  if (map_type == 'charlson') {
    charlson_score <- charlson_from_comorbid(
      icd_map,
      visit_name = "patient_num",
      scoring_system = "charlson",
      hierarchy = TRUE
    ) %>%
      data.frame(charlson_score = .) %>%
      tibble::rownames_to_column("patient_num")
  }

  # need to check If I've done this correctly - it seems to be that their are different versions
  # of the elixhuaser mapping and that in one HTN is combined. I think this is the version needed
  # to run the van_walraven_from_comorb function
  if (map_type == "elixhauser") {
    # combine hypertension into one category
    icd_map <- icd_map %>%
      mutate(HTN = pmax(HTN, HTNcx, na.rm = TRUE)) %>%
      select(-HTNcx)

    van_walraven_score <- icd::van_walraven_from_comorbid(
      icd_map,
      visit_name = 'patient_num',
      hierarchy = TRUE
    ) %>%
      data.frame(van_walraven_score = .) %>%
      tibble::rownames_to_column("patient_num")
  }

  if(map_type == 'charlson') {
    index_scores <- icd_map %>%
      full_join(charlson_score, by = "patient_num") %>%
      arrange(desc(charlson_score)) %>%
      select(
        patient_num,
        charlson_score,
        everything()
      )
  } else {
    index_scores <- icd_map %>%
      full_join(van_walraven_score, by = "patient_num") %>%
      arrange(desc(van_walraven_score)) %>%
      select(
        patient_num,
        van_walraven_score,
        everything()
      )
  }


  # Identify the specific codes that mapped
  # unlist the charlson mapping lists
  icd10$concept_code <- as.character(icd10$concept_code)
  icd9$concept_code <- as.character(icd9$concept_code)

  icd10_map <-
    purrr::map_df(icd10_comorb_map, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    mutate(concept_code = as.character(concept_code)) %>%
    distinct() %>%
    # merge the mapping dataframe to the patient level ICD codes
    # this will return all comorbidities that mapped to our patient data
    inner_join(icd10, by = "concept_code")

  icd9_map <-
    purrr::map_df(icd9_comorb_map, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    mutate(concept_code = as.character(concept_code)) %>%
    distinct() %>%
    inner_join(icd9, by = "concept_code")

  # explain_codes will add additional information regarding the code name
  # add if statements in order to handle sites that only have ICD 9 or 10 codes but not both
  if (nrow(icd10_map) > 0) {
    icd10_mapped_table <- unique(icd10_map$concept_code) %>%
      icd::explain_table() %>%
      left_join(icd10_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Abbreviation, long_desc) %>%
      distinct()
  }

  if (nrow(icd9_map) > 0) {
    icd9_mapped_table <- unique(icd9_map$concept_code) %>%
      icd::explain_table() %>%
      left_join(icd9_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Abbreviation, long_desc) %>%
      distinct()
  }

  # Bind both ICD 9 and 10 code tables together
  mapped_codes_table <-
    icd10_mapped_table %>%
    try(bind_rows(icd9_mapped_table), silent = TRUE) %>%
    # calculate how many patients had each unique comorbidity/concept_code
    count(concept_code, Abbreviation, long_desc, name = 'n_patients') %>%
    arrange(desc(n_patients))

  map_results <-
    list(
      index_scores = index_scores,
      mapped_codes_table = mapped_codes_table
    )

  map_results

}

# where df takes in the matrix for the initial mapping

get_table1 <- function(
  df, num_pats, comorbidities,
  pat_col = 'patients', ...
)
  {
  npat_col = sym(paste('n', pat_col, sep = '_'))
  proppat_col = sym(paste('prop', pat_col, sep = '_'))
  comorbidities_map = comorbidities$Abbreviation

  df %>%
    select(all_of(comorbidities_map)) %>%
    colSums() %>%
    data.frame(n_patients = .) %>%
    tibble::rownames_to_column("Abbreviation") %>%
    blur_it('n_patients', ...) %>%
    mutate(prop_patients = n_patients/num_pats) %>%
    rename(!!npat_col := n_patients,
           !!proppat_col := prop_patients) %>%
    right_join(comorbidities, ., by = "Abbreviation")
}

# df =scores_unique %>% filter(neuro_post == 'None')
list_table1 <- function(x, df, num_pats, comorb_names, ...){
  get_table1(
    df %>% filter(neuro_post == x),
    num_pats,
    comorbidities = comorb_names,
    pat_col = x, ...)
}

process_tables <- function(index_scores, ...) {
  get_table1(index_scores %>% filter(patient_num %in% neuro_pt_post), ...) %>%
    rename('n_neuro_pats' = n_patients) %>%
    left_join(get_table1(index_scores, ...),
              by = c("Comorbidity", "Abbreviation"))
}

get_charlson_names <- function(){
  data.frame(
    Comorbidity = do.call(rbind, icd::names_charlson),
    Abbreviation = do.call(rbind, icd::names_charlson_abbrev))
}

get_quan_elix_names <- function(){
  data.frame(
    Comorbidity = do.call(rbind, icd::names_quan_elix),
    Abbreviation = do.call(rbind, icd::names_quan_elix_abbrev))
}

first_3 <- function(x) {
  # retain first 3 characters of the ICD code
  substr(x, 1, 3) %>% unique()
}
