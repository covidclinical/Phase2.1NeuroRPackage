
## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    '.', ':=', 'Abbreviation', 'HTN', 'HTNcx',
    'Neurological Disease Category', 'Total',
    'admission_date', 'age_group', 'alive',
    'both_neuro', 'charlson_from_comorbid',
    'concept_code', 'concept_type', 'days_since_admission',
    'death_date', 'deceased', 'decimal_to_short',
    'delta_hospitalized', 'diagnosis_code', 'elixhauser_score',
    'explain_table', 'icd10_comorbid', 'icd10_map_quan_deyo',
    'icd10_map_quan_elix', 'icd9_comorbid', 'icd9_map_quan_deyo',
    'icd9_map_quan_elix', 'in_hospital', 'last_discharge_date',
    'long_desc', 'map_df', 'max_elix', 'max_readmis', 'max_time',
    'mean_elix', 'mean_readmis', 'mean_time', 'median_elix',
    'median_readmis', 'median_time', 'min_elix', 'min_readmis',
    'min_time', 'n_Both', 'n_Central', 'n_None',
    'n_Peripheral', 'n_Total', 'n_neuro_cond', 'n_no_neuro_cond',
    'n_patients', '', '', '', 'n_readmissions', 'n_stay', 'n_var',
    'names_charlson', 'names_charlson_abbrev', 'names_quan_elix',
    'names_quan_elix_abbrev', 'neuro_icds_10', 'neuro_post',
    'neuro_pt_post', 'neuro_type', 'non_severe', 'outdir', 'patient_num',
    'pns_cns', 'pres',  'prop', 'prop_patients', 'race',
    'replace_na', 'sd_elix', 'sd_readmis', 'sd_time', 'severe', 'severe_date',
    'sex', 'time_death', 'time_severe', 'time_to_death', 'time_to_severe',
    'value', 'van_walraven_from_comorbid', 'van_walraven_score'
  ))
}
