library(dplyr)

neuro_icds_10 <-
  readxl::read_excel('data-raw/2020-09-10_neuro-icd10_CNSvPNS.xlsx') %>%
  rename('icd' = `ICD-10`,
         'pns_cns' = `Nervous system Involvement (1=central, 2=peripheral)`) %>%
  mutate(pns_cns = as.factor(pns_cns) %>% fct_recode(Central = '1',
                                                     Peripheral = '2'))

usethis::use_data(loinc, neuro_icds_10, overwrite = TRUE)
