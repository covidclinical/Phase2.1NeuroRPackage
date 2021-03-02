library(dplyr)
library(forcats)

neuro_icds_10 <-
  readxl::read_excel('data-raw/2020-09-10_neuro-icd10_CNSvPNS.xlsx') %>%
  rename('icd' = `ICD-10`,
         'pns_cns' = `Nervous system Involvement (1=central, 2=peripheral)`) %>%
  mutate(pns_cns = as.factor(pns_cns) %>% fct_recode(Central = '1',
                                                     Peripheral = '2'))

neuro_icds_9 <- read.csv('data-raw/icd9_tab_CNSvPNS.csv') %>%
  rename('Neurological Disease Category' = 'Neurological.Disease.Category',
         'pns_cns' = `Nervous.system.Involvement..1.central..2.peripheral.`,
         'icd_description' = `icd9_desc`) %>%
  mutate(pns_cns = as.factor(pns_cns) %>% fct_recode(Central = '1',
                                                     Peripheral = '2'),
         concept_type = "DIAG-ICD9")

site_params <- googlesheets4::read_sheet(
  'https://docs.google.com/spreadsheets/d/1epcYNd_0jCUMktOHf8mz5v651zy1JALD6PgzobrGWDY/edit?usp=sharing'
)

usethis::use_data(neuro_icds_9, neuro_icds_10, site_params, overwrite = TRUE)
