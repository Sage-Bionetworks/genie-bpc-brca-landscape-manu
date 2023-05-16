output_brca_synid <- "syn51359343" #2023-04-17-BrCa-landscape-paper-outputs

library(synapser)

synLogin()
synapser::File(here("analysis", "reports", "brca_regimens.html"),
               parent = output_brca_synid) %>%
  synStore()


synpser::File(here("data-raw", "drug_map.csv"),
              parent = output_brca_synid) %>%
  synStore()