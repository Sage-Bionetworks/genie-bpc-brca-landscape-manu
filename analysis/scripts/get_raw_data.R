# Description: Grabs the raw data from Synapse and stores it in the data-raw folder.
#  We are using the genieBPC function in the hopes that it will support
#  PAT access in the future.
# Author: Alex Paynter

# Before running this file you will need to login to synapse with
#   genieBPC::set_synapse_credentials()
# This login will expire every 24h, and should be replaced with a PAT
#   access ASAP.

library(cli)
library(synapser)
library(here)
library(genieBPC)
library(fs)
library(purrr)
library(here)

dir.create(here("data-raw"), showWarnings = F)
dir.create(here("data"), showWarnings = F)

# Rdata version:
data_list <- genieBPC::pull_data_synapse(cohort = "BrCa",
                                         version = "v1.2-consortium")
# Save the Rdata version:
readr::write_rds(x = data_list,
                 file = here('data-raw', 'data_list.rds'))


# CSV versions.
# These may not match the processing done by MSK exactly (for example I think
#   they filter down to index cancer regimens only).
synid_clin_data <- 'syn39802381'
synid_gen_data <- 'syn32299078'
vec_gen_files <- c('data_CNA.txt', 
                   'data_fusions.txt', 
                   'data_mutations_extended.txt')
dft_dat_list <- synGetChildren(synid_clin_data) %>%
  as.list %>%
  purrr::map_dfr(.x = .,
                 .f = as_tibble)
dft_gen_list <- synGetChildren(synid_gen_data) %>%
  as.list %>%
  purrr::map_dfr(.x = .,
                 .f = as_tibble) %>%
  filter(name %in% vec_gen_files)
if (any(stringr::str_detect(df_clin_children$name, ".csv^"))) {
  warning("Non-CSV files unexpectedly contained in {synid_clin_data}.")
}
syn_store_in_dataraw <- function(sid) {
  synGet(entity = sid, downloadLocation = here("data-raw"))
}
purrr::walk(.x = df_clin_children$id, 
            .f = syn_store_in_dataraw)

genieBPC::pull_data_synapse(cohort = "BrCa",
                            version = "v1.2-consortium",
                            download_location = "data-raw")

# Move the CSV versions out of the subfolder:
purrr::walk(
  .x = here('data-raw', 'BrCa_v1.2', 
     dir(here('data-raw', 'BrCa_v1.2'))),
  .f = file_move,
  new_path = here('data-raw')
)
fs::dir_delete(here('data-raw', 'BrCa_v1.2'))



# Update: Also want to pull a sheet shared by Brooke M. @ MSK:
synLogin()
classified_meds_synid <- "syn51405609"
synGet(entity = classified_meds_synid,
       downloadLocation = here("data"))
