# Description: Grabs the raw data from Synapse and stores it in the data-raw folder.
#  We are using the genieBPC function in the hopes that it will support
#  PAT access in the future.
# Author: Alex Paynter

# Before running this file you will need to login to synapse with
#   genieBPC::set_synapse_credentials()
# This login will expire every 24h, and should be replaced with a PAT
#   access ASAP.

library(synapser)
library(here)
library(genieBPC)
library(fs)

dir.create(here("data-raw"), showWarnings = F)
dir.create(here("data"), showWarnings = F)

# Rdata version:
data_list <- genieBPC::pull_data_synapse(cohort = "BrCa",
                                         version = "v1.2-consortium")
# CSV versions just in case:
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

# Save the Rdata version:
readr::write_rds(x = data_list,
                 file = here('data-raw', 'data_list.rds'))

