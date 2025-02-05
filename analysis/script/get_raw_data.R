# Description: Grabs the raw data from Synapse and stores it in the data-raw folder.
# Author: Alex Paynter

library(purrr); library(here); library(fs);
# Load all helper functions
purrr::walk(.x = fs::dir_ls('R'), .f = source)

synLogin()



# CSV versions - no longer using {genieBPC} for various reasons.  
synid_clin_data <- 'syn39802381'
synid_gen_data <- 'syn32299078'
vec_gen_files <- c(
  'data_CNA.txt', 
  'data_fusions.txt', 
  'data_mutations_extended.txt',
  'data_gene_matrix.txt'
  #'data_gene_panel_DFCI-ONCOPANEL-1.txt'
  # Also includes all gene panels using the logic below.
)
dft_clin_list <- synGetChildren(synid_clin_data) %>%
  as.list %>%
  purrr::map_dfr(.x = .,
                 .f = as_tibble)

dft_gen_list <- synGetChildren(synid_gen_data) %>%
  as.list %>%
  purrr::map_dfr(.x = .,
                 .f = as_tibble) %>%
  filter(name %in% vec_gen_files | str_detect(name, "^data_gene_panel_.*"))

syn_store_help <- function(sid, loc = here('data-raw')) {
  synGet(
    entity = sid, 
    downloadLocation = loc,
    ifcollision = 'overwrite.local' # replaces local copy if it exists.
  )
}
purrr::walk(
  .x = dft_clin_list$id, 
  .f = syn_store_help
)

purrr::walk(
  .x = dft_gen_list$id, 
  .f = (function(x) {
    syn_store_help(x, loc = here('data-raw', 'genomic'))
  })
)

