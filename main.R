# Description:  Top level workflow for the project.  Eventually this should
#  be automated, clean.  For now it's a list of comments and source files.

library(purrr)
library(here)
library(fs)

# Important!
# Manually add the files in https://app.box.com/folder/200927784826 to
# data/msk_box_derived (╭ರ_⊙)

# Load all helper functions
purrr::walk(.x = fs::dir_ls('R'), .f = source)

source(here('analysis', 'script', 'get_raw_data.R'))
source(here('analysis', 'script', 'filter_data_for_cohort.R'))
source(here('analysis', 'script', 'process_data.R'))
rmarkdown::render(
  input = here('analysis', 'report', 'brca_regimens.Rmd'),
  output_file = '01-bpc-brca-regimen-overview.html',
  output_dir = here('output')
)


# No longer needed the following, using inputs from MSK bioinformatics team:
# source(here('analysis', 'script', 'merge_gene_panels.R'))
source(here('analysis', 'script', 'process_msk_gene.R'))
source(here('analysis', 'script', 'create_surv_dat_dmet.R'))

# Lots of stuff here...

# source(here('analysis', 'scripts', 'upload_outputs_synapse.R'))
