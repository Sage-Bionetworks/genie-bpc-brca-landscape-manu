# Description:  Top level workflow for the project.  Eventually this should
#  be automated, clean.  For now it's a list of comments and source files.

library(purrr)
library(here)
library(fs)

# Important!
# Some files must be manually added to build on other's work:
# - Files in https://app.box.com/folder/200927784826 to
# data/msk_box_derived (╭ರ_⊙)

# Load all helper functions
purrr::walk(.x = fs::dir_ls('R'), .f = source)

source(here('analysis', 'script', 'folder_setup.R'))
source(here('analysis', 'script', 'get_raw_data.R'))

#################
# Old Pipeline: #
#################
# source(here('analysis', 'script', 'filter_data_for_cohort.R'))
# source(here('analysis', 'script', 'save_rds_msk_gene.R'))
# source(here('analysis', 'script', 'combine_gene_feat.R'))
# source(here('analysis', 'script', 'prepare_data_for_oncokb_annotate.R'))
# # run annotate_oncokb.sh from the command line.  See comments on enviro vars.
# source(here('analysis', 'script', 'create_gene_panel_dat.R'))
# source(here('analysis', 'script', 'filter_oncogenic_create_features.R')) # some room to trim here.
# source(here('analysis', 'script', 'gene_feat_prep.R'))
# source(here('analysis', 'script', 'clin_feature_prep_dmet.R')) # Added since v1.


#################
# New Pipeline: #
#################
source(here('analysis', 'script', 'filter_data_for_cohort.R'))
source(here('analysis', 'script', 'create_gene_panel_dat.R'))
# new script
source(here('analysis', 'script', 'clin_feature_prep_dmet.R'))


##################
# Survival Part: #
##################
source(here('analysis', 'script', 'surv_prep_dmet_2.R'))
source(here('analysis', 'script', 'surv_fit_dmet_2.R'))
source(here('analysis', 'script', 'surv_process_results_dmet_2.R'))
rmarkdown::render(
  input = here('analysis', 'report', 'bpc-breast-surv-dmet_2.Rmd'),
  output_file = '02-bpc-breast-surv-dmet-v2.pdf',
  output_dir = here('output')
)

source(here('analysis', 'script', 'output_manu_figs.R'))



# A separate thread for Evan's analysis of age and genomics:
source(here('analysis', 'script', 'surv_age_bca_prep_plot.R'))







