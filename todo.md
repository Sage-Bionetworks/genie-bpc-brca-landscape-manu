Catastrophic:

- The surv_fit_dmet_wrap() function was using dft_dmet_os, which will be detrimental for the subset analyses.

For Aug 22:
- Slides explaining the simulation.
- Finish off the lasso run eval functions.
- Check the absolute bias calculation on lasso - seems too good to be true.

Todo:

- Change all files to load in the data in data/cohort_clin_data
- Create a makefile to automate the running of this code.
- Switch brca_regimens.rmd and process_drug_data.R away from using data lists,
   and use the individual datasets instead.
- It's not a rare event at all (>50 cases) for participants to have a report 
  date which is after their death or loss to followup.  Why is this?  Does
  this indeed challenge our assumptions?
- Update all reports to read in the distinctly filtered survival datasets.
- Fix the 5fcv method to return a coef_est with 'estimate' instead of 'log_hr' (just the names are off)

Curiousities:
- Are GENIE-MSK-P-0017706-T02-IM6 and GENIE-MSK-P-0017706-T01-IM6 totally identical?