
# GENIE BPC Breast Cancer Landscape Manuscript (Sage part)

# Overview

An initial descriptive analysis of AACR's [GENIE BPC](https://www.aacr.org/professionals/research/aacr-project-genie/bpc/) breast cancer cohort.  This repo contains only the analysis portions completed by Sage Bionetworks (originally Alex Paynter), destined for merger with an MSK repository later on.

# Installation and Setup

To clone this repository, run the following in the command line from a machine with git installed:

```
git clone https://github.com/Sage-Bionetworks/genie-bpc-brca-landscape-manu/tree/manuscript-minimum
```

## Reproducibility tools

This repository:
- Was tested and run on R version 4.4.2.
- Uses R projects.  When running any codes, please open the `.RProj` file first.  
- Does **not** use `renv` to manage package environments.
- Does **not** use `docker` or other containerization to manage deployment.

The code may work without appreciation of these tools, but no guarantees.

## Requirements

To run the code in this respository you will need:

- A Synapse account which has taken the quiz granting data download rights for GENIE (**link?**)
- The [synapser](https://r-docs.synapse.org/articles/synapser.html) R package, which will also require a python installation (follow instructions at that link).
	- *Note:*  This is only used to acquire the data.  It is technically possible to do this manually if needed.

# Code structure

The top-level workflow of the project is in `main.R`.  This calls the other analysis scripts in the correct sequence to reproduce my workflow.  Other top level folders include:

- `/analysis` - Scripts (`analysis/scripts`), quarto/rmarkdown files (`analysis/reports`) and any other analysis code excluding function definitions.
- `/data-raw` - Raw data, or as raw as it comes in release.
- `/data` - Processed data, saved at various stages in the analysis.
- `/output` - Figures, rendered reports, tables, etc.
- `/R` - Function definitions.  These are sometimes written with {roxygen}-style documentation.


# Data

We use GENIE BPC Breast release version 1.2, which is only available to GENIE consortium members.  It should closely match the forthcoming Breast 2.0-public release, expected [here](https://www.synapse.org/Synapse:syn27056172/wiki/616631) sometime in the year 2025.  

The structure, processing and flow of data is described in detail in the PDF data guide, which will come out with the release.


# Acknowledgments/References

We wish to thank the following groups for their upstream contributions to the data:

- [AACR Project GENIE team](https://www.aacr.org/professionals/research/aacr-project-genie/about-us/)
- Sage bionetworks GENIE team - processing and releases.
- [MSKCC biostatistics team](https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/project-genie-bpc-genomics-evidence-neoplasia-information-exchange-biopharma-collaborative)
- The patients and institutions who contributed data to the GENIE and GENIE BPC registries.

# License

The license for this material is [GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/).  If you want to convince me that allowing closed-source decendents of this project is a better policy, I'm interested in that topic.