library(purrr); library(here); library(fs);
purrr::walk(.x = fs::dir_ls(here('R')), .f = source)

read_wrap <- function(p) {
  read_rds(file = here("data", "clin_data_cohort", p))
}

dft_pt <- read_wrap("dft_pt.rds")
dft_ca_ind <- read_wrap("dft_ca_ind.rds")
dft_reg <- read_wrap("dft_reg.rds")




dft_drug <- dft_reg %>% 
  select(
    record_id, ca_seq, regimen_number,
    matches("_[1-5]$"),
    matches("tt_os_d[1-5]_days$")
  ) %>%
  mutate(
    across(
      .cols = -c(record_id, ca_seq, regimen_number),
      .fn = as.character # for now
    ) 
  ) 

dft_drug %<>%
  pivot_longer(
    cols = -c(record_id, ca_seq, regimen_number)
  ) %>%
  mutate(
    drug_number = readr::parse_number(name),
    # pattern for most variables:
    name = str_replace(name, "_[1-5]$", ""),
    # pattern for the tt_os_d[#]_days variables:
    name = str_replace(name, "^tt_os_d[1-5]_days$", "tt_os_days")
  ) %>%
  pivot_wider(
    names_from = name,
    values_from = value
  )

dft_drug %<>% 
  rename(drug = drugs_drug) %>%
  select(-drugs_drug_oth) # just no usable information there sadly.

dft_drug %<>%
  filter(!is.na(drug))

dft_drug %<>%
  mutate(
    across(
      .cols = matches("_int|_days$"),
      .fns = as.numeric
    )
  )

readr::write_rds(
  x = dft_drug,
  file = here('data', 'clin_data_cohort', 'drug.rds')
)
