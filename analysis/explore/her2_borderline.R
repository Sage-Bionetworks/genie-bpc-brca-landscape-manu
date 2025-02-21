library(purrr); library(here); library(fs)
purrr::walk(.x = fs::dir_ls('R'), .f = source)

ca_ind <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_ca_ind.rds')
)

her2_border <- ca_ind %>% 
  filter(ca_bca_her_summ %in% 'Borderline/equivocal/indeterminant')
her2_border %<>%
  select(record_id, ca_seq, contains('her'), bca_subtype)

path <- readr::read_rds(
  here('data', 'clin_data_cohort', 'dft_path.rds')
)

path_her2_border <- path %>%
  filter(record_id %in% her2_border$record_id)

path_her2_border_long <- path_her2_border %>% 
  select(record_id, path_proc_int,
         contains('path_her')) %>%
  rename(path_her2ihc_1 = path_herihc_1,
         path_her2ish_1 = path_herish_1) %>%
  pivot_longer(
    cols = -c(record_id, path_proc_int),
    names_to = 'measure',
    values_to = 'value'
  ) %>%
  mutate(value = case_when(
    str_trim(value) %in% "Test not done" ~ NA_character_,
    value %in% 'Negative (1+) IHC only' ~ "Negative",
    T ~  value
  )) %>%
  filter(!is.na(value)) %>%
  separate_wider_delim(
    cols = measure, delim = "_",
    names = c('trash', 'lab', 'arbitrary_number'),
    too_few = "error", too_many = "error"
  )

her2_lev <- c(
  "Negative", # Negative 1+ nonsense is gone.
  "Equivocal", 'Positive'
)

path_her2_border_long %<>%
  # filter(lab %in% 'her2ihc') %>%
  mutate(value = factor(value, levels = her2_lev, ordered = T)) %>%
  group_by(record_id, lab) %>% 
  arrange(path_proc_int) %>%
  summarize(
    max = max(value),
    first = first(value),
    .groups = 'drop'
  ) 

path_her2_border_merge <- path_her2_border_long %>%
  pivot_longer(
    cols = c(max, first),
    names_to = 'metric',
    values_to = 'result'
  ) %>%
  mutate(
    new_col = paste0(lab, '_', metric),
  ) %>%
  select(-c(lab, metric)) %>%
  pivot_wider(names_from = new_col, values_from = result) %>%
  rename_at(.vars = vars(-record_id), .funs = \(x) paste0('path_', x))

her2_border %<>% 
  left_join(
    ., path_her2_border_merge, by = c('record_id')
  ) 

readr::write_csv(
  her2_border,
  here('analysis', 'explore', 'her2_borderline_with_path.csv'),
  na = ''
)
