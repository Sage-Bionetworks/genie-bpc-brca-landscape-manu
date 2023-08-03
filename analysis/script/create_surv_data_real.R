# Description:  Permute the dmet surv data for simulations.

library(fs)
library(here)
library(purrr)

purrr::walk(.x = fs::dir_ls('R'), .f = source)

dft_surv_dmet <- readr::read_rds(
  file = here('data', 'survival', 'prepared_data',
              'surv_dmet_all.rds')
)

# grab the binary covariates for now:
sim_skel <- dft_surv_dmet %>% 
  select(
    record_id, 
    white, hispanic, 
    matches("_mut$"),
    matches("_cna$"),
    matches("_fus$")
  )

colnames(sim_skel) <- c(
  "sim_rec_id", 
  paste0(
    "x_", 
    str_pad(
      1:(ncol(sim_skel)-1), 
      width = 3,
      pad = 0, 
      side = "left"
    )
  )
)

sim_skel %<>% 
  mutate(
    sim_rec_id = paste0(
      "rec_", 
      str_pad(
        1:n(), 
        width = 4, 
        pad = 0, 
        side = "left"
      )
    )
  ) %>%
  # make sure all the non-record ID columns are doubles:
  mutate(
    across(
      .cols = -sim_rec_id,
      .fns = as.double
    )
  )

test_x <- select(sim_skel, -sim_rec_id)
test_beta <- test_x %>% 
  create_surv_beta_unif(., proportion_null = 0.8)
test_x %<>% as.matrix(.)

test_beta*test_x

-log(runif(nrow(test_x), 0, 1))/exp(t(test_x)%*%test_beta)

# A super simple example to make sure I'm lot losing it:
# The frist covariate increases risk, the rest do nothing.
test_beta = c(0.7, 0, 0, 0)
# Only participant #2 has the fist coviarate positive, so they should be at higher risk.
test_x = matrix(c(0, 1, 1, 1,
                  1, 0, 0, 0,
                  0, 1, 0, 0),
                nrow = 3,
                byrow = T)
test_x %*% test_beta # col vector output
t(test_beta) %*% t(test_x) # row vector output.
purrr::map_dfr(.x = 1:100, .f = (function(do_nothing) {
  vec <- -log(runif(nrow(test_x), 0, 1))/exp(test_x %*% test_beta)
  tibble(y1 = vec[1], y2 = vec[2], y3 = vec[3])
})) %>%
  colMeans(.)
# They tend to die sooner, great.

# https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring
# taking the parametrization h(t|a,p) = pa^pt^(p-1), so that
# H(t) = a^pt^p, and if p < 1 (shape < 1) the hazard is decreasing, p > 1 => increasing.
# See Marco Carone's course notes for this parametrization.
# p = shape, lambda = scale.
gen_times_cox_weib <- function(beta, x, shape, scale) {
  num <- log(runif(n = nrow(x), 0, 1))
  denom <- as.vector((scale^shape)*exp(x %*% beta)) # constant * rel haz
  draws <- (-num/denom)^(1/shape)

  return(draws)
}

# test one:
gen_times_cox_weib(test_beta, test_x, shape = 0.5, scale = 1)

gen_times_cox_weib(test_beta, test_x, shape = 0.7, scale = 0.25)
# test mnay times:
purrr::map_dfr(.x = 1:100, .f = (function(do_nothing) {
  vec <- gen_times_cox_weib(test_beta, test_x, shape = 0.9, scale = 0.25)
  tibble(y1 = vec[1], y2 = vec[2], y3 = vec[3])
})) %>%
  colMeans(.)
# Mean for the undisturbed (no covaritaes affect) vectors should be
gamma(1+1/0.9)/0.25
# So that seems about right.

# Next steps:  Generate data using the cox model weibull wrt betas for survival.  You could use this with beta = 0 for truncation times, or go uniform.  Uniform for censoring makes some sense too.




