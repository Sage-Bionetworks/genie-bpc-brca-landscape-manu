
str_filter <- function(vec, pattern, negate = F) {
  str_ind <- str_which(vec, pattern, negate = negate)
  vec[str_ind]
}
