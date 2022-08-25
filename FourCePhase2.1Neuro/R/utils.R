first_3 <- function(x) {
  # retains first 3 characters of the ICD code
  substr(x, 1, 3) %>% unique()
}

right_join0 <- function(x, y, fill = 0L, ...){
  z <- right_join(x, y, ...)
  tmp <- setdiff(names(z), names(y))
  tidyr::replace_na(z, setNames(as.list(rep(fill, length(tmp))), tmp))
}

concat <- function(x, y){
  paste0(x, ' (', round(y, 3)*100, '%)')
}

concat_median <- function(med, mi, ma){
  paste0(med, ' [', mi, ', ', ma, ']')
}

concat_mean <- function(mea, s, acc = 0){
  paste0(round(mea, acc), ' (', round(s, acc), ')')
}

nas_to_zeros <- function(x){
  x[is.na(x)] <- 0
  x
}

subtract_days <- function(x0, x1){
  lubridate::interval(x0, x1) %/% lubridate::days(1)
}

tidywhere <- function (fn){
  # https://github.com/r-lib/tidyselect/issues/201#issuecomment-650547846
  # where is not exported in tidyselect
  predicate <- purrr::as_mapper(fn)
  function(x, ...) {
    out <- predicate(x, ...)
    out
  }
}


