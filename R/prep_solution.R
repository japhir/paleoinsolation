prep_solution <- function(astronomical_solution = "PT-ZB18a(1,1)") {
  dat <- get_solution(astronomical_solution)
  sol <- stringr::str_extract(astronomical_solution, "^PT-(ZB\\d{2}[a-z])\\(", group = 1)
  if (sol == "ZB18a") {
    sol <- "ZB18a-100"
  }
  ecc <- get_solution(sol)

  # TODO: instead just apply snvec with output = "all"?
  dat <- dplyr::mutate(dat,
                       # interpolate ecc to this time scale
                       ecc = stats::approx(ecc$time, ecc$ecc, xout = .data$time)$y,
                       # inversely calculate varpi
                       # cp = e * sin(omegabar)
                       # cp / e = sin(omegabar)
                       # asin(cp / e) = omegabar
                       # NOTE: this does result in NaN sometimes!!
                       omegabar = asin(cp / ecc))

  select(dat, all_of(c("time", "ecc", "epl", "phi", "omegabar", "cp")))
}
