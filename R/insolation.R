# parameters and whole function idea taken from palinsol::Insol
#' Insolation
#'
# #' @param solution Result of prep_solution(). A dataframe with columns ecc, varpi, epl.
#' @param eccentricity The eccentricity. \eqn{e}
#' @param obliquity The axial tilt. \eqn{\epsilon}
#' @param lpx The longitude of perihelion from the moving equinox \eqn{bar omega}{\bar{\omega}}.
#' @param longitude True solar longitude (rad). Defaults to June solstice pi/2.
#' @param latitude Latitude on Earth (rad). Defaults to 65 degrees North.
#' @param S0 Total solar irradiance (W m^-2). Defaults to 1360.7 W/m^2.
#' @param H Sun hour angle (rad). Daily mean by default (NULL).
insolation <- function(#solution, # result of prep_solution
                       eccentricity, obliquity, lpx,
                       longitude = pi / 2,
                       latitude = 65 * pi / 180,
                       S0 = 1361, # different from palinsol 1365!
                       # based on Lunt et al., 2017 DeepMIP => Matthes et al., 2016
                       # "although the early Eocene (51 Ma) solar constant was ∼ 0.43 % less than this (Gough, 1981), i.e. ∼ 1355 W m−2"
                       H = NULL) {
  # now explicitly takes arguments so it's less confusing
  # make variables locally available
  # this needs: varpi, ecc, epl
  # note we have epl in stead of eps (EPsiLon vs EPSilon)
  # NOTE: varpi is actually pibar reframed so this is shit!
  ## for (i in names(solution)) {
  ##   assign(i, solution[[i]])
  ## }
  nu <- longitude - lpx # true anomaly
  rho <- (1 - eccentricity^2) / (1 + eccentricity * cos(nu))
  sindelta <- sin(obliquity) * sin(longitude)
  cosdelta <- sqrt(1 - sindelta^2)
  sinlatsindelta <- sin(latitude) * sindelta
  coslatcosdelta <- cos(latitude) * cosdelta

  if (is.null(H)) {
    cosH0 <- pmin(pmax(-1, -sinlatsindelta / coslatcosdelta), 1)
    sinH0 <- sqrt(1 - cosH0^2)
    H0 <- acos(cosH0)
    insol <- S0 / (pi * rho^2) * (H0 * sinlatsindelta + coslatcosdelta * sinH0)
  } else {
    insol <- max(0, S0 / (rho^2) * (sinlatsindelta + coslatcosdelta * cos(H)))
  }
  return(insol)
}
