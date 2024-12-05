prep_solution <- function(orbital_solution = "full-ZB18a", tend = -10e3, ed = 1, td = 1) {
  R2D <- 180.0 / pi
  OMT <- 75.5940
  refsln <- snvec(tend, ed, td, orbital_solution, output = "all") |>
    mutate(time = time * 1e3) |> # time in years so it's compatible with older stuff
    # is it omegabar or varpi* or varpi???
    # ber78 says 'varpi' true solar longitude of perihelion
    # so should just be lph ?
    # omegabar = varpi + | phi |
    # so varpi = lphui + OMT => rewrapped
    mutate(varpi = (lphi + OMT) / R2D) |>
    mutate(omegabar = varpi - phi) |> # what snvec does
    ## mutate(omegabar = varpi + abs(phi)) |> # what it says on the tin
    # crucifix calls omegabar or pibar varpi
    # maybe we need this line from https://github.com/mcrucifix/palinsol/blob/master/data-raw/LA04.R
    ## la04past['varpi'] <- (la04past['varpi'] - pi ) %% (2*pi)
    mutate(omegabar = (omegabar - pi) %% (2 * pi)) |> # what Crucifix does in his prep
    rename(ecc = eei)
  # calculate 65°N summer solstice insolation
  refsln <- refsln |> mutate(insol = insolation(refsln$ecc, refsln$epl, refsln$omegabar,
                                                longitude = pi / 2, latitude = 65 * pi / 180, S0 = 1361))
  refsln |> select(all_of(c("time",
                            "sx", "sy", "sz",
                            "nnx", "nny", "nnz",
                            "ecc", "epl", "phi", "cp", "varpi", "omegabar", "insol")))
}
