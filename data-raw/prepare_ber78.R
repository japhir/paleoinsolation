library(readr)
library(tibble)

# from 0 to -100 Myr (in kyr)
times_kyr <- seq(0, -1e5, -1)
# get the Ber78 solution using palinsol
# note time in years is after 1950
ber_xss <- Map(\(x) palinsol::ber78(x), times_kyr * 1e3)
ber <- tibble(Reduce(rbind, ber_xss))
ber <- ber |>
  mutate(time = times_kyr, .before = eps)

ber |>
  write_csv("dat/Ber78_from_palinsol.csv ")
