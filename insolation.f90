module insol
  use kind
  implicit none
contains

real(dp) elemental function insolation(ecc, obl, omegabar, long, lat, Sz) result(ins)

  real(dp), intent(in) :: ecc, obl, omegabar
  real(dp), intent(in) :: long, lat, Sz
  real(dp) :: pi, nu, rho, sindelta, cosdelta, sinlatsindelta, coslatcosdelta, cosHz, sinHz, Hz

  pi = 3.1415926535897932

  nu = long - omegabar ! true anomaly
  rho = (1._dp - ecc**2._dp) / (1._dp + ecc * cos(nu))
  sindelta = sin(obl) * sin(long)
  cosdelta = sqrt(1._dp - sindelta**2._dp)
  sinlatsindelta = sin(lat) * sindelta
  coslatcosdelta = cos(lat) * cosdelta

  ! if is.null(H) ! for now just default
  ! palinsol::Insol uses pmin, pmax here
  ! but because we're using elemental, that shouldn't matter?
  cosHz = min(max(-1._dp, -sinlatsindelta / coslatcosdelta), 1._dp)
  sinHz = sqrt(1._dp - cosHz**2._dp)
  Hz = acos(cosHz)
  ins = Sz / (pi * rho**2._dp) * (Hz * sinlatsindelta + coslatcosdelta * sinHz)

end function insolation

end module insol
