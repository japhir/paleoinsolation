!> Calculate Insolation.
module insol
  use kind, only : dp
  implicit none
  public

contains

!> Calculate Insolation at the top of the Atmosphere
real(dp) elemental function insolation(ecc, obl, lpx, long, lat, Sz) result(ins)

  real(dp), intent(in) :: ecc, obl, lpx
  real(dp), intent(in), optional :: long, lat, Sz
  ! local place-holders for default values
  real(dp) :: long_, lat_, Sz_
  real(dp) :: pi, nu, rho, sindelta, cosdelta, sinlatsindelta, coslatcosdelta, cosHz, sinHz, Hz

  pi = 3.1415926535897932
  ! default values
  long_ = pi/2._dp
  lat_ = 65*pi/180._dp
  Sz_ = 1360.7_dp ! PMIP4 transient simulations default Otto-Bliesner et al., 2017, Menviel et al., 2019
  ! overwrite if passed
  if (present(long)) long_ = long
  if (present(lat)) lat_ = lat
  if (present(Sz)) Sz_ = Sz

  nu = long_ - lpx ! true anomaly
  rho = (1._dp - ecc**2._dp) / (1._dp + ecc * cos(nu))
  sindelta = sin(obl) * sin(long_)
  cosdelta = sqrt(1._dp - sindelta**2._dp)
  sinlatsindelta = sin(lat_) * sindelta
  coslatcosdelta = cos(lat_) * cosdelta

  ! if is.null(H) ! for now just default
  ! palinsol::Insol uses pmin, pmax here
  ! but because we're using elemental, that shouldn't matter?
  cosHz = min(max(-1._dp, -sinlatsindelta / coslatcosdelta), 1._dp)
  sinHz = sqrt(1._dp - cosHz**2._dp)
  Hz = acos(cosHz)
  ins = Sz_ / (pi * rho**2._dp) * (Hz * sinlatsindelta + coslatcosdelta * sinHz)

end function insolation

end module insol
