module insol
  implicit none
contains

real elemental function insolation(ecc, obl, lpx, long, lat, Sz) result(ins)
  real, intent(in) :: ecc, obl, lpx
!!$  real, dimension(249481), intent(in) :: ecc, obl, lpx
  real, intent(in) :: long, lat, Sz
  real :: pi, lpxrot, nu, rho, sindelta, cosdelta, sinlatsindelta, coslatcosdelta, cosHz, sinHz, Hz

  lpxrot = mod(lpx - pi, 2. * pi)

  nu = long - lpxrot
  rho = (1 - ecc**2) / (1 + ecc * cos(nu))
  sindelta = sin(obl) * sin(long)
  cosdelta = sqrt(1 - sindelta**2)
  sinlatsindelta = sin(lat) * sindelta
  coslatcosdelta = cos(lat) * cosdelta

  ! if is.null(H)
  cosHz = min(max(-1., -sinlatsindelta / coslatcosdelta), 1.)
  sinHz = sqrt(1. - cosHz**2.)
  Hz = acos(cosHz)
  ins = Sz / (pi * rho**2.)

end function insolation

end module insol
