program paleoinsolation
  use data
  use insol
  use kind
  implicit none

  real(dp), dimension(249481) :: time, ecc, obl, prec, climprec, lpx
  real(dp) :: long, lat, Sz

  real(dp) :: pi
  real(dp) :: OMT, R2D

  real(dp), dimension(249481) :: sixtyfive ! 65°N summer insolation

  pi = 3.1415926535897932
  OMT = 75.594_dp
  R2D = 180._dp / pi

  call readdata(time, ecc, obl, prec, lpx, climprec)
  lpx = modulo(lpx - pi, 2._dp * pi)

  ! call insolation with default 65°N summer insolation
  long = pi / 2._dp
  lat = 65._dp / R2D !pi / 180._dp
  Sz = 1361._dp

  sixtyfive = insolation(ecc, obl, lpx, long, lat, Sz)
  call writedata(time,ecc,obl,prec,lpx,climprec,sixtyfive)

end program paleoinsolation
