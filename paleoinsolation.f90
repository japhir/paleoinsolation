program paleoinsolation
  use data
  use insol

  implicit none
  real, dimension(249481) :: time, ecc, obl, prec, climprec, lpx
  real :: long, lat, Sz
  real :: pi
  real, dimension(249481) :: sixtyfive
  real :: tz

  call readdata(time, ecc, obl, prec, climprec, lpx)

  print *, 'called readdata.'
  print *, 'at t0'
  print *, 'Time = ', time(1)
  print *, 'Eccentricity = ', ecc(1)
  print *, 'Obliquity = ', obl(1)
  print *, 'Precession = ', prec(1)
  print *, 'Climatic Precession = ', climprec(1)
  print *, 'Longitude of Perihelion of the Ascending Node = ', lpx(1)

  print *, 'at t1'
  print *, 'Time = ', time(2)
  print *, 'Eccentricity = ', ecc(2)
  print *, 'Obliquity = ', obl(2)
  print *, 'Precession = ', prec(2)
  print *, 'Climatic Precession = ', climprec(2)
  print *, 'Longitude of Perihelion of the Ascending Node = ', lpx(2)


  print *, 'at tfin line 249481'
  print *, 'Time = ', time(249481)
  print *, 'Eccentricity = ', ecc(249481)
  print *, 'Obliquity = ', obl(249481)
  print *, 'Precession = ', prec(249481)
  print *, 'Climatic Precession = ', climprec(249481)
  print *, 'Longitude of Perihelion of the Ascending Node = ', lpx(249481)

  ! call insolation with default 65°N summer insolation
  long = pi / 2
  lat = 65 * pi / 180
  Sz = 1361

  tz = insolation(ecc(1), obl(1), lpx(1), long, lat, Sz)
  print *, '65°N summer insolation at t0 = ', tz

  sixtyfive = insolation(ecc, obl, lpx, long, lat, Sz)
  print *, '65°N summer insolation at t0 = ', sixtyfive(1)
  print *, '65°N summer insolation at tend = ', sixtyfive(249481)

  ! TODO
  ! call writedata(time,ecc,obl,prec,climprec,lpx,ins)

end program paleoinsolation
