!> The paleoinsolation program makes recent astronomical solutions
!> readily available for use in palaeoclimate simulations, such as the CESM.
!>
!> the Makefile:
!> downloads the astronomical solution ZB18a to ems-plan3.dat
!> available on https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html
!> OR uses the files for the full 300 Myr (default)
!> downloads, patches, and compiles the snvec c-program
!> https://github.com/rezeebe/snvec
!> calculates the PT-solution with Ed = 1 and Td = 1, resulting in dat/PT-ZB18a_1-1.dat
!> compiles this program and runs it to generate ins.dat for 65°N summer insolation
!> we also show how to linearly interpolate our astronomical solutions
!> and how to generate a grid of Solar longitudes, Earth latitudes
program paleoinsolation
  use kind, only : dp
  use data, only : readdata, readbindata, writedata
  use interp, only : locate
  use orb, only : orbpar
  use insol, only : insolation
  use shr_kind_mod, only : SHR_KIND_R8, SHR_KIND_IN
  use shr_orb_mod, only : shr_orb_params, SHR_ORB_UNDEF_INT
  implicit none

  ! some constants
  ! pi = 3.1415926535897932_dp
  real(dp), parameter :: pi = 4.0_dp*datan(1.0_dp)
  real(dp), parameter :: R2D = 180._dp / pi ! radians to degrees

  ! the variables that hold ZB18a(1,1) input
  real(dp), allocatable :: &
       time(:), &
       ecc(:), &
       obl(:), &
       prec(:), &
       lpx(:), &
       climprec(:)
  ! insolation output
  real(dp), allocatable :: &
       interpolate_time(:), & 
       interpolate_ecc(:), & 
       interpolate_obl(:), & 
       interpolate_lpx(:), & 
       interpolate_insolation(:)

  ! a grid of lat/lon to hold insolation
  real(dp), allocatable :: latlons(:,:,:)
  ! the length of time(:) etc.
  integer :: n, io
  ! desired true Solar longitude, Earth's latitude, Solar constant
  ! input parsing
  ! in degrees
  real(dp) :: S0
  ! converted to in radians
  real(dp) :: long, lat

  ! get orbital parameters at a specific calendar year
  real(dp) :: yearBP
  real(dp) :: yearCE
  ! output a single value in radians
  real(dp) :: ecc1, obl1, lpx1
  ! also in degrees
  real(dp) :: obl1_deg, lpx1_deg
  real(dp) :: yAD

  integer(SHR_KIND_IN) :: iyear_AD
  real(SHR_KIND_R8) :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp

  ! desired grid of Solar longitudes and Earth's latitudes to calcualte insolation for
  integer :: i, j
  real(dp), dimension(5) :: longs
  real(dp), dimension(7) :: lats

  ! globally shared defaults for tests
  long = pi / 2._dp
  lat = 65._dp / R2D !pi / 180._dp
  S0 = 1360.7_dp ! the input total insolation

  print *,'using orbpar'
  ! interpolate astronomical solution to a single calendar year
  ! set desired year
  !yearBP =    -800000.0_dp ! i.e. -8 Ma = 8 kyr into the future, should throw
  yearBP = 0._dp ! J2000.0
!!$  yearBP = -(yearCE - 2000._dp) ! J2000.0
  ! convert from year before present to calendar year
  ! the model output has t0 = J2000.0
  ! (but note that model years are 365.25 days long)
  yearCE = -yearBP + 2000.0_dp
  ! or set years directly
!!$  yearCE = 1990
!!$  yearCE = 2100 ! this will throw, because we did not calculate for the future.

  call orbpar(yearCE,ecc1,obl1,lpx1)

  print *,'--------------------------------------------------------------------------------'
!!$!$  print *,'linearly interpolated astronomical solution at ',yearBP,' BP'
  print *,'linearly interpolated astronomical solution at ',yearBP*1.0e-6_dp,' Ma'
  ! convert from radians to degrees
  obl1_deg = obl1 * R2D
  lpx1_deg = lpx1 * R2D
  print *, 'eccentricity: ', ecc1, '[-]'
  print *, 'obliquity:    ', obl1_deg, '°'
  print *, 'longitude of perihelion w/ respect to moving equinox − 180°: ', lpx1_deg, '°'

  print *,'--------------------------------------------------------------------------------'

  print *,'using shr_orb_params'
  ! same but implemented with the ESCOMP/CDEPS API
  iyear_AD = int(-yearBP) + 2000_SHR_KIND_IN
  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
  print *,'linearly interpolated astronomical solution at ',iyear_AD,' CE'
  print *,'using shr_orb_params'

  ! tests
  ! test just passing the parameters for J2000.0
  iyear_AD = SHR_ORB_UNDEF_INT
  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.) ! just specify first parameters, then calculate remainder
  ! check if unreasonable checks work
!!$  obliq = SHR_ORB_UNDEF_REAL
!!$  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
!!$  obliq = -91.0_SHR_KIND_R8 ! too low
!!$  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
!!$  obliq = 91.0_SHR_KIND_R8 ! too high
!!$  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
!!$  eccen = -0.1_SHR_KIND_R8
!!$  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
!!$  eccen = 0.11_SHR_KIND_R8
!!$  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
!!$  mvelp = -0.1_SHR_KIND_R8
!!$  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
!!$  mvelp = 360.1_SHR_KIND_R8
!!$  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.)
  ! these should error:
!!$  call shr_orb_params(-300002001_SHR_KIND_IN,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.) ! past >300 Ma
!!$  call shr_orb_params(      2001_SHR_KIND_IN,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.) ! future
  ! these should warn
  iyear_AD = -65998000_SHR_KIND_IN
  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.) ! 66 Ma, warn ZB18a
  iyear_AD = -82998000_SHR_KIND_IN
  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.) ! 83 Ma, warn ZB18a and ZB20a
  ! this should run without warnings/errors
  iyear_AD = -8000_SHR_KIND_IN
  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.) ! 10 ka, within bounds

  print *,'--------------------------------------------------------------------------------'
  print *, 'calculate 65°N peak summer insolation for last shr_orb_param timestep:'
  print *, 'year: ', iyear_AD
  print *, 'insolation: ', insolation(eccen, obliqr, mvelpp, long, lat, S0)
  print *,'--------------------------------------------------------------------------------'
  ! see if linear interpolation results in weird kinks
  ! to make sure we're OK without fancier interpolation
  !
  ! first datapoint in solution is approx. -0.379 kyr
  ! so that's 1621 ish

  print *, 'interpolating solution to years around a kink'
  n=20

  allocate( &
       interpolate_time(n), &
       interpolate_ecc(n), &
       interpolate_obl(n), &
       interpolate_lpx(n))

  yAD = 1609_dp
  do i=1,n
     interpolate_time(i) = yAD+i
     call orbpar(real(yAD+i, kind = 8), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          interpolate_lpx(i))
  end do
  print *, 'interpolated to times between', yAD, ' and ', yAD+n

  open(newunit=io, file = 'out/interp_1620.dat', status="replace", action="write")
  do i=1,n
     write(io,*) interpolate_time(i), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          interpolate_lpx(i)
  enddo
  close(io)

  deallocate( &
       interpolate_time, &
       interpolate_ecc, &
       interpolate_obl, &
       interpolate_lpx)

  print *, 'wrote linearly interpolated orbital forcing to file'

  print *,'--------------------------------------------------------------------------------'
  print *, 'interpolating solution to years around a kink'
  n=20

  allocate( &
       interpolate_time(n), &
       interpolate_ecc(n), &
       interpolate_obl(n), &
       interpolate_lpx(n))

  yAD = -9820_dp
  do i=1,n
     interpolate_time(i) = yAD+i
     call orbpar(real(yAD+i, kind = 8), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          interpolate_lpx(i))
  end do
  print *, 'interpolated to times between', yAD, ' and ', yAD+n

  open(newunit=io, file = 'out/interp_-9849.dat', status="replace", action="write")
  do i=1,n
     write(io,*) interpolate_time(i), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          interpolate_lpx(i)
  enddo
  close(io)

  deallocate( &
       interpolate_time, &
       interpolate_ecc, &
       interpolate_obl, &
       interpolate_lpx)

  print *, 'wrote linearly interpolated orbital forcing to file'

  print *,'--------------------------------------------------------------------------------'

  print *, 'interpolating solution for past 40 kyr to check'
  n = 421

  allocate(interpolate_time(n), &
       interpolate_ecc(n), &
       interpolate_obl(n), &
       interpolate_lpx(n), &
       interpolate_insolation(n))

  yAD = -40000._dp
  do i=1,n
     interpolate_time(i) = yAD+(i-1)*100
     call orbpar(real(yAD+(i-1)*100, kind = 8), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          interpolate_lpx(i))
  end do
  print *, 'interpolated to times between', yAD, ' and ', yAD+(n-1)*100, ' AD, in ', 100, ' yr steps'
  interpolate_insolation = insolation(interpolate_ecc, &
       interpolate_obl, &
       interpolate_lpx, &
       long, lat, S0)

  open(newunit=io, file = 'out/interp_-40000.dat', status="replace", action="write")
  do i=1,n
     write(io,*) interpolate_time(i), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          modulo(interpolate_lpx(i)-pi,2.0_dp*pi), &
          interpolate_insolation(i)
  enddo
  close(io)

  deallocate( &
       interpolate_time, &
       interpolate_ecc, &
       interpolate_obl, &
       interpolate_lpx, &
       interpolate_insolation)

  print *, 'wrote linearly interpolated orbital forcing to file'

  print *,'--------------------------------------------------------------------------------'

  print *, 'interpolating solution from 250 kyr'
  n = 421

  allocate( &
       interpolate_time(n), &
       interpolate_ecc(n), &
       interpolate_obl(n), &
       interpolate_lpx(n), &
       interpolate_insolation(n))

  yAD = -250000.0_dp
  do i=1,n
     interpolate_time(i) = yAD+(i-1)*100
     call orbpar(real(yAD+(i-1)*100, kind = 8), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          interpolate_lpx(i))
  end do
  print *, 'interpolated to times between', yAD, ' and ', yAD+(n-1)*100, ' AD, in ', 100, ' yr steps'
  interpolate_insolation = insolation(interpolate_ecc, &
       interpolate_obl, &
       interpolate_lpx, &
       long, lat, S0)

  open(newunit=io, file = 'out/interp_-250_-170.dat', status="replace", action="write")
  do i=1,n
     write(io,*) interpolate_time(i), &
          interpolate_ecc(i), &
          interpolate_obl(i), &
          modulo(interpolate_lpx(i)-pi,2.0_dp*pi), &
          interpolate_insolation(i)
  enddo
  close(io)

  deallocate( &
       interpolate_time, &
       interpolate_ecc, &
       interpolate_obl, &
       interpolate_lpx, &
       interpolate_insolation)

  print *, 'wrote linearly interpolated orbital forcing to file'

  print *,'--------------------------------------------------------------------------------'

  ! calculate insolation for a grid of latitudes and longitudes
  ! this is how we can use it for multiple longitudes and latitudes
  longs = [0._dp, pi / 2._dp, pi, 1.5_dp * pi, 2._dp * pi]
  lats = [-90._dp, -60._dp, -30._dp, 0._dp, 30._dp, 60._dp, 90._dp] / R2D

  call readbindata("dat/PT-ZB18a_1-1.bin", time, ecc, obl, prec, lpx, climprec)
  n=size(time)
  ! for all the timesteps in the astronomical solution

  allocate(latlons(n, 5, 7))

  do i = 1, 5
     do j = 1, 7
        latlons(:,i,j) = insolation(ecc, obl, lpx, longs(i), lats(j), S0)
     end do
  end do
  print *, 'calculated insolation at a grid of solar longitudes'
  print *,' and Earth’s latitudes for all timesteps in astronomical solution'

  !> write only the insolation at t0 to file as a matrix
  open(newunit=io, file='out/insgrid.dat', status = "replace", action = "write")
  do j=1, 7
        write(io,*) latlons(1,1,j), latlons(1,2,j), latlons(1,3,j), latlons(1,4,j), latlons(1,5,j)
  enddo
  close(io)

  deallocate(time,ecc,obl,lpx)
  deallocate(latlons)

  print *, 'wrote the grid of insolation at t0 to out/insgrid.dat'
  ! TODO: write this to a netCDF file instead?
  print *,'--------------------------------------------------------------------------------'

end program paleoinsolation
