!> The paleoinsolation program makes recent astronomical solutions
!> readily available for use in palaeoclimate simulations, such as the CESM.
!>
!> the Makefile:
!> downloads the astronomical solution ZB18a to ems-plan3.dat
!> available on https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html
!> downloads, patches, and compiles the snvec c-program
!> https://github.com/rezeebe/snvec
!> calculates the PT-solution with Ed = 1 and Td = 1, resulting in out.dat
!> compiles this program and runs it to generate ins.dat for 65°N summer insolation
!> we also show how to linearly interpolate our astronomical solutions
!> and how to generate a grid of Solar longitudes, Earth latitudes
program paleoinsolation
  use kind, only : dp
  use data, only : readdata, writedata
  use orb, only : orbpar
  use insol, only : insolation
  use shr_kind_mod, only : SHR_KIND_R8, SHR_KIND_IN
  use shr_orb_mod, only : shr_orb_params, SHR_ORB_UNDEF_INT
  implicit none

  ! some constants
  ! pi = 3.1415926535897932_dp
  real(dp), parameter :: pi = 4.0_dp*datan(1.0_dp)
  real(dp), parameter :: OMT = 75.594_dp
  real(dp), parameter :: R2D = 180._dp / pi ! radians to degrees


  ! the variables that hold ZB18a(1,1) input
  real(dp), allocatable :: time(:), ecc(:), obl(:), prec(:), lpx(:), climprec(:)
  ! 65°N summer insolation output
  real(dp), allocatable :: sixtyfive(:)
  ! a grid of lat/lon to hold insolation
  real(dp), allocatable :: latlons(:,:,:)
  ! the length of time(:) etc.
  integer :: n
  ! desired true Solar longitude, Earth's latitude, Solar constant
  real(dp) :: long, lat, S0

  ! get orbital parameters at a specific calendar year
  real(dp) :: yearBP
  real(dp) :: yearCE
  ! output a single value in radians
  real(dp) :: ecc1, obl1, lpx1
  ! also in degrees
  real(dp) :: obl1_deg, lpx1_deg

  integer(SHR_KIND_IN) :: iyear_AD
  real(SHR_KIND_R8) :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp

  ! desired grid of Solar longitudes and Earth's latitudes to calcualte insolation for
  integer :: i, j
  real(dp), dimension(5) :: longs
  real(dp), dimension(7) :: lats

  ! the readdata function also allocates these variables
  call readdata("dat/PT-ZB18a_1-1.dat", time, ecc, obl, prec, lpx, climprec)
  print *, 'read snvec ZB18a(1,1) astronomical solution'
  ! alternatively:
  ! call readdata('dat/PT-ZB20a(1,1).dat', time, ecc, obl, prec, lpx, climprec)
  ! print *, 'read snvec ZB20a(1,1) astronomical solution'

  ! re-wrap LPX (it was 'unwrapped' to prevent numerical glitches near the 2 pi jumps)
  lpx = modulo(lpx - pi, 2._dp * pi)
  ! print *, 'wrapped LPX'

  ! interpolate astronomical solution to a single calendar year
  ! set desired year
  !yearBP =    -800000.0_dp ! i.e. -8 Ma = 8 kyr into the future, should throw
  !yearBP =  101000000.0_dp ! i.e. 101 Ma , should throw because we're using 100-0 Ma input solution
  !yearBP =     800000.0_dp ! i.e. 8 ka
  yearBP =   66000000.0_dp ! i.e. 66 Ma should warn
  yearCE = -yearBP + 1950.0_dp
  print *,'using orbpar'
  call orbpar(yearCE,ecc1,obl1,lpx1)
  print *,'linearly interpolated astronomical solution at ',yearBP*1.0e-6_dp,' Ma'
  ! convert from radians to degrees
  obl1_deg = obl1 * R2D
  lpx1_deg = lpx1 * R2D
  print *, 'eccentricity: ', ecc1
  print *, 'obliquity: ', obl1_deg
  print *, 'longitude of perihelion w/ respect to moving equinox: ', lpx1_deg

  yearCE = 1950
  print *,'using orbpar'
  call orbpar(yearCE,ecc1,obl1,lpx1)
print *,'linearly interpolated astronomical solution at ',yearCE,' CE'
  ! convert from radians to degrees
  obl1_deg = obl1 * R2D
  lpx1_deg = lpx1 * R2D
  print *, 'eccentricity: ', ecc1
  print *, 'obliquity: ', obl1_deg
  print *, 'longitude of perihelion w/ respect to moving equinox: ', lpx1_deg

  ! same but implemented with the ESCOMP/CDEPS API
  iyear_AD = yearBP + 1950_SHR_KIND_IN
  call shr_orb_params(iyear_AD,eccen,obliq,mvelp,obliqr,lambm0,mvelpp,.true.) ! the thing that does the work
  print *,'linearly interpolated astronomical solution at ',yearBP*1.0e-6_dp,' Ma'
  print *,'using shr_orb_params'
  print *, 'eccen: ', eccen
  print *, 'obliq: ', obliq
  print *, 'mvelp: ', mvelp
  print *, 'obliqr: ', obliqr
  print *, 'lambm0: ', lambm0
  print *, 'mvelpp: ', mvelpp

  ! a check to see that it still returns correct values
  iyear_AD = SHR_ORB_UNDEF_INT
  eccen = 0.0255_SHR_KIND_R8
  obliq = 23.52866_SHR_KIND_R8
  mvelp = 223.677_SHR_KIND_R8
  call shr_orb_params(iyear_AD, eccen, obliq, mvelp,obliqr,lambm0,mvelpp,.true.)
  print *,'test if shr_orb_params works with pre-set values'
  print *, 'eccen: (0.0255) = ', eccen
  print *, 'obliq: (23.52866) = ', obliq
  print *, 'mvelp: (223.677) = ', mvelp
  print *, 'obliqr: ', obliqr
  print *, 'lambm0: ', lambm0
  print *, 'mvelpp: ', mvelpp

  ! calculate 65°N summer insolation for all timesteps in the astronomical solution
  long = pi / 2._dp
  lat = 65._dp / R2D !pi / 180._dp
  S0 = 1360.7_dp ! the input total insolation

  allocate(sixtyfive(n))
  sixtyfive = insolation(ecc, obl, lpx, long, lat, S0)
  print *, 'calculated 65°N summer insolation at all timesteps in astronomical solution'
  call writedata("dat/ZB18a_insolation.dat", time,ecc,obl,prec,lpx,climprec,sixtyfive)
  print *, 'wrote 65°N summer insolation to file'

  ! calculate insolation for a grid of latitudes and longitudes
  ! this is how we can use it for multiple longitudes and latitudes
  longs = [0._dp, pi / 2._dp, pi, 1.5_dp * pi, 2._dp * pi]
  lats = [-90._dp, -60._dp, -30._dp, 0._dp, 30._dp, 60._dp, 90._dp] / R2D

  n=size(time)
  ! for all the timesteps in the astronomical solution
  allocate(latlons(n, 5, 7))
  do i = 1, 5
     do j = 1, 7
        latlons(:,i,j) = insolation(ecc, obl, lpx, longs(i), lats(j), S0)
     end do
  end do
  print *, 'calculated insolation at a grid of solar longitudes and Earth’s latitudes for all timesteps in astronomical solution'

  !> write only the insolation at t0 to file as a matrix
  open(unit=44, file='insgrid.dat')
  do j=1, 7
        write(44,*) latlons(1,1,j), latlons(1,2,j), latlons(1,3,j), latlons(1,4,j), latlons(1,5,j)
  enddo
  close(44)
  print *, 'wrote the grid of insolation at t0 to insgrid.dat'

  ! TODO: write this to a netCDF file instead?

  deallocate(time, ecc, obl, prec, lpx, climprec, sixtyfive, latlons)

end program paleoinsolation
