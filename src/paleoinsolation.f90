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
  use shr_orb_mod, only : shr_orb_params, SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
  implicit none

  ! some constants
  ! pi = 3.1415926535897932_dp
  real(dp), parameter :: pi = 4.0_dp*datan(1.0_dp)
  real(dp), parameter :: OMT = 75.594_dp
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

  ! the length of time(:) etc.
  integer :: i, n, io, interpolate_n
  ! desired true Solar longitude, Earth's latitude, Solar constant
  ! input parsing
  real(dp) :: time_start_myr, time_end_myr, time_res_kyr
  real(dp) :: time_start_kyr, time_end_kyr
  ! TODO: implement solution? But then orbpar no longer has the same interface...
!!$  character(len=5) :: orbital_solution = "ZB18a"
  ! in degrees
  real(dp) :: longr, latr, S0
  ! converted to in radians
  real(dp) :: long, lat

  ! get orbital parameters at a specific calendar year
  real(dp) :: yearBP
  real(dp) :: yearCE

  print *,'--------------------------------------------------------------------------------'
  print *, 'reading input file input.txt'
  open(newunit=io, file="input.txt", status="old", action="read")
  read(io, *) time_start_myr, time_end_myr, time_res_kyr, longr, latr, S0!, orbital_solution
  close(io)
  print *, 'Time start (Myr in the past, t < 0): ', time_start_myr, '(Myr)'
  print *, 'Time end (Myr in the past, t <= 0):  ', time_end_myr, '(Myr)'
  print *, 'Time resolution:                     ', time_res_kyr, '(kyr)'
  print *, 'True solar longitude:                ', longr, '(°)'
  print *, 'Earth’s latitude:                    ', latr, '(°)'
  print *, 'S0:                                  ', S0, '(Wm¯²)'
!!$  print *, 'Solution:                               ', orbital_solution 

  ! convert to radians
  long = longr / R2D
  lat = latr / R2D
  time_start_kyr = time_start_myr * 1e3_dp
  time_end_kyr = time_end_myr * 1e3_dp

  print *,'--------------------------------------------------------------------------------'

  if (time_res_kyr < 0.0_dp) then
     print *, 'time_res_kyr must be >= 0.0'
     error stop
  else if (time_res_kyr > 0.0_dp) then ! resample to desired time resolution
     ! calculate number of steps
     interpolate_n = int(abs(time_start_kyr - time_end_kyr) / time_res_kyr)
     print *, 'Number of timesteps:     ', interpolate_n
     allocate( &
          interpolate_time(interpolate_n), &
          interpolate_ecc(interpolate_n), &
          interpolate_obl(interpolate_n), &
          interpolate_lpx(interpolate_n), &
          interpolate_insolation(interpolate_n) &
          )

     do i=1,interpolate_n
        interpolate_time(i) = time_end_kyr - (i-1)*time_res_kyr
        yearCE = interpolate_time(i) * 1e3_dp + 2000.0_dp
        call orbpar(yearCE,interpolate_ecc(i),interpolate_obl(i),interpolate_lpx(i))
     end do
     interpolate_insolation = insolation( &
          interpolate_ecc, &
          interpolate_obl, &
          interpolate_lpx, &
          long, lat, S0)
     print *, 'linearly interpolated solution'
     open(newunit=io, file = 'out/ZB18a_insolation.dat', status="replace", action="write")
     do i=1,interpolate_n
        write(io,*) interpolate_time(i), interpolate_ecc(i), interpolate_obl(i)*R2D, &
             modulo(interpolate_lpx(i)-pi,2.0_dp*pi)*R2D, interpolate_insolation(i)
     enddo
     close(io)
     print *, 'wrote linearly interpolated orbital forcing to out/ZB18a_insolation.dat'
     print *, 'columns are:'
     print *, '- time (kyr)'
     print *, '- eccentricity (-)'
     print *, '- obliquity (°)'
     print *, '- longitude of perihelion with respect to the moving equinox (°)'
     print *, '- insolation for specified input parameters (Wm¯²)'

     deallocate( &
          interpolate_time, &
          interpolate_ecc, &
          interpolate_obl, &
          interpolate_lpx, &
          interpolate_insolation)
  else ! if time_res_kyr is 0, we take the times in the solution

     ! print nice errors if time outside orbital solution is specified
     if((time_start_kyr .lt. -300.0e3_dp) .or. (time_end_kyr .gt. 0.0_dp)) then
        print *, 'ERROR: Orbital solution is provided between -300 Myr and 0 Myr.'
        print *, 'time_start_kyr = ',time_start_kyr
        print *, 'time_end_kyr = ',time_end_kyr
        error stop
     end if

     ! print warning if time is older then well-constrained eccentricity
     if(time_start_kyr .lt. -58.0e3_dp .or. time_end_kyr .lt. -58.0e3_dp) then
        print *, 'Caution: For ZB18a, the interval -300 Myr to -58 Myr'
        print *, 'is unconstrained due to solar system chaos.'
       !!$ print *, 'Caution: For ZB20a, the interval -300 Myr to -71 Myr is unconstrained due to solar system chaos.'
     end if

     ! the readdata function also allocates these variables
     ! so make sure to deallocate at the end!
!!$     call readdata("dat/PT-ZB18a_1-1.dat", time, ecc, obl, prec, lpx, climprec)
     ! instead, we can read from binary files if speed is important
     call readbindata("dat/PT-ZB18a_1-1.bin", time, ecc, obl, prec, lpx, climprec)
     ! alternatively, for ZB20a
     ! call readdata('dat/PT-ZB20a_1-1.dat', time, ecc, obl, prec, lpx, climprec)
     ! print *, 'read snvec ZB20a(1,1) astronomical solution'

     ! note that the longitude of perihelion with respect to the moving equinox is
     ! "unwrapped", in order to prevent interpolation artefacts!
     ! furthermore, it is 180° off from what we need for the insolation calculation.
     ! the functions orbpar takes care of this conversion!
     n = size(time)
     allocate(interpolate_insolation(n))

     interpolate_insolation = insolation( &
          ecc, &
          obl, &
          modulo(lpx - pi, 2.0_dp*pi), &
          long, lat, S0)
     print *, 'calculated insolation for each timestep in the solution'
     open(newunit=io, file = 'out/ZB18a_insolation.dat', status="replace", action="write")
     do i=1,n
        if ((time(i) < time_start_kyr) .or. (time(i) > time_end_kyr)) then
           cycle ! only output desired timesteps
        end if
        write(io,*) time(i), ecc(i), obl(i)*R2D, lpx(i)*R2D, interpolate_insolation(i)
     enddo
     close(io)

     print *, 'wrote orbital forcing to out/ZB18a_insolation.dat'
     print *, 'columns are:'
     print *, '- time (kyr)'
     print *, '- eccentricity (-)'
     print *, '- obliquity (°)'
     print *, '- longitude of perihelion with respect to the moving equinox (°)'
     print *, '- insolation for specified input parameters (Wm¯²)'

     deallocate(time, ecc, obl, prec, lpx, climprec)
     deallocate(interpolate_insolation)
  end if ! done with main program

  print *,'--------------------------------------------------------------------------------'

end program paleoinsolation
