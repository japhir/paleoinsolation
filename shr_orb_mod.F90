module SHR_ORB_MOD
  use shr_kind_mod, only: SHR_KIND_R8, SHR_KIND_IN
  use shr_log_mod, only: shr_log_unit, shr_log_getLogUnit
  use shr_const_mod, only: shr_const_pi

  ! extra dependencies from this project
  use data, only: readdata
  use interp, only: locate

  implicit none

  public :: shr_orb_params

  real   (SHR_KIND_R8),public,parameter :: SHR_ORB_UNDEF_REAL = 1.e36_SHR_KIND_R8 ! undefined real
  integer(SHR_KIND_IN),public,parameter :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int

  private

  real   (SHR_KIND_R8),parameter :: pi                 = SHR_CONST_PI
  real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MIN  =   0.0_SHR_KIND_R8 ! min value for eccen
  real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MAX  =   0.1_SHR_KIND_R8 ! max value for eccen
  real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0_SHR_KIND_R8 ! min value for obliq
  real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0_SHR_KIND_R8 ! max value for obliq
  real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MIN  =   0.0_SHR_KIND_R8 ! min value for mvelp
  real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MAX  = 360.0_SHR_KIND_R8 ! max value for mvelp

  contains


!> this should be a drop-in replacement for the shr_orb_params function
!> available on https://github.com/ESCOMP/CDEPS/blob/main/share/shr_orb_mod.F90#L236
!> (last accessed on 2024-12-19)
subroutine shr_orb_params( iyear_AD, eccen, obliq, mvelp, &
     & obliqr, lambm0, mvelpp, log_print )
  !-------------------------------------------------------------------------------
  !
  ! Calculate earths orbital parameters by interpolating recent astronomical
  ! solution ZB18a(1,1) from
  !
  !  Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
  !  Paleocene–Eocene boundary age constrained by geology and
  !  astronomy. _Science_, 365(6456), 926–929.
  !  doi:10.1126/science.aax0612
  !  <https://doi.org/10.1126/science.aax0612>.'
  !
  !  Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
  !  astronomical solutions for the Cenozoic era. _Earth and Planetary
  !  Science Letters_. doi:10.1016/j.epsl.2022.117595
  !  <https://doi.org/10.1016/j.epsl.2022.117595>
  !
  ! made easily available in this code in
  !
  ! Kocken, I. J. & Zeebe, R. E. (2025) New Insolation Forcing for Paleoclimate Models
  !
  !------------------------------Code history-------------------------------------
  !
  ! Original Author: Ilja J. Kocken
  ! Date:            2024-12-19
  ! Adaptation of API by Erik Kluzek, 1997
  !
  !-------------------------------------------------------------------------------

  !----------------------------- Arguments ------------------------------------
  integer(SHR_KIND_IN),intent(in)    :: iyear_AD  ! Year to calculate orbit for
  real   (SHR_KIND_R8),intent(inout) :: eccen     ! orbital eccentricity
  real   (SHR_KIND_R8),intent(inout) :: obliq     ! obliquity in degrees
  real   (SHR_KIND_R8),intent(inout) :: mvelp     ! moving vernal equinox long
  real   (SHR_KIND_R8),intent(out)   :: obliqr    ! Earths obliquity in rad
  real   (SHR_KIND_R8),intent(out)   :: lambm0    ! Mean long of perihelion at
  ! vernal equinox (radians)
  real   (SHR_KIND_R8),intent(out)   :: mvelpp    ! moving vernal equinox long
  ! of perihelion plus pi (rad)
  logical             ,intent(in)    :: log_print ! Flags print of status/error

  !------------------------------ Parameters ----------------------------------
!!$  integer(SHR_KIND_IN),parameter :: poblen =47 ! # of elements in series wrt obliquity
!!$  integer(SHR_KIND_IN),parameter :: pecclen=19 ! # of elements in series wrt eccentricity
!!$  integer(SHR_KIND_IN),parameter :: pmvelen=78 ! # of elements in series wrt vernal equinox
!!$  real   (SHR_KIND_R8),parameter :: psecdeg = 1.0_SHR_KIND_R8/3600.0_SHR_KIND_R8 ! arc sec to deg conversion
!!$
  real   (SHR_KIND_R8) :: degrad = pi/180._SHR_KIND_R8   ! degree to radian conversion factor
  real   (SHR_KIND_R8) :: yb4_1950AD         ! number of years before 1950 AD

  character(len=*),parameter :: subname = '(shr_orb_params)'
  integer              :: s_logunit

  !-------------------------- Formats -----------------------------------------
  character(len=*),parameter :: F00 = "('(shr_orb_params) ',4a)"
  character(len=*),parameter :: F01 = "('(shr_orb_params) ',a,i9)"
  character(len=*),parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
  character(len=*),parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

  !-------------------------- data interpolation ------------------------------
  real(SHR_KIND_R8), dimension(:), allocatable :: times, eccs, obls, precs, lpxs, climprecs
  real   (SHR_KIND_R8) :: frac, time_kyr
  integer(SHR_KIND_IN) :: n, ipos

  !----------------------------------------------------------------------------
  ! radinp and algorithms below will need a degree to radian conversion factor
  !
  call shr_log_getLogUnit(s_logunit)
  if ( log_print ) then
     write(s_logunit,F00) 'Calculate characteristics of the orbit:'
  end if

  ! Check for flag to use input orbit parameters

  ! Check for flag to use input orbit parameters

    IF ( iyear_AD == SHR_ORB_UNDEF_INT ) THEN

       ! Check input obliq, eccen, and mvelp to ensure reasonable

       if( obliq == SHR_ORB_UNDEF_REAL )then
          write(s_logunit,F00) trim(subname)//' Have to specify orbital parameters:'
          write(s_logunit,F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
          write(s_logunit,F00) 'iyear_AD is the year to simulate orbit for (ie. 1950): '
          write(s_logunit,F00) 'obliq, eccen, mvelp specify the orbit directly:'
          write(s_logunit,F00) 'The AMIP II settings (for a 1995 orbit) are: '
          write(s_logunit,F00) ' obliq =  23.4441'
          write(s_logunit,F00) ' eccen =   0.016715'
          write(s_logunit,F00) ' mvelp = 102.7'
          ! commented out so that I don't have to introduce dependency on shr_sys_mod and shr_abort_mod
!!$          call shr_sys_abort(subname//' ERROR: unreasonable obliq')

          error stop
       else if ( log_print ) then
          write(s_logunit,F00) 'Use input orbital parameters: '
       end if
       if( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
          write(s_logunit,F03) 'Input obliquity unreasonable: ', obliq
!!$          call shr_sys_abort(subname//' ERROR: unreasonable obliq')
          error stop
       end if
       if( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
          write(s_logunit,F03) 'Input eccentricity unreasonable: ', eccen
          error stop
!!$          call shr_sys_abort(subname//' ERROR: unreasonable eccen')
       end if
       if( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
          write(s_logunit,F03) 'Input mvelp unreasonable: ' , mvelp
          error stop
!!$          call shr_sys_abort(subname//' ERROR: unreasonable mvelp')
       end if

       !eccen2 = eccen*eccen
       !eccen3 = eccen2*eccen

    ELSE  ! Otherwise calculate based on years before present

       if ( log_print ) then
          write(s_logunit,F01) 'Calculate orbit for year: ' , iyear_AD
       end if
       yb4_1950AD = 1950.0_SHR_KIND_R8 - real(iyear_AD,SHR_KIND_R8)
       if ( yb4_1950AD .lt. -100000000.0_SHR_KIND_R8 )then
          write(s_logunit,F00) 'orbit only available for years -100.000.000'
          write(s_logunit,F00) 'Relative to 1950 AD'
          write(s_logunit,F00) 'ZB18a has been verified with Geological data up to 58 Ma.'
          write(s_logunit,F03) '# of years before 1950: ',yb4_1950AD
          write(s_logunit,F01) 'Year to simulate was  : ',iyear_AD
          error stop
!!$          call shr_sys_abort(subname//' ERROR: unreasonable year')
       end if
       if ( yb4_1950AD .lt. -58000000.0_SHR_KIND_R8)then
          write(s_logunit,F00) 'Caution: For ZB18a, the interval -100 Myr to -58 Myr is unconstrained due to solar system chaos.'
       end if

       ! get orbital solution ZB18a(1,1)
       call readdata(times, eccs, obls, precs, lpxs, climprecs)
       n = size(times)
       ! re-wrap
       lpxs = modulo(lpxs - pi, 2.0_SHR_KIND_R8*pi)

       ! the DE431 ephimerides used for the ZB18a solution
       ! has t0 = Julian day 2443144.5003725
       ! this is approximately 1977-01-01 at 00:00:32
       ! so 0 model years = 1977 CE
       ! so to convert from CE to model years, subtract 1977
       ! also: model years are negative, but yearCE is positive

       ! this is the model year in kyr for interpolation
       time_kyr = -(real(iyear_AD,SHR_KIND_R8) - 1977.0_SHR_KIND_R8)*1.0e-3_SHR_KIND_R8

       ipos = locate(times,time_kyr)
       if(ipos == -1) then
          write(s_logunit,F00) 'interpolation of orbital solution failed.'
          write(s_logunit,F03) 'requested time_kyr: ',time_kyr
          write(s_logunit,F03) 'time at t0: ',times(1)
          write(s_logunit,F03) 'time at tfinal',times(n)
          error stop
!!$          call shr_sys_abort(subname//' ERROR: Interpolation out of bounds')
       end if

       frac = ( time_kyr - times(ipos) ) / ( times(ipos+1) - times(ipos) )

       eccen = eccs(ipos) + frac * (eccs(ipos+1) - eccs(ipos))
       obliqr = obls(ipos) + frac * (obls(ipos+1) - obls(ipos))
       ! prec = linear_interpolation(times,precs,time_kyr)
       mvelp = lpxs(ipos) + frac * (lpxs(ipos+1) - lpxs(ipos))

       obliq = obliqr * degrad

       lambm0 = SHR_ORB_UNDEF_REAL
       mvelpp = mvelp + pi

    END IF

end subroutine shr_orb_params


end module shr_orb_mod
