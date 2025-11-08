MODULE shr_orb_mod

  use shr_kind_mod, only: SHR_KIND_R8, SHR_KIND_IN
  use shr_sys_mod, only: shr_sys_abort
  use shr_const_mod, only: shr_const_pi
  use shr_log_mod, only: shr_log_getLogUnit
  ! extra dependencies from PaleoInsolation
  use data, only: readdata
  use interp, only: locate

  IMPLICIT none

  !----------------------------------------------------------------------------
  ! PUBLIC: Interfaces and global data
  !----------------------------------------------------------------------------
  public :: shr_orb_azimuth
  public :: shr_orb_cosinc
  public :: shr_orb_cosz
  public :: shr_orb_params
  public :: shr_orb_decl
  public :: shr_orb_print
  public :: set_constant_zenith_angle_deg

  real   (SHR_KIND_R8),public,parameter :: SHR_ORB_UNDEF_REAL = 1.e36_SHR_KIND_R8 ! undefined real
  integer(SHR_KIND_IN),public,parameter :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int

  !----------------------------------------------------------------------------
  ! PRIVATE: by default everything else is private to this module
  !----------------------------------------------------------------------------
  private

  real   (SHR_KIND_R8),parameter :: pi                 = SHR_CONST_PI
  real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MIN  =   0.0_SHR_KIND_R8 ! min value for eccen
  real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MAX  =   0.1_SHR_KIND_R8 ! max value for eccen
  real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0_SHR_KIND_R8 ! min value for obliq
  real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0_SHR_KIND_R8 ! max value for obliq
  real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MIN  =   0.0_SHR_KIND_R8 ! min value for mvelp
  real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MAX  = 360.0_SHR_KIND_R8 ! max value for mvelp

  ! This variable overrides the behavior of shr_orb_cosz() when >=0
  ! this is be set by calling set_constant_zenith_angle_deg()
  real   (SHR_KIND_R8) :: constant_zenith_angle_deg = -1  ! constant, uniform zneith angle [degrees]

  !===============================================================================
CONTAINS
  !===============================================================================

  SUBROUTINE set_constant_zenith_angle_deg(angle_deg)
    real(SHR_KIND_R8),intent(in) :: angle_deg
    constant_zenith_angle_deg = angle_deg
  END SUBROUTINE set_constant_zenith_angle_deg

  !=======================================================================
  real(SHR_KIND_R8) pure function shr_orb_azimuth(jday,lat,lon,declin,z)

    !----------------------------------------------------------------------------
    !
    ! function returns the solar azimuth angle.
    ! azimuth angle is defined with respect to north, positive to east
    ! based on Sproul, Renewable Energy, 2007
    !
    !----------------------------------------------------------------------------
    real   (SHR_KIND_R8),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
    real   (SHR_KIND_R8),intent(in) :: lat    ! Centered latitude  (radians)
    real   (SHR_KIND_R8),intent(in) :: lon    ! Centered longitude (radians)
    real   (SHR_KIND_R8),intent(in) :: declin ! Solar declination  (radians)
    real   (SHR_KIND_R8),intent(in) :: z      ! Solar zenith angle (radians)
    
    real(SHR_KIND_R8) :: hour_angle
    !----------------------------------------------------------------------------
   
    hour_angle = 2.0_SHR_KIND_R8*pi*((jday-floor(jday)) - 0.5_SHR_KIND_R8) + lon

    ! constrain hour_angle to [-pi,pi] to determine east/west below
    if(hour_angle > pi) hour_angle = hour_angle - 2.0_SHR_KIND_R8*pi

    shr_orb_azimuth = (sin(declin)*cos(lat) - cos(declin)*sin(lat)*cos(hour_angle))/ sin(z)
    
    shr_orb_azimuth = max(-1._SHR_KIND_R8, min(shr_orb_azimuth, 1._SHR_KIND_R8))

    shr_orb_azimuth = acos(shr_orb_azimuth)

    ! azimuth is east for times between midnight and noon (h < 0)
    ! azimuth is west for times between noon and midnight (h > 0)
    if(hour_angle > 0.) shr_orb_azimuth = 2.0_SHR_KIND_R8*pi - shr_orb_azimuth
    
  end function shr_orb_azimuth
  
  !=======================================================================
  real(SHR_KIND_R8) pure function shr_orb_cosinc(z,azimuth,beta,aspect)
    
    !----------------------------------------------------------------------------
    !
    ! Returns incidence angle to local surface
    !
    !----------------------------------------------------------------------------
    real   (SHR_KIND_R8),intent(in) :: z            ! Solar zenith angle (radians)
    real   (SHR_KIND_R8),intent(in) :: azimuth      ! Solar azimuth angle (radians)
    real   (SHR_KIND_R8),intent(in) :: beta         ! Surface slope angle (radians)
    real   (SHR_KIND_R8),intent(in) :: aspect       ! Surface azimuth angle (radians)
    
    !----------------------------------------------------------------------------
    
    shr_orb_cosinc = sin(beta)*sin(z)*cos(aspect - azimuth) + cos(beta)*cos(z)
    
    shr_orb_cosinc = max(-1._SHR_KIND_R8, min(shr_orb_cosinc, 1._SHR_KIND_R8))
    
  end function shr_orb_cosinc

  !=======================================================================

  real(SHR_KIND_R8) pure FUNCTION shr_orb_cosz(jday,lat,lon,declin,dt_avg,uniform_angle)

    !----------------------------------------------------------------------------
    !
    ! FUNCTION to return the cosine of the solar zenith angle.
    ! Assumes 365.0 days/year.
    !
    !--------------- Code History -----------------------------------------------
    !
    ! Original Author: Brian Kauffman
    ! Date:            Jan/98
    ! History:         adapted from statement FUNCTION in share/orb_cosz.h
    !
    !----------------------------------------------------------------------------

    real   (SHR_KIND_R8),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
    real   (SHR_KIND_R8),intent(in) :: lat    ! Centered latitude (radians)
    real   (SHR_KIND_R8),intent(in) :: lon    ! Centered longitude (radians)
    real   (SHR_KIND_R8),intent(in) :: declin ! Solar declination (radians)
    real   (SHR_KIND_R8),intent(in), optional   :: dt_avg ! if present and set non-zero, then use in the
    real   (SHR_KIND_R8),intent(in), optional   :: uniform_angle ! if present and true, apply uniform insolation
    ! average cosz calculation
    logical :: use_dt_avg

    !----------------------------------------------------------------------------

    if ( constant_zenith_angle_deg >= 0 ) then
      shr_orb_cosz = cos( constant_zenith_angle_deg * SHR_CONST_PI/180. )
      return
    end if

    if (present(uniform_angle)) then
      shr_orb_cosz = cos(uniform_angle)
      return
    end if

    ! perform the calculation of shr_orb_cosz
    use_dt_avg = .false.
    if (present(dt_avg)) then
      if (dt_avg /= 0.0_shr_kind_r8) use_dt_avg = .true.
    end if
    ! If dt for the average cosz is specified, then call the shr_orb_avg_cosz
    if (use_dt_avg) then
      shr_orb_cosz = shr_orb_avg_cosz(jday, lat, lon, declin, dt_avg)
    else
      shr_orb_cosz = sin(lat)*sin(declin) - cos(lat)*cos(declin) * &
                     cos((jday-floor(jday))*2.0_SHR_KIND_R8*pi + lon)
    end if

  END FUNCTION shr_orb_cosz

  !=======================================================================
  ! A New Algorithm for Calculation of Cosine Solar Zenith Angle
  ! Author: Linjiong Zhou
  ! E-mail: linjiongzhou@hotmail.com
  ! Date  : 2015.02.22
  ! Ref.  : Zhou et al., GRL, 2015
  !=======================================================================

  real (SHR_KIND_R8) pure function shr_orb_avg_cosz(jday, lat, lon, declin, dt_avg)

    use shr_const_mod, only : pi => shr_const_pi

    implicit none

    !-----------------------------------------------------------------------
    ! In/Out Arguements

    real(SHR_KIND_R8), intent(in) :: jday   ! Julian calendar day (1.xx to 365.xx)
    real(SHR_KIND_R8), intent(in) :: lat    ! latitude (radian)
    real(SHR_KIND_R8), intent(in) :: lon    ! longitude (radian)
    real(SHR_KIND_R8), intent(in) :: declin ! solar declination (radian)
    real(SHR_KIND_R8), intent(in) :: dt_avg ! dt for averaged cosz calculation

    !-----------------------------------------------------------------------
    ! Local Arguments

    real(SHR_KIND_R8),parameter :: piover2 = pi/2.0_SHR_KIND_R8
    real(SHR_KIND_R8),parameter :: twopi   = pi*2.0_SHR_KIND_R8

    real(SHR_KIND_R8) :: aa, bb
    real(SHR_KIND_R8) :: del, phi
    real(SHR_KIND_R8) :: cos_h, h
    real(SHR_KIND_R8) :: t1, t2, dt
    real(SHR_KIND_R8) :: tt1, tt2, tt3, tt4

    !-----------------------------------------------------------------------
    ! Compute Half-day Length

    ! adjust latitude so that its tangent will be defined
    if (lat ==  piover2) then
       del = lat - 1.0e-05_SHR_KIND_R8
    else if (lat ==  -piover2) then
       del = lat + 1.0e-05_SHR_KIND_R8
    else
       del = lat
    end if

    ! adjust declination so that its tangent will be defined
    if (declin == piover2) then
       phi = declin - 1.0e-05_SHR_KIND_R8
    else if (declin == -piover2) then
       phi = declin + 1.0e-05_SHR_KIND_R8
    else
       phi = declin
    end if

    ! define the cosine of the half-day length
    ! adjust for cases of all daylight or all night
    cos_h = - tan(del) * tan(phi)
    if (cos_h <= -1.0_SHR_KIND_R8) then
       h = pi
    else if (cos_h >= 1.0_SHR_KIND_R8) then
       h = 0.0_SHR_KIND_R8
    else
       h = acos(cos_h)
    end if

    !-----------------------------------------------------------------------
    ! Define Local Time t and t + dt

    ! adjust t to be between -pi and pi
    t1 = (jday - int(jday)) * twopi + lon - pi

    if (t1 >=  pi) then
       t1 = t1 - twopi
    else if (t1 < -pi) then
       t1 = t1 + twopi
    end if

    dt = dt_avg / 86400.0_SHR_KIND_R8 * twopi
    t2 = t1 + dt

    !-----------------------------------------------------------------------
    ! Compute Cosine Solar Zenith angle

    ! define terms needed in the cosine zenith angle equation
    aa = sin(lat) * sin(declin)
    bb = cos(lat) * cos(declin)

    ! define the hour angle
    ! force it to be between -h and h
    ! consider the situation when the night period is too short
    if (t2 >= pi .and. t1 <= pi .and. pi - h <= dt) then
       tt2 = h
       tt1 = min(max(t1, -h)       ,         h)
       tt4 = min(max(t2, twopi - h), twopi + h)
       tt3 = twopi - h
    else if (t2 >= -pi .and. t1 <= -pi .and. pi - h <= dt) then
       tt2 = - twopi + h
       tt1 = min(max(t1, -twopi - h), -twopi + h)
       tt4 = min(max(t2, -h)        ,          h)
       tt3 = -h
    else
       if (t2 > pi) then
          tt2 = min(max(t2 - twopi, -h), h)
       else if (t2 < - pi) then
          tt2 = min(max(t2 + twopi, -h), h)
       else
          tt2 = min(max(t2 ,        -h), h)
       end if
       if (t1 > pi) then
          tt1 = min(max(t1 - twopi, -h), h)
       else if (t1 < - pi) then
          tt1 = min(max(t1 + twopi, -h), h)
       else
          tt1 = min(max(t1        , -h), h)
       end if
       tt4 = 0.0_SHR_KIND_R8
       tt3 = 0.0_SHR_KIND_R8
    end if

    ! perform a time integration to obtain cosz if desired
    ! output is valid over the period from t to t + dt
    if (tt2 > tt1 .or. tt4 > tt3) then
       shr_orb_avg_cosz = (aa * (tt2 - tt1) + bb * (sin(tt2) - sin(tt1))) / dt + &
            (aa * (tt4 - tt3) + bb * (sin(tt4) - sin(tt3))) / dt
    else
       shr_orb_avg_cosz = 0.0_SHR_KIND_R8
    end if

  end function shr_orb_avg_cosz

  !===============================================================================
  ! New Insolation Forcing for Paleoclimate Studies
  ! Kocken, I.J. and Zeebe, R.E. (2025). ESS Open Archive.
  ! doi: 10.22541/essoar.175511741.18639670
  ! https://github.com/japhir/paleoinsolation
  !===============================================================================
  SUBROUTINE shr_orb_params( iyear_AD, eccen, obliq, mvelp, &
     & obliqr, lambm0, mvelpp, log_print )
    !-------------------------------------------------------------------------------
    ! Calculate earths orbital parameters by interpolating recent astronomical
    ! solution ZB18a(1,1) or ZB20a(1,1).
    !
    ! ZB18a:
    !  Zeebe, R. E., & Lourens, L. J. (2019). Solar System chaos and the
    !  Paleocene–Eocene boundary age constrained by geology and
    !  astronomy. _Science_, 365(6456), 926–929.
    !  doi:10.1126/science.aax0612
    !  <https://doi.org/10.1126/science.aax0612>.'
    !
    ! ZB20a:
    !  Zeebe, R. E. and Lourens, L. J. (2022). Geologically constrained
    !  astronomical solutions for the Cenozoic era. _Earth and Planetary
    !  Science Letters_. doi:10.1016/j.epsl.2022.117595
    !  <https://doi.org/10.1016/j.epsl.2022.117595>
    !
    ! This code is hosted on:
    ! https://github.com/japhir/paleoinsolation
    ! And is permanently archived on:
    ! Kocken, I. J. (2025) paleoinsolation: New Insolation Forcing for Paleoclimate Models
    ! Zenodo. https://doi.org/10.5281/zenodo.17478418
    !
    !------------------------------Code history-------------------------------------
    !
    ! Original Author: Erik Kluzek
    ! Date:            Oct/97
    ! Major overhaul to support ZB18a(1,1) and ZB20a(1,1) orbital solutions:
    ! Author:          Ilja J. Kocken
    ! Date:            2025-11-06
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
    real   (SHR_KIND_R8) :: degrad = pi/180._SHR_KIND_R8   ! degree to radian conversion factor
    real   (SHR_KIND_R8) :: yb4_J2000         ! number of years before J2000.0

    character(len=*),parameter :: subname = '(shr_orb_params)'
    !---------------------------Local variables----------------------------------
    real(SHR_KIND_R8), dimension(:), allocatable :: times, eccs, obls, precs, lpxs, climprecs
    real   (SHR_KIND_R8) :: frac, time_kyr
    integer(SHR_KIND_IN) :: n, ipos
    real   (SHR_KIND_R8) :: beta    ! Intermediate argument for lambm0
    real   (SHR_KIND_R8) :: eccen2  ! eccentricity squared
    real   (SHR_KIND_R8) :: eccen3  ! eccentricity cubed
    integer              :: s_logunit
    !-------------------------- Formats -----------------------------------------
    character(len=*),parameter :: F00 = "('(shr_orb_params) ',4a)"
    character(len=*),parameter :: F01 = "('(shr_orb_params) ',a,i9)"
    character(len=*),parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
    character(len=*),parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

    !----------------------------------------------------------------------------
    ! radinp and algorithms below will need a degree to radian conversion factor
    call shr_log_getLogUnit(s_logunit)
    if ( log_print ) then
       write(s_logunit,F00) 'Calculate characteristics of the orbit:'
    end if

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
          call shr_sys_abort(subname//' ERROR: unreasonable obliq')
       else if ( log_print ) then
          write(s_logunit,F00) 'Use input orbital parameters: '
       end if
       if( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
          write(s_logunit,F03) 'Input obliquity unreasonable: ', obliq
          call shr_sys_abort(subname//' ERROR: unreasonable obliq')
       end if
       if( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
          write(s_logunit,F03) 'Input eccentricity unreasonable: ', eccen
          call shr_sys_abort(subname//' ERROR: unreasonable eccen')
       end if
       if( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
          write(s_logunit,F03) 'Input mvelp unreasonable: ' , mvelp
          call shr_sys_abort(subname//' ERROR: unreasonable mvelp')
       end if
       ! calculate obliquity in radians
       obliqr = obliq*degrad
    ELSE  ! Otherwise calculate based on years before present

       ! this is the model year in kyr for interpolation
       yb4_J2000 = (real(iyear_AD,SHR_KIND_R8) - 2000.0_SHR_KIND_R8)
       time_kyr = yb4_J2000*1.0e-3_SHR_KIND_R8

       if ( log_print ) then
          write(s_logunit,F01) 'Calculate orbit for year: ' , iyear_AD
          write(s_logunit,F03) 'Model time in kyr: ' , time_kyr
       end if

       if ((time_kyr .lt. -300000.0_SHR_KIND_R8) .or. (time_kyr .gt. 0.0_SHR_KIND_R8))then
          write(s_logunit,F00) 'orbit only available for years -300,000,000 to 0'
          write(s_logunit,F00) 'Relative to J2000.0'
          write(s_logunit,F03) '# of years before J2000: ',yb4_J2000
          write(s_logunit,F01) 'Year to simulate was  : ',iyear_AD
          call shr_sys_abort(subname//' ERROR: unreasonable year')
       end if
       ! TODO: you could comment out either of these checks, depending on which solution you use
       if ( time_kyr .lt. -58000.0_SHR_KIND_R8)then
          write(s_logunit,F00) 'Caution: For ZB18a, the interval -300 Myr to -58 Myr is unconstrained due to solar system chaos.'
       end if
       if ( time_kyr .lt. -71000.0_SHR_KIND_R8)then
          write(s_logunit,F00) 'Caution: For ZB20a, the interval -300 Myr to -71 Myr is unconstrained due to solar system chaos.'
       end if
       ! get orbital solution ZB18a(1,1) or ZB20a(1,1)
       ! TODO: update this line to your liking!
       call readdata('dat/PT-ZB18a_1-1.dat', times, eccs, obls, precs, lpxs, climprecs)
       !call readdata('/path/to/my_cesm_sandbox/orb/dat/PT-ZB18a_1-1.dat', times, eccs, obls, precs, lpxs, climprecs) ! update path
       !call readdata('dat/PT-ZB20a_1-1.dat', times, eccs, obls, precs, lpxs, climprecs) ! update solution
       !call readbindata('dat/PT-ZB20a_1-1.bin', times, eccs, obls, precs, lpxs, climprecs) ! use binary files
       n = size(times)

       ! the DE431 ephimerides used for the ZB18a solution
       ! has t0 = Julian day 2443144.5003725 (Folkner et al., 2014)
       ! this is approximately 1977-01-01 at 00:00:32 (web tool conversion)
       ! However, t0 was set to J2000.0 when initial conditions were calculated.
       ! also: model years are negative, but yearCE is positive

       ipos = locate(times,time_kyr)
       if(ipos == -1) then
          write(s_logunit,F00) 'interpolation of orbital solution failed.'
          write(s_logunit,F03) 'requested time_kyr: ',time_kyr
          write(s_logunit,F03) 'time at t0: ',times(1)
          write(s_logunit,F03) 'time at tfinal',times(n)
          call shr_sys_abort(subname//' ERROR: Interpolation out of bounds')
       end if

       frac = ( time_kyr - times(ipos) ) / ( times(ipos+1) - times(ipos) )

       eccen = eccs(ipos) + frac * (eccs(ipos+1) - eccs(ipos))
       obliqr = obls(ipos) + frac * (obls(ipos+1) - obls(ipos))
       ! prec = linear_interpolation(times,precs,time_kyr)
       ! snvec currently provides this in radians, convert to degrees for consistency
       mvelp = (lpxs(ipos) + frac * (lpxs(ipos+1) - lpxs(ipos))) / degrad

       ! Cases to make sure mvelp is between 0 and 360.

       do while (mvelp .lt. 0.0_SHR_KIND_R8)
          mvelp = mvelp + 360.0_SHR_KIND_R8
       end do
       do while (mvelp .ge. 360.0_SHR_KIND_R8)
          mvelp = mvelp - 360.0_SHR_KIND_R8
       end do
       ! calculate obliquity in degrees
       obliq = obliqr / degrad

    END IF  ! end of test on whether to calculate or use input orbital params

    ! 180 degrees must be added to mvelp since observations are made from the
    ! earth and the sun is considered (wrongly for the algorithm) to go around
    ! the earth. For a more graphic explanation see Appendix B in:
    !
    ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
    ! Periods.  J. of Geophysical Research 98:10,341-10,362.
    !
    ! Additionally, orbit will need this value in radians. So mvelp becomes
    ! mvelpp (mvelp plus pi)

    mvelpp = (mvelp + 180._SHR_KIND_R8)*degrad

    ! Set up arguments used several times in lambm0 calculation ahead.

    eccen2 = eccen*eccen
    eccen3 = eccen2*eccen
    beta = sqrt(1._SHR_KIND_R8 - eccen2)

    ! The mean longitude at the vernal equinox (lambda m nought in Berger
    ! 1978; in radians) is calculated from the following formula given in
    ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    ! 1978) is 0.

    lambm0 = 2._SHR_KIND_R8*((.5_SHR_KIND_R8*eccen + .125_SHR_KIND_R8*eccen3)*(1._SHR_KIND_R8 + beta)*sin(mvelpp)  &
         &      - .250_SHR_KIND_R8*eccen2*(.5_SHR_KIND_R8    + beta)*sin(2._SHR_KIND_R8*mvelpp)            &
         &      + .125_SHR_KIND_R8*eccen3*(1._SHR_KIND_R8/3._SHR_KIND_R8 + beta)*sin(3._SHR_KIND_R8*mvelpp))

    if ( log_print ) then
       write(s_logunit,F03) '------ Computed Orbital Parameters ------'
       write(s_logunit,F03) 'Eccentricity      = ',eccen
       write(s_logunit,F03) 'Obliquity (deg)   = ',obliq
       write(s_logunit,F03) 'Obliquity (rad)   = ',obliqr
       write(s_logunit,F03) 'Long of perh(deg) = ',mvelp
       write(s_logunit,F03) 'Long of perh(rad) = ',mvelpp
       write(s_logunit,F03) 'Long at v.e.(rad) = ',lambm0
       write(s_logunit,F03) '-----------------------------------------'
    end if

  END SUBROUTINE shr_orb_params

  !===============================================================================

  SUBROUTINE shr_orb_decl(calday ,eccen ,mvelpp ,lambm0 ,obliqr ,delta ,eccf)

    !-------------------------------------------------------------------------------
    !
    ! Compute earth/orbit parameters using formula suggested by
    ! Duane Thresher.
    !
    !---------------------------Code history----------------------------------------
    !
    ! Original version:  Erik Kluzek
    ! Date:              Oct/1997
    !
    !-------------------------------------------------------------------------------

    !------------------------------Arguments--------------------------------
    real   (SHR_KIND_R8),intent(in)  :: calday ! Calendar day, including fraction
    real   (SHR_KIND_R8),intent(in)  :: eccen  ! Eccentricity
    real   (SHR_KIND_R8),intent(in)  :: obliqr ! Earths obliquity in radians
    real   (SHR_KIND_R8),intent(in)  :: lambm0 ! Mean long of perihelion at the
    ! vernal equinox (radians)
    real   (SHR_KIND_R8),intent(in)  :: mvelpp ! moving vernal equinox longitude
    ! of perihelion plus pi (radians)
    real   (SHR_KIND_R8),intent(out) :: delta  ! Solar declination angle in rad
    real   (SHR_KIND_R8),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)

    !---------------------------Local variables-----------------------------
    real   (SHR_KIND_R8),parameter :: dayspy = 365.0_SHR_KIND_R8  ! days per year
    real   (SHR_KIND_R8),parameter :: ve     = 80.5_SHR_KIND_R8   ! Calday of vernal equinox
    ! assumes Jan 1 = calday 1

    real   (SHR_KIND_R8) ::   lambm  ! Lambda m, mean long of perihelion (rad)
    real   (SHR_KIND_R8) ::   lmm    ! Intermediate argument involving lambm
    real   (SHR_KIND_R8) ::   lamb   ! Lambda, the earths long of perihelion
    real   (SHR_KIND_R8) ::   invrho ! Inverse normalized sun/earth distance
    real   (SHR_KIND_R8) ::   sinl   ! Sine of lmm
    ! Compute eccentricity factor and solar declination using
    ! day value where a round day (such as 213.0) refers to 0z at
    ! Greenwich longitude.
    !
    ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
    ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
    ! 35:2362-2367.
    !
    ! To get the earths true longitude (position in orbit; lambda in Berger
    ! 1978) which is necessary to find the eccentricity factor and declination,
    ! must first calculate the mean longitude (lambda m in Berger 1978) at
    ! the present day.  This is done by adding to lambm0 (the mean longitude
    ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
    ! an increment (delta lambda m in Berger 1978) that is the number of
    ! days past or before (a negative increment) the vernal equinox divided by
    ! the days in a model year times the 2*pi radians in a complete orbit.

    lambm = lambm0 + (calday - ve)*2._SHR_KIND_R8*pi/dayspy
    lmm   = lambm  - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl  = sin(lmm)
    lamb  = lambm  + eccen*(2._SHR_KIND_R8*sinl + eccen*(1.25_SHR_KIND_R8*sin(2._SHR_KIND_R8*lmm)  &
         &     + eccen*((13.0_SHR_KIND_R8/12.0_SHR_KIND_R8)*sin(3._SHR_KIND_R8*lmm) - 0.25_SHR_KIND_R8*sinl)))

    ! Using the obliquity, eccentricity, moving vernal equinox longitude of
    ! perihelion (plus), and earths true longitude, the declination (delta)
    ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
    ! rho will be used), and thus the eccentricity factor (eccf), can be
    ! calculated from formulas given in Berger 1978.

    invrho = (1._SHR_KIND_R8 + eccen*cos(lamb - mvelpp)) / (1._SHR_KIND_R8 - eccen*eccen)

    ! Set solar declination and eccentricity factor

    delta  = asin(sin(obliqr)*sin(lamb))
    eccf   = invrho*invrho

    return

  END SUBROUTINE shr_orb_decl

  !===============================================================================

  SUBROUTINE shr_orb_print( iyear_AD, eccen, obliq, mvelp )

    !-------------------------------------------------------------------------------
    !
    ! Print out the information on the Earths input orbital characteristics
    !
    !---------------------------Code history----------------------------------------
    !
    ! Original version:  Erik Kluzek
    ! Date:              Oct/1997
    !
    !-------------------------------------------------------------------------------

    !---------------------------Arguments----------------------------------------
    integer(SHR_KIND_IN),intent(in) :: iyear_AD ! requested Year (AD)
    real   (SHR_KIND_R8),intent(in) :: eccen    ! eccentricity (unitless)
    ! (typically 0 to 0.1)
    real   (SHR_KIND_R8),intent(in) :: obliq    ! obliquity (-90 to +90 degrees)
    ! typically 22-26
    real   (SHR_KIND_R8),intent(in) :: mvelp    ! moving vernal equinox at perhel
    ! (0 to 360 degrees)
    integer                         :: s_logunit
    logical                         :: debug = .false.
    !-------------------------- Formats -----------------------------------------
    character(len=*),parameter :: F00 = "('(shr_orb_print) ',4a)"
    character(len=*),parameter :: F01 = "('(shr_orb_print) ',a,i9.4)"
    character(len=*),parameter :: F02 = "('(shr_orb_print) ',a,f6.3)"
    character(len=*),parameter :: F03 = "('(shr_orb_print) ',a,es14.6)"
    !----------------------------------------------------------------------------
#ifdef DEBUG
    debug = .true.
#endif
    call shr_log_getLogUnit(s_logunit)
    if(s_logunit .ne. 6 .or. debug) then
       if ( iyear_AD .ne. SHR_ORB_UNDEF_INT ) then
          if ( iyear_AD > 0 ) then
             write(s_logunit,F01) 'Orbital parameters calculated for year: AD ',iyear_AD
          else
             write(s_logunit,F01) 'Orbital parameters calculated for year: BC ',iyear_AD
          end if
       else if ( obliq /= SHR_ORB_UNDEF_REAL ) then
          write(s_logunit,F03) 'Orbital parameters: '
          write(s_logunit,F03) 'Obliquity (degree):              ', obliq
          write(s_logunit,F03) 'Eccentricity (unitless):         ', eccen
          write(s_logunit,F03) 'Long. of moving Perhelion (deg): ', mvelp
       else
          write(s_logunit,F03) 'Orbit parameters not set!'
       end if
    endif

  END SUBROUTINE shr_orb_print
  !===============================================================================

END MODULE shr_orb_mod
