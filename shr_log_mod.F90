!BOP ===========================================================================
!
! !MODULE: shr_log_mod -- variables and methods for logging
!
! !DESCRIPTION:
!    Low-level shared variables for logging.
!
!    Also, routines for generating log file messages.
!
! !INTERFACE: ------------------------------------------------------------------

!> IJK:
!> WARNING: this is a trimmed-down version to cut back dependencies
!> so that I could tweak the module shr_orb_mod, do not use this version of the shr_log_mod module!

module shr_log_mod

! !USES:

  use shr_kind_mod, only: shr_kind_in, shr_kind_cx
  !use shr_strconvert_mod, only: toString

  use, intrinsic :: iso_fortran_env, only: output_unit

  implicit none
  private

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

  !public :: shr_log_errMsg
  !public :: shr_log_OOBMsg
  !public :: shr_log_setLogUnit
  public :: shr_log_getLogUnit

! !PUBLIC DATA MEMBERS:

  public :: shr_log_Level
  public :: shr_log_Unit

!EOP

  ! low-level shared variables for logging, these may not be parameters
  integer(SHR_KIND_IN) :: shr_log_Level = 0
  integer(SHR_KIND_IN) :: shr_log_Unit  = output_unit

contains

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_log_errMsg -- Return an error message containing file & line info
!
! !DESCRIPTION:
!     Return an error message containing file & line info
!     \newline
!     errMsg = shr\_log\_errMsg(__FILE__, __LINE__)
!
! This is meant to be used when a routine expects a string argument for some message,
! but you want to provide file and line information.
!
! However: Note that the performance of this function can be very bad. It is currently
! maintained because it is used by old code, but you should probably avoid using this
! in new code if possible.
!
! !REVISION HISTORY:
!     2013-July-23 - Bill Sacks
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_log_getLogUnit(unit)
    integer, intent(out) :: unit

     unit = shr_log_unit

  end subroutine shr_log_getLogUnit

end module shr_log_mod
