module interp
  use kind, only : dp

  implicit none

  private
  public :: locate

contains

! Below interpolation function adapted from
! https://github.com/astrofrog/fortranlib/blob/master/src/lib_array.f90
! repeat of their copyright notice follows:
!
! MD5 of template: 13982be359376d65d4f966553ef27636
! Array related routines (Integration, Interpolation, etc.)
!
! ------------------------------------------------------------------------------
! Copyright (c) 2009-13, Thomas P. Robitaille
!
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ------------------------------------------------------------------------------
! locate_dp
integer function locate(xx,x)
  ! Locate a value in a sorted array
  real(dp), dimension(:), intent(in) :: xx
  real(dp), intent(in) :: x
  integer :: n,jl,jm,ju
  logical :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
     if (ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm
     else
        ju=jm
     end if
  end do

  if (x == xx(1)) then
     locate = 1
  else if (x == xx(n)) then
     locate = n-1
  else if(ascnd.and. (x > xx(n) .or. x < xx(1))) then
     locate = -1
  else if(.not.ascnd.and. (x < xx(n) .or. x > xx(1))) then
     locate = -1
  else
     locate = jl
  end if
end function locate

end module interp
