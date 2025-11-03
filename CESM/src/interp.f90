module interp
  use kind, only : dp

  implicit none

  private
  public :: locate

contains

! adapted from https://github.com/astrofrog/fortranlib/blob/master/src/lib_array.f90
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
