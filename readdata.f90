module data
  implicit none
contains

subroutine readdata(time, ecc, obl, prec, climprec, lpx)
  integer :: i !,n
  !!$  n = 249481

  real, dimension(249481) :: time, ecc, obl, prec, climprec, lpx

  open (unit=42,file='out.dat')
  do i=1,249481
  read(42,*) time(i), ecc(i), obl(i), prec(i), climprec(i), lpx(i)
  enddo
  close(42)
end subroutine readdata

!! TODO
! subroutine writedata(time, ecc, obl, prec, climprec, lpx, ins)
! end subroutine writedata

end module data
