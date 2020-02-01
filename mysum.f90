function mysum(x,n)
  implicit none
  integer(2) :: i,n
  real(8) :: x(n), mysum
  mysum = 0.0d0
  do i = 1,n
    mysum =  mysum +x(i)
  end do
end function mysum

