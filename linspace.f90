subroutine linspace(a,b,n,E)
  implicit none
  real(8) :: a, b, h           ! input variables
  integer(2) :: i
  integer(2) :: n
  real(8) :: E(n)
! For array E
  h = abs(b-a)/dble(n-1)
  do i =1,n
     E(i)=a+h*dble(i-1)
  end do
end subroutine linspace

    

