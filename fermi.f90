subroutine fermi(E,mu,kT,f,n) !Fermi Function
  implicit none
  integer(2) :: i, n
  real(8) :: E(n)  !Energy, main variable
  real(8) :: mu !Chemical potential
  real(8) :: kT ! product of Boltzman constant and absolute temp
  real(8) :: f(n)
  do i= 1, n
    f(i) = 1.0d0/(1.0d0+exp((E(i)-mu)/kT))  
  end do
end subroutine fermi
