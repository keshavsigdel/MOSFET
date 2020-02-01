function fermi(E,mu,kT) !Fermi Function
  implicit none
  real(8) :: fermi
  real(8) :: E  !Energy, main variable
  real(8) :: mu !Chemical potential
  real(8) :: kT ! product of Boltzman constant and absolute temp
  fermi = 1.0d0/(1.0d0+exp((E-mu)/kT))  
end function fermi
