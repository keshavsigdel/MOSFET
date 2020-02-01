program main
  implicit none
  real(8) :: a, b
  integer(2) :: n, i, j, IV
  real(8) :: D0
  real(8), allocatable :: E(:), D(:), VV(:), II(:), NN(:), f0(:),f1(:), f2(:)
  real(8), parameter :: pi = 3.141590d0
  real(8), parameter :: hbar = 1.055d-34
  real(8), parameter :: q = 1.602d-19
  real(8), parameter :: eps0 = 8.854d-12
  real(8), parameter :: epsr = 4
  real(8), parameter :: m = 0.25d0*9.11d-31
  real(8) :: v,I0, dE, N0, N1,N2, Vg, Vd
  real(8), external :: mysum 
  !parameters
  real(8) :: W = 1.0d-06
  real(8) :: L = 10.0d-09
  real(8) :: t = 1.50d-09
  real(8) :: Cg, Cs, Cd, CE, U0, alphag, alphad, kT, mu,mu1,mu2, ep, g1, g2, g
  real(8) :: U, UL, dU, Unew
  Cg = epsr*eps0*W*L/t
  Cs = 0.05d0*Cg
  Cd = 0.05d0*Cg
  U0 = q/CE
  alphag = Cg/CE
  alphad = Cd/CE
  kT = 0.025d0
  mu = 0.0d0
  ep = 0.2d0
  v =1.0d+05
  g1 = hbar*v/(q*L)
  g2 = g1
  g = g1+g2
  
! Energy grid
  a = -1.0d0
  b = 1.0d0
  n = 501
  I0 = q*q/hbar
  dE = E(2)-E(1)
  open(1, file='IV.dat',action='write')
  allocate(E(n),D(n),VV(n),II(n), NN(n),f0(n),f1(n),f2(n))
  call linspace(a,b,n,E)
  D0 = m*q*W*L/(pi*hbar*hbar) ! step density of states per eV
  D = 0.0d0
  do i = (n+1)/2+1, n
    D(i) = 1.0d0
  end do
  D=D0*D

  ! refereence number of electrons
  call fermi(E+ep,mu,kT,f0,n)
  N0 = 2.0d0*dE*mysum(D*f0,n)
  !Bias
  IV = 61
  allocate(VV(IV))
  call linspace(0.0d0, 0.6d0, IV, VV)
  do i = 1, IV
    Vg = 0.5d0
    Vd = VV(i)
    mu1 = mu
    mu2 = mu1-Vd
    UL = -(alphag*Vg)-(alphad*Vd)
    U = 0.0d0
    dU = 1.0d0
   
   do while(dU .gt. 1.0d-06)
     call fermi(E+UL+U+ep, mu1, kT, f1, n)
     call fermi(E+UL+U+ep, mu2, kT, f2, n)
     NN(i) = dE*mysum(D*f1*g1/g+f2*g2/g,n)
     Unew = U0*(NN(i)-N0)
     dU = abs(U-Unew)
     U = U+0.10d0*(Unew-U)
   end do
   II(i) = dE*I0*mysum(D*(f1-f2))*g1*g2/g
  write(1,*) VV(i),II(i)
end do
  close(1)
deallocate(E,D,VV,II,NN,f0,f1,f2)
end program main


