module CF_int
!Purpose: 
!         1) Calculation of dN/dy and dN/dp_t using Cooper-Frye formula with several kinds of distributions 
!            selected by options: Fermi, Bose, Boltzmann and also all of them * theta(p^mu dSigma_mu)
!         2) dN/dy and dN/dp_t in histograms
!Uses: 
!         gauss_legendre module to calculate weights and points in gauss-legendre integration
!         Histogram_module
!Author:
!         Dmytro Oliinychenko (oliiny@fias.uni-frankfurt.de)

use Histogram_module

implicit none



 !Number of integration points in one dimension
 integer, parameter :: nG = 36

 !Parameters for Gauss-Legendre integration
 double precision, parameter, private :: pi = 3.141592653589793238462643d0
 double precision, parameter, private :: phil = 0.d0, phiu= 2.d0*pi
 double precision, parameter, private :: xl = 0.d0, xu= 0.5d0*pi
 double precision, parameter, private :: tl = -0.5d0*pi, tu= 0.5d0*pi
 double precision, dimension(nG), private :: xpts,phipts, xw, phiw, tpts, tw



contains

!===========================================================================================
subroutine init_CF_int()
 use gauss_legendre
 implicit none

 !Get points and coefficients for Gauss-Legendre integration
 call gauss_legendre_quadrature(nG,phil,phiu,phipts,phiw)
 call gauss_legendre_quadrature(nG,xl,xu,xpts,xw)
 call gauss_legendre_quadrature(nG,tl,tu,tpts,tw)

end subroutine

!===========================================================================================

subroutine get_dNdy_histo(hist, m, T, mu, u, dSigma, t_flag, bose_fermi_boltz, degen)
!Creates histogram with dNdy distribution
 implicit none
 type (histo) hist
 double precision, intent(in) :: m, T, mu, u(0:3), dSigma(0:3), bose_fermi_boltz, degen
 logical, intent(in) :: t_flag
 double precision ybin
 integer i

 do i=-hist%Nbin, hist%Nbin
  ybin = BinToValue(hist,i)
  hist%h(i) = dNdy(ybin, m, T, mu, u, dSigma, t_flag, bose_fermi_boltz, degen)
 end do

end subroutine

!===========================================================================================

subroutine get_dNptdpt_histo(hist, m, T, mu, u, dSigma, t_flag, bose_fermi_boltz, degen)
!Creates histogram with dNdy distribution
 implicit none
 type (histo) hist
 double precision, intent(in) :: m, T, mu, u(0:3), dSigma(0:3), bose_fermi_boltz, degen
 logical, intent(in) :: t_flag
 double precision ptbin
 integer i

 do i=-hist%Nbin, hist%Nbin
  ptbin = BinToValue(hist,i)
  hist%h(i) = dNptdpt(ptbin, m, T, mu, u, dSigma, t_flag, bose_fermi_boltz, degen)
 end do

end subroutine

!===========================================================================================
double precision function dNdy(y, m, T, mu, u, dSigma, t_flag, bose_fermi_boltz, degen) result (res)
 implicit none
 ! y - rapidity
 ! m - mass [GeV]
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! Returns dN/dy for pions
 ! assuming Juttner distribution function f = 1/(exp(p^mu u_mu / T) - 1) if t_flag is .FALSE.
 ! assuming Bugaev  distribution function f = 1/(exp(p^mu u_mu / T) - 1)*Theta(p^mu * dSigma_mu) if t_flag is .TRUE.
 ! bose_fermi_boltz = 1 (Fermi distribution), 0 (Boltzmann), -1 (Bose)
 ! degen - degeneracy = (2S+1)(2I+1)
 double precision, intent(in) :: y, m, T, mu, u(0:3), dSigma(0:3)
 logical, intent(in) :: t_flag
 double precision intres, pt, mt, pmu(0:3), pmuumu, pmusigmamu, hbarc, bose_fermi_boltz, degen
 integer i,j

 intres=0.d0
 do i=1,nG; do j=1,nG
  pt=tan(xpts(j))                  !x integration variable: tg x = p_T/T
  mt = sqrt(pt**2 +(m/T)**2)       !m_T/T
  pmu(0) = mt*cosh(y)              !p(0)/T
  pmu(1) = pt*cos(phipts(i))       !p(1)/T
  pmu(2) = pt*sin(phipts(i))       !p(2)/T
  pmu(3) = mt*sinh(y)              !p(3)/T
  pmusigmamu = pmu(0)*dSigma(0) + pmu(1)*dSigma(1) + pmu(2)*dSigma(2) + pmu(3)*dSigma(3) !p^mu*dSigma_mu/T
  pmuumu = pmu(0)*u(0) + pmu(1)*u(1) + pmu(2)*u(2) + pmu(3)*u(3)                         !p^mu*u_mu/T
  if ((.not.t_flag) .OR. pmusigmamu > 0.d0) then
   intres = intres + phiw(i)*xw(j)* &      !Gauss weights
            (1.d0 + pt**2)* &              !d(p_T/T) = (1+(tg x)^2) dx
            pt*&                           !d^3p/p0/dy
            pmusigmamu/(exp(pmuumu - mu/T) + bose_fermi_boltz) !integrated function itself
  endif
 end do; end do

 hbarc = 0.19732d0 !hbar*c [GeV*fm]

 res = degen* (T/(2.d0*pi*hbarc))**3 * intres

end function dNdy

!======================================================================================================
double precision function dNptdpt(pT, m, T, mu, u, dSigma, t_flag, bose_fermi_boltz, degen) result (res)
 implicit none
 ! pT - transverse momentum
 ! m - mass [GeV]
 ! T - temperature [GeV]
 ! u - 4-velocity with lower index
 ! dSigma - normal 4-vector with lower index [fm^3]
 ! Returns dN/dy for pions
 ! assuming Juttner distribution function f = 1/(exp(p^mu u_mu / T) - 1) if t_flag is .FALSE.
 ! assuming Bugaev  distribution function f = 1/(exp(p^mu u_mu / T) - 1)*Theta(p^mu * dSigma_mu) if t_flag is .TRUE.
 ! bose_fermi_boltz = 1 (Fermi distribution), 0 (Boltzmann), -1 (Bose)
 ! degen - degeneracy = (2S+1)(2I+1)
 double precision, intent(in) :: pT, m, T, mu, u(0:3), dSigma(0:3)
 logical, intent(in) :: t_flag
 double precision intres, y, mt, pmu(0:3), pmuumu, pmusigmamu, hbarc, bose_fermi_boltz, degen
 integer i,j

 intres=0.d0
 do i=1,nG; do j=1,nG
  y=tan(tpts(j))                         !t integration variable: tg t = y
  mt = sqrt(pt**2 +m**2)                 !m_T
  pmu(0) = mt*cosh(y)                  !p(0)
  pmu(1) = pt*cos(phipts(i))           !p(1)
  pmu(2) = pt*sin(phipts(i))           !p(2)
  pmu(3) = mt*sinh(y)                  !p(3)
  pmusigmamu = pmu(0)*dSigma(0) + pmu(1)*dSigma(1) + pmu(2)*dSigma(2) + pmu(3)*dSigma(3) !p^mu*dSigma_mu/T
  pmuumu = pmu(0)*u(0) + pmu(1)*u(1) + pmu(2)*u(2) + pmu(3)*u(3)                         !p^mu*u_mu/T
  if ((.not.t_flag) .OR. pmusigmamu > 0.d0) then
   intres = intres + phiw(i)*tw(j)* &      !Gauss weights
            (1.d0 + y**2)* &               !dy = (1+(tg t)^2) dt = (1+y^2)dt
            pmusigmamu/(exp(pmuumu/T - mu/T) + bose_fermi_boltz) !integrated function itself
  endif
 end do; end do

 hbarc = 0.19732d0 !hbar*c [GeV*fm]

 res = degen /(2.d0*pi*hbarc)**3 * intres

end function dNptdpt


end module CF_int
