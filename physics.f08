! wiscobolt: A deterministic, finite-element photon-electron Boltzmann transport solver
! Copyright (C) 2023.  Muhsin H. Younis
!
! wiscobolt is free software: you can redistribute it and/or modify it under the terms of the
! GNU General Public License as published by the Free Software Foundation, either version 3 of
! the License, or (at your option) any later version. wiscobolt is distributed in the hope that
! it will be useful, but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose. See the GNU General Public License,
! located at LICENSE/gpl-3.0.txt in the main wiscobolt directory or https://www.gnu.org/licenses/
! for more details.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

module physics
    use math
implicit none
public
real, parameter :: kB = 1.0 ! J/K
real, parameter :: lightspeed = 2.99792458E+10 ! cm/s
real, parameter :: elem_q = 1.60217663E-19 ! C
real, parameter :: epsilon_0 = 8.85418781E-12 ! F/m
real, parameter :: e_mass = 9.1093837E-28 ! g
real, parameter :: e_mass_E = 0.51099895 ! MeV
real, parameter :: e_radius = 2.817941E-13 ! cm
real, parameter :: n_mass = 1.6749275E-24 ! g
real, parameter :: n_mass_E = 1. ! MeV
real, parameter :: p_mass = 1.6726219E-24 ! g
real, parameter :: p_mass_E = 1. ! MeV
real, parameter :: fsa = 0.0072973525693
real, parameter :: Avogadro = 6.0221408E+23
real, parameter :: Planck = 6.62607015E-34 ! J*s
real, parameter :: hbar = Planck/twopi ! J*s/rad
real, parameter :: Bohr = 5.29177210903E-9 ! cm

! Unit conversions
real, parameter :: m_to_cm = 1.0E+2
real, parameter :: g_to_kg = 1.0E-3
real, parameter :: g_to_amu = 1.66053906E-24
real, parameter :: g_to_MeV = 1.

! Public variables for functions
real :: XSa, XSb, XSc, XSd, XSe, XScut, XSbeta_U2, XSbeta_W2, XSw, XSDelta, XSEgp, XSEg, XSR, XSmu, XSx, XSI
real :: XSBP, XSCP, XSx1, XSx0, XSeta
integer :: XSlp, XSk, XSZ
logical :: XSlogic

contains
real function beta(T) ! NOTE: THIS IS FOR ELECTRON ONLY. SHOULD MAKE IT beta(T,M)
    implicit none
    real, intent(in) :: T ! MeV

    beta = sqrt(1-(T/e_mass_E+1)**(-2))
end function beta

real function unitless_e_beta(x)
    implicit none
    real, intent(in) :: x ! Electron rest energies

    unitless_e_beta = sqrt(1-(x+1)**(-2))
end function unitless_e_beta

real function Boltzmann_factor(E1,E2,T)
    implicit none
    real, intent(in) :: E1 ! J
    real, intent(in) :: E2 ! J
    real, intent(in) :: T ! K

    Boltzmann_factor = exp((E1-E2)/(kB*T))
end function Boltzmann_factor

real function Partition_function(E,T)
    implicit none
    real, dimension(:), intent(in) :: E ! J
    real, intent(in) :: T ! K

    Partition_function = sum(exp(-E/(kB*T)))
end function Partition_function

! XS functions
! ----- Compton
real function unitless_sigmapct(x)
    implicit none
    real, intent(in) :: x ! Energy in electron rest masses of incident photon

    unitless_sigmapct = (2*(x+1)**2/(x**2*(2*x+1)))+(0.5/x-(x+1)/x**3)*log(1+2*x)-(3*x+1)/(2*x+1)**2
end function unitless_sigmapct

real function dsigmapcdE(Eprime,E)
! Incoming energy: Eprime
! Outgoing energy: E
    implicit none
    real, intent(in) :: Eprime ! MeV
    real, intent(in) :: E ! MeV

    dsigmapcdE = (pi*(e_radius**2)*e_mass_E/(Eprime**2))*(E/Eprime+Eprime/E + &
    2*e_mass_E*(1/Eprime-1/E)+(e_mass_E**2)*(1/Eprime-1/E)**2)
end function dsigmapcdE

real function Eminpc(Eprime)
! Incoming energy: Eprime
    implicit none
    real, intent(in) :: Eprime ! MeV

    Eminpc = Eprime/(1+2*Eprime/e_mass_E)
end function Eminpc

real function Mupc(Eprime,E)
! Incoming energy: Eprime
! Outgoing energy: E
    implicit none
    real, intent(in) :: Eprime ! MeV
    real, intent(in) :: E ! MeV

    Mupc = 1+e_mass_E*(1/Eprime-1/E)
end function Mupc

real function unitless_dsigmapcdE(x,y)
    implicit none
    real, intent(in) :: x
    real, intent(in) :: y

    unitless_dsigmapcdE = (x/y+y/x+2*(1/x-1/y)+(1/x-1/y)**2)/(x**2)
end function unitless_dsigmapcdE

real function unitless_Eminpc(x)
    implicit none
    real, intent(in) :: x

    unitless_Eminpc = x/(1+2*x)
end function unitless_Eminpc

real function unitless_Mupc(x,y)
    implicit none
    real, intent(in) :: x
    real, intent(in) :: y

    unitless_Mupc = 1+1/x-1/y
end function unitless_Mupc

real function MGXS_photon_photon_Compton_integrand(u,v)
    implicit none
    ! Requires manually defining public variables:
    ! XSa - lower bound of g' group
    ! XSb - upper bound of g' group
    ! XSc - lower bound of g group
    ! XSd - upper bound of g group
    ! XSlp - degree of mu
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,x),XSc)
    l_gamma = min(max(XSc,unitless_Eminpc(x)),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_photon_photon_Compton_integrand = &
    0.25*(XSb-XSa)*(l_delta-l_gamma)*&
    (unitless_Mupc(x,y)**XSlp)*&
    unitless_dsigmapcdE(x,y)
end function MGXS_photon_photon_Compton_integrand

real function FEXS_photon_photon_Compton_integrand(u,v)
    implicit none
    ! Requires manually defining public variables:
    ! XSa - E^g'_2/M
    ! XSb - E^g'_1/M
    ! XSc - E^g_2/M
    ! XSd - E^g_1/M
    ! XSEgp - in Lambda, this is E^g'_(3-n')/M or E_(g'-n'+2)/M
    ! XSEg - in Lambda, this is E^g_(3-n)/M or E_(g-n+2)/M
    ! XSlp - degree of mu
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,x),XSc)
    l_gamma = min(max(XSc,unitless_Eminpc(x)),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    FEXS_photon_photon_Compton_integrand = &
    0.25*(XSEgp-x)*(XSEg-y)*&
    (XSb-XSa)*(l_delta-l_gamma)*&
    (unitless_Mupc(x,y)**XSlp)*&
    unitless_dsigmapcdE(x,y)
end function FEXS_photon_photon_Compton_integrand

real function FEXS_photon_total_Compton_integrand(u)
    implicit none
    ! Requires manually defining public variables:
    ! XSa - E_g/M
    ! XSb - E_(g+1)/M
    ! XSEgp - in Lambda, this is E^g'_(3-n')/M or E_(g'-n'+2)/M
    ! XSEg - in Lambda, this is E^g_(3-n)/M or E_(g-n+2)/M
    real, intent(in) :: u

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)

    FEXS_photon_total_Compton_integrand = (XSEgp-x)*(XSEg-x)*&
    0.5*(XSb-XSa)*&
    unitless_sigmapct(x)
end function FEXS_photon_total_Compton_integrand

real function Emaxec(Eprime)
! Incoming (photon) energy: Eprime
    implicit none
    real, intent(in) :: Eprime

    Emaxec = 2*Eprime**2/(e_mass_E+2*Eprime)
end function Emaxec

real function Muec(Eprime,Ee)
! Incoming (photon) energy: Eprime
! Outgoing (electron) energy: Ee
    implicit none
    real, intent(in) :: Eprime
    real, intent(in) :: Ee

    Muec = (1+e_mass_E/Eprime)/sqrt(1+2*e_mass_E/Ee)
end function Muec

real function unitless_Emaxec(x)
    implicit none
    real, intent(in) :: x ! Incoming (photon) energy

    unitless_Emaxec = 2*x**2/(1+2*x)
end function unitless_Emaxec

real function unitless_Muec(x,y)
    implicit none
    real, intent(in) :: x ! Incoming (photon) energy
    real, intent(in) :: y ! Outgoing (electron) energy

    unitless_Muec = (1+1/x)/sqrt(1+2/y)
end function unitless_Muec

real function MGXS_photon_electron_Compton_integrand(u,v)
    implicit none
    ! Requires manually defining public variables:
    ! XSa - lower bound of g' group
    ! XSb - upper bound of g' group
    ! XSc - lower bound of g group
    ! XSd - lower bound of g group
    ! XSlp - degree of mu
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,unitless_Emaxec(x)),XSc)
    y = 0.5*(l_delta-XSc)*v + 0.5*(l_delta+XSc)

    MGXS_photon_electron_Compton_integrand = &
    0.25*(XSb-XSa)*(l_delta-XSc)*&
    (unitless_Muec(x,y)**XSlp)*&
    unitless_dsigmapcdE(x,x-y)
end function MGXS_photon_electron_Compton_integrand

real function FEXS_photon_electron_Compton_integrand(u,v)
    implicit none
    ! Requires manually defining public variables:
    ! XSa - E^g'_2/M
    ! XSb - E^g'_1/M
    ! XSc - E^g_2/M
    ! XSd - E^g_1/M
    ! XSEgp - in Lambda, this is E^g'_(3-n')/M or E_(g'-n'+2)/M
    ! XSEg - in Lambda, this is E^g_(3-n)/M or E_(g-n+2)/M
    ! XSlp - degree of mu
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,unitless_Emaxec(x)),XSc)
    y = 0.5*(l_delta-XSc)*v + 0.5*(l_delta+XSc)

    FEXS_photon_electron_Compton_integrand = &
    0.25*(XSEgp-x)*(XSEg-y)*&
    (XSb-XSa)*(l_delta-XSc)*&
    (unitless_Muec(x,y)**XSlp)*&
    unitless_dsigmapcdE(x,x-y)
end function FEXS_photon_electron_Compton_integrand

! ----- Inelastic Electron
real function inel_mu(x,eps)
    implicit none
    real, intent(in) :: x ! incident electron energy
    real, intent(in) :: eps ! outgoing electron energy

    inel_mu = sqrt(eps*(x+2)/(x*(eps+2)))
end function inel_mu

real function TEMP_inel_mu_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,x),XSc)
    l_gamma = min(max(XSc,x/2),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    TEMP_inel_mu_integrand = 0.25*(XSb-XSa)*(l_delta-l_gamma)*&
        inel_mu(x,y)
end function TEMP_inel_mu_integrand

real function TEMP_inel_mu_area_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,x),XSc)
    l_gamma = min(max(XSc,x/2),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    TEMP_inel_mu_area_integrand = 0.25*(XSb-XSa)*(l_delta-l_gamma)
end function TEMP_inel_mu_area_integrand

real function inel_zmin(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g'+2
    ! XSlogic: True if FEXS, false if MGXS

    inel_zmin = merge(XScut, x-XScut, XSlogic)
end function inel_zmin

real function inel_zmax(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g'+2
    ! XSw: The shell binding energy
    ! XSlogic: True if FEXS, false if MGXS

    inel_zmax = merge(max(0.5*(x-XSw), XScut), max(0.5*(x-XSw), x-XScut), XSlogic)
end function inel_zmax

real function inel_ymin(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! XSw: The shell binding energy

    inel_ymin = 0.5*(x-XSw)
end function inel_ymin

real function inel_ymax(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g'+2
    ! XSw: The shell binding energy

    inel_ymax = merge(max(x-XSw-XScut,0.5*(x-XSw)), max(XScut-XSw,0.5*(x-XSw)), XSlogic)
end function inel_ymax

real function unitless_Moller(x,z)
    implicit none
    real, intent(in) :: x ! incident electron energy
    real, intent(in) :: z ! outgoing electron energy

    unitless_Moller = &
        (1/z**2+1/(x-z)**2+1/(x+1)**2-(2*x+1)/(z*(x-z)*(x+1)**2))/unitless_e_beta(x)**2
end function unitless_Moller

real function unitless_Moller_F(x,z)
    ! This function is the antiderivative of unitless_Moller(x,z) w.r.t. z. It is meaningless unless
    ! a difference is taken in z.
    implicit none
    real, intent(in) :: x ! incident electron energy
    real, intent(in) :: z ! outgoing electron energy

    unitless_Moller_F = &
        ((2*x+1)*log((x-z)/z)/(x*(x+1)**2) + z/(x+1)**2-1/z+1/(x-z))/unitless_e_beta(x)**2
end function unitless_Moller_F

real function unitless_Moller_R_Sigmat(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! Requires specification of:
    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g'+2
    ! XSw: Shell binding energy (but, always zero for Moller)

    unitless_Moller_R_Sigmat = &
        merge(0.0,unitless_Moller_F(x,inel_zmax(x)) - unitless_Moller_F(x,inel_zmin(x)),&
            inel_zmax(x) .eq. inel_zmin(x))
end function unitless_Moller_R_Sigmat

real function unitless_Moller_G(x,z)
    ! This function is the antiderivative of z*unitless_Moller(x,z) w.r.t. z. It is meaningless unless
    ! a difference is taken in z.
    implicit none
    real, intent(in) :: x ! incident electron energy
    real, intent(in) :: z ! outgoing electron energy

    unitless_Moller_G = &
        (((2*x+1)/(x+1)**2+1)*log(x-z) + log(z) + x/(x-z) + z**2/(2*(x+1)**2))/unitless_e_beta(x)**2
end function unitless_Moller_G

real function unitless_Moller_RLET(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! XSw: Must be zero
    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g'+2
    ! XSI: Mean excitation energy

    unitless_Moller_RLET = merge(&
        merge(unitless_Bethe(x), &
        unitless_Bethe(x) - (unitless_Moller_G(x,inel_zmax(x)) - unitless_Moller_G(x,inel_zmin(x))), &
        x .lt. XScut), &
        merge(unitless_Bethe(x) ,unitless_Bethe(x) - (unitless_Moller_G(x,inel_zmax(x)) - unitless_Moller_G(x,inel_zmin(x))), &
        XScut .eq. 0), &
        XSlogic)
end function unitless_Moller_RLET

real function unitless_Gryzinski(x)
    implicit none
    real, intent(in) :: x ! incident electron energy in units of electron rest masses

    real :: beta2
    real :: betaW2
    real :: f

    ! XSw: Shell binding energy in units of electron rest masses

    beta2 = unitless_e_beta(x)**2
    betaW2 = unitless_e_beta(XSw)**2

    f = betaW2*(beta2/(beta2+betaW2-beta2*betaW2))**(3.0/2)/beta2

    unitless_Gryzinski = merge(0.5*f*(1-XSw/x)**(3.0/2)*&
        (1+2*(1-0.5*XSw/x)*log(2.7+sqrt(x/XSw-1))/3)/XSw**2, 0.0, x .gt. XSw)
end function unitless_Gryzinski

real function MOLLER_TEST_PRIMARY(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_gamma
    real :: l_delta
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = XSd
    l_gamma = min(max(XSc,x/2),l_delta)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MOLLER_TEST_PRIMARY = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*(inel_mu(x,y)**XSlp)*&
        unitless_Moller(x,y)
end function MOLLER_TEST_PRIMARY

real function MOLLER_TEST_SECONDARY(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_gamma
    real :: l_delta
    real :: z

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = min(XSd,x/2)
    l_gamma = min(XSc,l_delta)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MOLLER_TEST_SECONDARY = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*(inel_mu(x,z)**XSlp)*&
        unitless_Moller(x,z)
end function MOLLER_TEST_SECONDARY

real function MOLLER_TEST_TOTAL(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_gamma
    real :: l_delta
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = XScut
    l_gamma = min(XScut,x/2)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MOLLER_TEST_TOTAL = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*&
        unitless_Moller(x,y)
end function MOLLER_TEST_TOTAL

real function MGXS_electron_electron_Moller_primary_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    ! XSw: Must be zero
    ! XSa: x_g'+1
    ! XSb: x_g'
    ! XSc: x_g+1
    ! XSd: x_g
    ! XSlp: Power of mu
    ! XScut: x_g'+2

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_ymax(x)),XSc)
    l_gamma = min(max(XSc,inel_ymin(x)),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_electron_electron_Moller_primary_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,y)**XSlp)*unitless_Moller(x,y)
end function MGXS_electron_electron_Moller_primary_integrand

real function MGXS_electron_electron_Moller_secondary_integrand(u,v) ! Use merge with XSlogic2 to use either zmin or ymin, i.e., have primary and secondary in same block?? If using XSlogic isn't slow, then I will.
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: z

    ! XSw: Must be zero
    ! XSa: x_g'+1
    ! XSb: x_g'
    ! XSc: x_g+1
    ! XSd: x_g
    ! XSlp: Power of mu
    ! XScut: x_g'+2

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_zmax(x)),XSc)
    l_gamma = min(max(XSc,inel_zmin(x)),XSd)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_electron_electron_Moller_secondary_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,z)**XSlp)*unitless_Moller(x,z)
end function MGXS_electron_electron_Moller_secondary_integrand

real function MGXS_electron_Moller_R_total_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: l_delta
    real :: l_gamma
    real :: x

    ! XSw: Must be zero
    ! XSa: x_g+1
    ! XSb: x_g
    ! XScut: x_g'+2

    l_delta = max(XSb,inel_zmin(x))
    l_gamma = max(XSa,inel_zmin(x))
    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)

    MGXS_electron_Moller_R_total_integrand = &
        0.5*(l_delta-l_gamma)*unitless_Moller_R_Sigmat(x)
end function MGXS_electron_Moller_R_total_integrand

real function FEXS_electron_electron_Moller_primary_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    ! XSw: Must be zero
    ! XSa: x^g'_2
    ! XSb: x^g'_1
    ! XSc: x^g_2
    ! XSd: x^g_1
    ! XSEgp: x^g'_3-n'
    ! XSEg: x^g_3-n
    ! XSlp: Power of mu
    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g'+2

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_ymax(x)),XSc)
    l_gamma = min(max(XSc,inel_ymin(x)),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    FEXS_electron_electron_Moller_primary_integrand = &
        0.25*(XSEgp-x)*(XSEg-y)*&
        (XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,y)**XSlp)*unitless_Moller(x,y)
end function FEXS_electron_electron_Moller_primary_integrand

real function FEXS_electron_electron_Moller_secondary_integrand(u,v) ! Use merge with XSlogic2 to use either zmin or ymin, i.e., have primary and secondary in same block?? If using XSlogic isn't slow, then I will.
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: z

    ! XSw: Must be zero
    ! XSa: x^g'_2
    ! XSb: x^g'_1
    ! XSc: x^g_2
    ! XSd: x^g_1
    ! XSEgp: x^g'_3-n'
    ! XSEg: x^g_3-n
    ! XSlp: Power of mu
    ! XScut: RCSDA cutoff energy

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_zmax(x)),XSc)
    l_gamma = min(max(XSc,inel_zmin(x)),XSd)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    FEXS_electron_electron_Moller_secondary_integrand = &
        0.25*(XSEgp-x)*(XSEg-z)*&
        (XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,z)**XSlp)*unitless_Moller(x,z)
end function FEXS_electron_electron_Moller_secondary_integrand

real function FEXS_electron_Moller_R_total_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: l_delta
    real :: l_gamma
    real :: x

    ! XSw: Must be zero
    ! XSa: x^g_2
    ! XSb: x^g_1
    ! XSEgp: x^g_3-n'
    ! XSEg: x^g_3-n
    ! XScut: RCSDA cutoff energy

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)

    FEXS_electron_Moller_R_total_integrand = &
        0.5*(XSb-XSa)*(XSEgp-x)*(XSEg-x)*unitless_Moller_R_Sigmat(x)
end function FEXS_electron_Moller_R_total_integrand

real function FEXS_electron_electron_Moller_RCSDA_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: x

    ! XSw: Must be zero
    ! XSI: Mean excitation energy
    ! XScut: RCSDA cutoff energy
    ! XSa: x^g_2
    ! XSb: x^g_1
    ! XSEgp: x^g_3-n'
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)

    FEXS_electron_electron_Moller_RCSDA_integrand = &
        0.5*(XSEgp-x)*(XSb-XSa)*unitless_Moller_RLET(x)
end function FEXS_electron_electron_Moller_RCSDA_integrand

real function RCSDA_MMS_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: x

    ! XSw: Must be zero
    ! XSI: Mean excitation energy
    ! XScut: RCSDA cutoff energy
    ! XSa: x^g_2
    ! XSb: x^g_1
    ! XSEgp: x^g_3-n'
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)

    RCSDA_MMS_integrand = &
        0.5*(XSEgp-x)*(XSb-XSa)*RCSDA_MMS_LET(x)
end function RCSDA_MMS_integrand

    real function RCSDA_MMS_LET(x)
        implicit none
        real, intent(in) :: x

        RCSDA_MMS_LET = sqrt(x)
    end function RCSDA_MMS_LET

    real function RCSDA_MMS_fluence(x)
        implicit none
        real, intent(in) :: x

        real :: a = 0.1
        real :: b = 1
        real :: c = -3

        RCSDA_MMS_fluence = (8*a+4*b+2*c)-(a*x**3+b*x**2+c*x)
    end function RCSDA_MMS_fluence

    real function RCSDA_MMS_source(x)
        implicit none
        real, intent(in) :: x

        real :: a = 0.1
        real :: b = 1
        real :: c = -3

        RCSDA_MMS_source = -(7*a*x**3+5*b*x**2+3*c*x-2*c-4*b-8*a)/(2*sqrt(x))
    end function RCSDA_MMS_source

    real function RCSDA_MMS_grouped_fluence(AA,BB)
        implicit none
        real, intent(in) :: AA
        real, intent(in) :: BB

        real :: a = 0.1
        real :: b = 1
        real :: c = -3

        RCSDA_MMS_grouped_fluence = -(BB-AA)*((6*BB+6*AA-24)*c+(4*BB**2+4*AA*BB+4*AA**2-48)*b + &
            (3*BB**3+3*AA*BB**2+3*AA**2*BB+3*AA**3-96)*a)/12
    end function RCSDA_MMS_grouped_fluence

    real function RCSDA_MMS_grouped_source(AA,BB)
        implicit none
        real, intent(in) :: AA
        real, intent(in) :: BB

        real :: a = 0.1
        real :: b = 1
        real :: c = -3

        RCSDA_MMS_grouped_source = -sqrt(BB)*((2*BB-4)*c+(2*BB**2-8)*b+(2*BB**3-16)*a)/2 - &
            sqrt(AA)*((4-2*AA)*c+(8-2*AA**2)*b+(16-2*AA**3)*a)/2
    end function RCSDA_MMS_grouped_source

real function unitless_RBED(x,z)
    implicit none
    real, intent(in) :: x ! incident electron energy
    real, intent(in) :: z ! outgoing electron energy

    real :: l_beta2

    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy

    l_beta2 = unitless_e_beta(x)**2

    unitless_RBED = &
        (-4*(2*x+1)/((x+2)**2*(x-z)*(z+XSw)) &
        + 1/(x-z)**2 + 1/(z+XSw)**2 + 4/(x+2)**2 &
        + XSw*(1/(z+XSw)**3+1/(x-z)**3)*(log(l_beta2/(1-l_beta2))-l_beta2-log(2*XSw))) &
        /(l_beta2 + XSbeta_U2 + XSbeta_W2)
end function unitless_RBED

real function RBED_TEST(v)
    implicit none
    real, intent(in) :: v

    real :: z

    z = 0.5*(inel_zmax(XSx)-inel_zmin(XSx))*v + 0.5*(inel_zmax(XSx)+inel_zmin(XSx))

    RBED_TEST = 0.5*(inel_zmax(XSx)-inel_zmin(XSx))*unitless_RBED(XSx,z)
end function RBED_TEST

real function unitless_RBEB(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    real :: l_beta2

    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy

    l_beta2 = unitless_e_beta(x)**2

    unitless_RBEB = &
        merge(0.0,(0.5*(1-XSw**2/x**2)*(log(l_beta2/(1-l_beta2)) - l_beta2 - log(2*XSw)) + &
        1 - XSw/x - log(x/XSw)*(1+2*x)/((x/XSw + 1)*(1+x/2)**2) + XSw**2*(x/XSw-1)/(2*(1+x/2)**2))/&
        ((l_beta2 + XSbeta_U2 + XSbeta_W2)*XSw), x .lt. XSw)
end function unitless_RBEB

real function unitless_RBEBav(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    real :: l_beta2

    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy

    l_beta2 = unitless_e_beta(x)**2

    unitless_RBEBav = 0.5*(1+(l_beta2 + XSbeta_U2 + XSbeta_W2)/l_beta2)*unitless_RBEB(x)
end function unitless_RBEBav

real function unitless_RBED_F(x,z)
    ! This function is the antiderivative of unitless_RBED(x,z) w.r.t. z. It is meaningless unless
    ! a difference is taken in z.
    implicit none
    real, intent(in) :: x ! incident electron energy
    real, intent(in) :: z ! outgoing electron energy

    real :: l_beta2

    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy

    l_beta2 = unitless_e_beta(x)**2

    unitless_RBED_F = &
        (-4*(2*x+1)*log((z+XSw)/(x-z))/((x+XSw)*(x+2)**2) + 1/(x-z) - 1/(z+XSw) + 4*z/(x+2)**2 - &
        0.5*XSw*(1/(z+XSw)**2-1/(x-z)**2)*(log(l_beta2/(1-l_beta2))-l_beta2-log(2*XSw)))/&
        (l_beta2 + XSbeta_U2 + XSbeta_W2)
end function unitless_RBED_F

real function unitless_RBED_G(x,z)
    ! This function is the antiderivative of (z+w)*unitless_RBED(x,z) w.r.t. z. It is meaningless unless
    ! a difference is taken in z.
    implicit none
    real, intent(in) :: x ! incident electron energy
    real, intent(in) :: z ! outgoing electron energy

    real :: l_beta2

    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy

    l_beta2 = unitless_e_beta(x)**2

    unitless_RBED_G = &
        (4*(2*x+1)*(x*log(x-z)+XSw*log(z+XSw))/((x+XSw)*(x+2)**2) + log((x-z)*(z+XSw)) - &
        x/(z-x) + XSw/(z+XSw) + 2*z**2/(x+2)**2 - &
        0.5*XSw*((2*z+XSw)/(z+XSw)**2 - (2*z-x)/(x-z)**2)*&
        (log(l_beta2/(1-l_beta2))-l_beta2-log(2*XSw)))/(l_beta2 + XSbeta_U2 + XSbeta_W2) + &
        XSw*unitless_RBED_F(x,z)
end function unitless_RBED_G

real function unitless_RBED_Sigmat(x)
    ! This function may in some cases reduce to RBEB (case is XScut = 0).
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy
    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g'+2

    unitless_RBED_Sigmat = &
        merge(0.0,unitless_RBED_F(x,inel_zmax(x)) - unitless_RBED_F(x,inel_zmin(x)),&
            inel_zmax(x) .eq. inel_zmin(x))
end function unitless_RBED_Sigmat

real function unitless_RBED_LET(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! Must have XScut = 0
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy

    unitless_RBED_LET = &
        merge(0.0,unitless_RBED_G(x,inel_zmax(x)) - unitless_RBED_G(x,inel_zmin(x)),&
            inel_zmax(x) .eq. inel_zmin(x))
end function unitless_RBED_LET

real function unitless_RBED_RLET(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    ! XScut: FEXS: RCSDA cutoff energy. MGXS: x_g+2
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)
    ! XSw: Shell binding energy

    if (XSlogic) then
        ! RLET is zero when x - w - 2*d < 0
        unitless_RBED_RLET = &
            merge(unitless_RBED_G(x,XScut) - unitless_RBED_G(x,0.0),&
                0.0, x-XSw-2*XScut .gt. 0)
    else
        ! RLET is zero when -x - w + 2*x_g+2 < 0
        unitless_RBED_RLET = &
            unitless_RBED_G(x,inel_zmin(x)) - unitless_RBED_G(x,0.0)
        !    max(unitless_RBED_G(x,min(x-XScut,0.5*(x-XSw))) - unitless_RBED_G(x,0.0), 0.0)
            !merge(unitless_RBED_G(x,x-XScut) - unitless_RBED_G(x,0.0),&
            !    0.0, .true.)!-x-XSw+2*XScut .gt. 0)
        !if (unitless_RBED_RLET .eq. 0) then
        !    print *, x, XSw, XScut
        !end if
    end if
end function unitless_RBED_RLET

real function MGXS_electron_electron_RBED_primary_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    ! XSa: x_g'+1
    ! XSb: x_g'
    ! XSc: x_g+1
    ! XSd: x_g
    ! XSlp: Power of mu
    ! XScut: x_g'+2

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_ymax(x)),XSc)
    l_gamma = min(max(XSc,inel_ymin(x)),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_electron_electron_RBED_primary_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,y)**XSlp)*unitless_RBED(x,y)
end function MGXS_electron_electron_RBED_primary_integrand

real function MGXS_electron_electron_RBED_secondary_integrand(u,v) ! Use merge with XSlogic2 to use either zmin or ymin, i.e., have primary and secondary in same block?? If using XSlogic isn't slow, then I will.
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: z

    ! XSa: x_g'+1
    ! XSb: x_g'
    ! XSc: x_g+1
    ! XSd: x_g
    ! XSlp: Power of mu
    ! XScut: x_g'+2
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_zmax(x)),XSc)
    l_gamma = min(max(XSc,inel_zmin(x)),XSd)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_electron_electron_RBED_secondary_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,z)**XSlp)*unitless_RBED(x,z)
end function MGXS_electron_electron_RBED_secondary_integrand

real function MGXS_electron_RBEB_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: l_delta
    real :: l_gamma
    real :: x

    ! XScut must be 0
    ! XSa: x_g+1
    ! XSb: x_g
    ! XSw: Shell binding energy
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    l_delta = max(XSb,XSw)
    l_gamma = max(XSa,XSw)
    x = 0.5*(l_delta-l_gamma)*u + 0.5*(l_delta+l_gamma)

    MGXS_electron_RBEB_integrand = &
        0.5*(l_delta-l_gamma)*unitless_RBEB(x)
end function MGXS_electron_RBEB_integrand

real function MGXS_electron_RBED_total_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: l_delta
    real :: l_gamma
    real :: x

    ! XSa: x_g+1
    ! XSb: x_g
    ! XScut: x_g+2
    ! XSw: Shell binding energy
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    l_delta = max(XSb,XSw) !!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTE: So, I removed zmin(x). Had to because x defined after. But is this valid? I am basically assuming that I will not feed an x such that x < x_g+2, right? If this is my assumption then I am fine. If not, I may not be fine.
    l_gamma = max(XSa,XSw) !!!!!!!!!!!!!!!!!!!!!!!!!!!! NOTE:
    x = 0.5*(l_delta-l_gamma)*u + 0.5*(l_delta+l_gamma)

    MGXS_electron_RBED_total_integrand = &
        0.5*(l_delta-l_gamma)*unitless_RBED_Sigmat(x)
end function MGXS_electron_RBED_total_integrand

real function FEXS_electron_electron_RBED_primary_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    ! XSa: x^g'_2
    ! XSb: x^g'_1
    ! XSc: x^g_2
    ! XSd: x^g_1
    ! XSEgp: x^g'_3-n'
    ! XSEg: x^g_3-n
    ! XSlp: Power of mu
    ! XScut: RCSDA cutoff energy
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_ymax(x)),XSc)
    l_gamma = min(max(XSc,inel_ymin(x)),XSd)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    FEXS_electron_electron_RBED_primary_integrand = &
        0.25*(XSEgp-x)*(XSEg-y)*&
        (XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,y)**XSlp)*unitless_RBED(x,y)
end function FEXS_electron_electron_RBED_primary_integrand

real function FEXS_electron_electron_RBED_secondary_integrand(u,v) ! Use merge with XSlogic2 to use either zmin or ymin, i.e., have primary and secondary in same block?? If using XSlogic isn't slow, then I will.
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: z

    ! XSa: x_g'+1
    ! XSb: x_g'
    ! XSc: x_g+1
    ! XSd: x_g
    ! XSlp: Power of mu
    ! XScut: RCSDA cutoff energy
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,inel_zmax(x)),XSc)
    l_gamma = min(max(XSc,inel_zmin(x)),XSd)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    FEXS_electron_electron_RBED_secondary_integrand = &
        0.25*(XSEgp-x)*(XSEg-z)*&
        (XSb-XSa)*(l_delta-l_gamma)*&
        (inel_mu(x,z)**XSlp)*unitless_RBED(x,z)
end function FEXS_electron_electron_RBED_secondary_integrand

real function FEXS_electron_RBEB_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: l_delta
    real :: l_gamma
    real :: x

    ! XScut must be 0
    ! XSa: x^g_2
    ! XSb: x^g_1
    ! XSEgp: x^g_3-n'
    ! XSEg: x^g_3-n
    ! XSw: Shell binding energy
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    l_delta = max(XSb,XSw)
    l_gamma = max(XSa,XSw)
    x = 0.5*(l_delta-l_gamma)*u + 0.5*(l_delta+l_gamma)

    FEXS_electron_RBEB_integrand = &
        0.5*(XSEgp-x)*(XSEg-x)*(l_delta-l_gamma)*unitless_RBEB(x)
end function FEXS_electron_RBEB_integrand

real function FEXS_electron_RBED_total_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: l_delta
    real :: l_gamma
    real :: x

    ! XScut: RCSDA cutoff energy
    ! XSa: x^g_2
    ! XSb: x^g_1
    ! XSEgp: x^g_3-n'
    ! XSEg: x^g_3-n
    ! XSw: Shell binding energy
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    l_delta = max(XSb, 2*XScut + XSw)
    l_gamma = max(XSa, 2*XScut + XSw)
    x = 0.5*(l_delta-l_gamma)*u + 0.5*(l_delta+l_gamma)

    FEXS_electron_RBED_total_integrand = &
        0.5*(XSEgp-x)*(XSEg-x)*(l_delta-l_gamma)*unitless_RBED_Sigmat(x)
end function FEXS_electron_RBED_total_integrand

real function FEXS_electron_electron_RBED_RCSDA_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: l_delta
    real :: l_gamma
    real :: x

    ! XScut: RCSDA cutoff energy
    ! XSa: x^g_2
    ! XSb: x^g_1
    ! XSEgp: x^g_3-n'
    ! XSw: Shell binding energy
    ! XSbeta_W2: Electron speed for energy W (shell binding energy)
    ! XSbeta_U2: Electron speed for energy U (shell kinetic energy)

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)

    FEXS_electron_electron_RBED_RCSDA_integrand = &
        0.5*(XSEgp-x)*(XSb-XSa)*unitless_RBED_RLET(x)
end function FEXS_electron_electron_RBED_RCSDA_integrand

! ----- Elastic Electron
real function MGXS_electron_electron_elastic_fitted_integrand(b) ! CAN OPTIMIZE?
    implicit none
    real, intent(in) :: b

    real :: A
    real :: C

    ! XSZ is Z
    ! XSmu is mu
    ! XSk is k

    A = (fsa*XSZ**(1.0/3)/1.77608)**2*(1-b**2)*(1.13+3.76*(fsa*XSZ/b)**2)/b**2
    C = 1.0/(2*A+1-XSmu)**2

    MGXS_electron_electron_elastic_fitted_integrand = &
    (b-0.7181287)**(XSk-1)*C/sqrt(b**6 - b**8)
end function MGXS_electron_electron_elastic_fitted_integrand

real function Wentzel_Moliere_screening_parameter(x)
    implicit none
    real, intent(in) :: x

    ! XSZ: atomic number

    Wentzel_Moliere_screening_parameter = 0.25*((XSZ**(1.0/3)*fsa/0.885)**2)*&
        (1.13+3.76*(XSZ**2)*(fsa**2)*(x+1)**2/(x*(x+2)))/(x*(x+2))
end function Wentzel_Moliere_screening_parameter

real function Wentzel_Moliere_no_angular(x)
    implicit none
    real, intent(in) :: x

    ! XSZ: atomic number

    Wentzel_Moliere_no_angular = XSZ*(XSZ+1)*(x+1)**2/((x**2)*((x+2)**2))
end function Wentzel_Moliere_no_angular

real function Wentzel_Moliere_Legendre_angular(A)
! Moments of 1/(A-mu)^2. NOT Legendre moments. Don't know why I named it that
    implicit none
    real, intent(in) :: A

    integer :: k

    ! XSlp: power of moment

    Wentzel_Moliere_Legendre_angular = 0
    do k = 0, XSlp
        Wentzel_Moliere_Legendre_angular = Wentzel_Moliere_Legendre_angular + &
            ((-1)**k)*binomial_int(XSlp,k)*(A**(XSlp-k))*&
            merge(log((A+1)/(A-1)), (A+1)**k/((A+1)*k-A-1) - (A-1)**k/((A-1)*k-A+1), k .eq. 1)
    end do
end function Wentzel_Moliere_Legendre_angular

real function Wentzel_Moliere_total_attn(x)
    implicit none
    real, intent(in) :: x

    real :: eta

    ! XSZ: atomic number

    eta = Wentzel_Moliere_screening_parameter(x)

    Wentzel_Moliere_total_attn = pi*XSZ*(XSZ+1)*(x+1)**2/((x**2)*((x+2)**2)*(eta*(eta+1)))
end function Wentzel_Moliere_total_attn

real function temp_Wentzel_Moliere(mu)
    implicit none
    real, intent(in) :: mu

    real :: eta

    !eta = 8.176E-4
    eta = 5.699E-5

    temp_Wentzel_Moliere = 2*eta*(1+eta)/(1+2*eta-mu)**2
end function temp_Wentzel_Moliere

real function temp_Wentzel_Moliere_moments(mu)
    implicit none
    real, intent(in) :: mu

    real :: eta

    !eta = 8.176E-4
    eta = 5.699E-5

    temp_Wentzel_Moliere_moments = (mu**XSlp)*2*eta*(1+eta)/(1+2*eta-mu)**2
end function temp_Wentzel_Moliere_moments

real function MGXS_electron_electron_elastic_Wentzel_Moliere_integrand(u)
    implicit none
    real, intent(in) :: u

    real :: x
    real :: eta
    real :: A
    real :: scrfa

    ! XSZ: atomic number
    ! XSb: upper bound of energy group
    ! XSa: lower bound
    ! XSlp: power of moment

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    eta = Wentzel_Moliere_screening_parameter(x)
    A = 1 + 2*eta
    scrfa = Wentzel_Moliere_Legendre_angular(A)

    MGXS_electron_electron_elastic_Wentzel_Moliere_integrand = &
        0.5*(XSb-XSa)*Wentzel_Moliere_no_angular(x)*scrfa
end function MGXS_electron_electron_elastic_Wentzel_Moliere_integrand

! Bremsstrahlung
real function rad_zmin(x)
    implicit none
    real, intent(in) :: x

    ! XScut: MGXS: x_g+2, FEXS: delta

    rad_zmin = merge(XScut, x - XScut, XSlogic)
end function rad_zmin

real function rad_zmax(x)
    implicit none
    real, intent(in) :: x

    ! XScut: MGXS: x_g+2, FEXS: delta

    rad_zmax = max(x - XScut, XScut) ! Actually the same expression for MGXS and FEXS.
end function rad_zmax

real function Bethe_Heitler_DBREM(x,z)
    implicit none
    real, intent(in) :: x ! Incoming electron energy in electron rest energies
    real, intent(in) :: z ! Outgoing electron energy

    real :: g
    real :: e
    real :: b
    real :: p1
    real :: p2

    ! XSR is reduced screening radius

    g = 1 + x
    e = z/g
    b = XSR*e/(2*g*(1-e))
    p1 = 4*log(XSR) + 2 - 2*log(1+b**2) - 4*b*atan(1/b)
    p2 = 4*log(XSR) + 7.0/3 - 2*log(1+b**2) - 6*b*atan(1/b) - &
        b**2*(4 - 4*b*atan(1/b) - 3*log(1+1/b**2))

    Bethe_Heitler_DBREM = (e**2*p1 + 4*(1-e)*p2/3)/z
end function Bethe_Heitler_DBREM

real function MGXS_electron_electron_Bremsstrahlung_integrand(u,v)
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_gamma
    real :: l_delta
    real :: z

    ! XSR is reduced screening radius
    ! XScut: MGXS: x_g+2, FEXS: d
    ! XSa, XSb, XSc, XSd

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = max(min(XSd,rad_zmax(x)),XSc)
    l_gamma = min(max(XSc,rad_zmin(x)),XSd)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_electron_electron_Bremsstrahlung_integrand = &
        0.25*(l_delta-l_gamma)*(XSb-XSa)*Bethe_Heitler_DBREM(x,z)
end function MGXS_electron_electron_Bremsstrahlung_integrand

real function MGXS_electron_Bremsstrahlung_total_integrand(u,v)
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_gamma
    real :: l_delta
    real :: z

    ! XSR is reduced screening radius
    ! XScut: MGXS: x_g+2, FEXS: d
    ! XSa, XSb

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = rad_zmax(x)
    l_gamma = rad_zmin(x)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_electron_Bremsstrahlung_total_integrand = &
        0.25*(l_delta-l_gamma)*(XSb-XSa)*Bethe_Heitler_DBREM(x,z)
end function MGXS_electron_Bremsstrahlung_total_integrand

! ----- Stopping Power / LET
real function unitless_Bethe(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    real :: beta2
    real :: B0
    real :: x0
    real :: dB0dx
    real :: g

    ! XSI: Mean excitation energy


    !if (x .gt. 0.01957) then
        beta2 = unitless_e_beta(x)**2

        unitless_Bethe = 2*(log(2*beta2/(XSI*(1-beta2)))/beta2 - 1)
    !else
    !    x0 = 0.01956951183
    !    beta2 = unitless_e_beta(x0)**2
    !    B0 = 2*(log(2*beta2/(XSI*(1-beta2)))/beta2 - 1)
    !    dB0dx = -(x0+1)**3*(B0 - 2/(beta2*(1-beta2)) + 2)/(2*beta2)
    !    g = -x0*dB0dx/B0
    !
    !    unitless_Bethe = B0*(x0/x)**g
    !end if
    !
    !unitless_Bethe = unitless_Bethe !- unitless_density_effect(x)
end function unitless_Bethe

real function TEST_CSDA_BETA(x)
    implicit none
    real, intent(in) :: x ! incident electron energy

    real :: beta2

    ! XSI: Mean excitation energy

    beta2 = unitless_e_beta(x)**2

    TEST_CSDA_BETA = (log(x**2*(x+2)/(2*XSI**2))+1-beta2-(2*x+1)*log(2.0)/(x+1)**2 + x**2/(8*(x+1)**2))
end function TEST_CSDA_BETA

real function TEST_Moller_primary_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y


    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_gamma = max(XSc,x/2)
    l_delta = max(XSd,l_gamma)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    TEST_Moller_primary_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*(inel_mu(x,y)**XSlp)*unitless_Moller(x,y)
end function TEST_Moller_primary_integrand

real function TEST_Moller_secondary_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: z


    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_gamma = XSc
    l_delta = max(min(XSd,x/2),l_gamma)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    TEST_Moller_secondary_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*(inel_mu(x,z)**XSlp)*unitless_Moller(x,z)
end function TEST_Moller_secondary_integrand

real function TEST_Moller_total_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: z

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = XScut
    l_gamma = min(x/2,XScut)
    z = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    TEST_Moller_total_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*unitless_Moller(x,z)
end function TEST_Moller_total_integrand

real function TEST_Moller_FO_RCSDA_integrand(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: l_delta
    real :: l_gamma
    real :: y

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_delta = XScut
    l_gamma = min(x/2,XScut)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    TEST_Moller_FO_RCSDA_integrand = &
        0.25*(XSb-XSa)*(l_delta-l_gamma)*(x-y)*unitless_Moller(x,y)
end function TEST_Moller_FO_RCSDA_integrand

real function unitless_density_effect(x)
    implicit none
    real, intent(in) :: x

    real :: xp

    ! XSBP: B parameter
    ! XSCP: C parameter
    ! XSx0: x0 parameter
    ! XSx1: x1 parameter

    xp = log10(sqrt(2*x+x**2))

    if (xp .le. XSx0) then
        unitless_density_effect = 0
    else if (xp .gt. XSx0 .and. xp .le. XSx1) then
        unitless_density_effect = 4.606*xp + XSCP + XSBP*(XSx1-xp)**3
    else
        unitless_density_effect = 4.606*xp + XSCP
    end if
end function unitless_density_effect

subroutine mean_excitation_energy(elements, Numatoms, Ibar)
    implicit none
    integer, dimension(:), intent(in) :: elements
    integer, dimension(:), intent(in) :: Numatoms

    real, intent(inout) :: Ibar ! eV

    integer :: e
    integer :: Z
    real :: l_I
    real :: normterm
    real :: sumterm

    normterm = 0
    sumterm = 0

    do e = 1, size(elements)
        Z = elements(e)
        if (Z .eq. 0) exit
        if (Z .gt. 12) then
            l_I = Z*(9.76+58.8*Z**(-1.19))
        else if (Z .eq. 1) then
            l_I = 19.2
        else if (Z .eq. 8) then
            l_I = 95.0
        else
            print *, "WIP: Need Wbar values for all Z .le. 12."
            stop
        end if

        normterm = normterm + Z*Numatoms(e)
        sumterm = sumterm + Z*Numatoms(e)*log(l_I)
    end do

    Ibar = exp(sumterm/normterm)
    !Ibar = 75.0
end subroutine mean_excitation_energy

subroutine set_density_effect_parameters(rho, GAMmat, elements, Numatoms, Ibar, is_it_gas)
    implicit none
    real, intent(in) :: rho ! Material's density
    real, intent(in) :: GAMmat ! Material's GAM
    integer, dimension(:), intent(in) :: elements
    integer, dimension(:), intent(in) :: Numatoms
    real, intent(in) :: Ibar ! eV, right??
    logical, intent(in) :: is_it_gas

    real :: plasma_freq

    plasma_freq = 28.8*sqrt(rho*sum(elements*Numatoms)/GAMmat)

    XSCP = -2*log(Ibar/plasma_freq)-1

    if (is_it_gas) then
        print *, "GAAAAAASSSSSSSSSSSSSSS"
        stop
    else
        if (Ibar .ge. 100 .and. -XSCP .ge. 5.215) then
            XSx0 = -0.326*XSCP-1.5
        else if (Ibar .lt. 100 .and. -XSCP .ge. 3.681) then
            XSx0 = -0.326*XSCP-1.0
        else
            XSx0 = 2
        end if

        if (Ibar .ge. 100) then
            XSx1 = 3
        else
            XSx1 = 2
        end if
    end if

    XSBP = -(XSCP + 4.606*XSx0)/(XSx1 - XSx0)**3
end subroutine set_density_effect_parameters

! ----- Photoelectric Effect
real function Sauter(mu,y)
    implicit none
    real, intent(in) :: mu
    real, intent(in) :: y

    real :: g
    real :: nu
    real :: N
    real :: A

    g = 1 + y
    nu = 1 - mu
    A = g/sqrt(g**2-1)-1
    N = 0.5*sqrt(g+1)*(g-2)*(g-1)**(1.5)*log((g-sqrt(g**2-1))/(g+sqrt(g**2-1))) + &
        (1.0/3)*(g+1)*(g-1)**2*(3*g**2-2*g+4)

    Sauter = (1/(A+nu) + 0.5*sqrt(g+1)*(g-2)*(g-1)**(1.5))*nu*(2-nu)/(N*(A+nu)**3)
end function Sauter

! ----- Pair Production
real function Bethe_Heitler_DED(x,z)
    implicit none
    real, intent(in) :: x ! Photon energy in electron rest energies
    real, intent(in) :: z ! Fraction of photon energy kept as kinetic energy by electron

    real :: a
    real :: b
    real :: fC
    real :: F0
    real :: g0
    real :: g1
    real :: g2
    real :: phi1
    real :: phi2

    ! XSR is reduced screening radius
    ! XSZ is atomic number

    a = fsa*XSZ
    b = XSR/(2*x*z*(1-z))
    fC = a**2*(1/(1+a**2) + 0.202059 - 0.03693*a**2 + 0.00835*a**4 &
       - 0.00201*a**6 + 0.00049*a**8 - 0.00012*a**10 + 0.00003*a**12)
    F0 = (-0.1774 - 12.10*a + 11.18*a**2)*sqrt(2/x) &
         + (8.523 + 73.26*a - 44.41*a**2)*(2/x) &
         - (13.52 + 121.1*a - 96.41*a**2)*(2/x)**(3.0/2) &
         + (8.946 + 62.05*a - 63.41*a**2)*(2/x)**2
    g0 = 4*(log(XSR) - fC) + F0
    g1 = 7.0/3 - 2*log(1+b**2) - 6*b*atan(1/b) - b**2*(4 - 4*b*atan(1/b) - 3*log(1+1/b**2))
    g2 = 11.0/6 - 2*log(1+b**2) - 3*b*atan(1/b) + 0.5*b**2*(4 - 4*b*atan(1/b) - 3*log(1+1/b**2))
    phi1 = max(g1 + g0, 0.0)
    phi2 = max(g2 + g0, 0.0)

    Bethe_Heitler_DED = 2*(0.5-z)**2*phi1 + phi2
end function Bethe_Heitler_DED

! ----- Reading Data
subroutine read_atomic_parameters(Z, g, W, K)
    implicit none
    integer, intent(in) :: Z

    real, dimension(:), allocatable, intent(inout) :: g
    real, dimension(:), allocatable, intent(inout) :: W ! MeV
    real, dimension(:), allocatable, intent(inout) :: K ! MeV

    integer :: s, cascade
    integer :: Nshells
    character(50) :: Zstr
    real, dimension(16) :: trash16
    real, dimension(2) :: trash2
    real, dimension(:), allocatable :: B_g

    write(Zstr,*) Z

    allocate(g(0))
    open(1, file = "Physics data/Atomic parameters/EADL.dat")
    !   Locate section corresponding to Z
    do
        read(1,*) trash16
        if (trash16(1) .eq. Z*1000) exit
        do
            read(1,*) trash2
            if (all(trash2 .eq. -1137.0)) exit
        end do
    end do

    !   Begin reading degeneracy, count number of shells
    do
        read(1,*) trash2
        if (all(trash2 .eq. -1137.0)) exit
        allocate(B_g(size(g)+1))
        B_g(1:size(g)) = g
        B_g(size(B_g)) = trash2(2)
        call move_alloc(B_g, g)
    end do
    Nshells = size(g)
    allocate(W(Nshells))
    allocate(K(Nshells))

    !   Begin reading binding energies
    read(1,*) trash16
    do s = 1, Nshells+1
        read(1,*) trash2
        if (all(trash2 .eq. -1137.0)) exit
        W(s) = trash2(2)
    end do

    !   Begin reading kinetic energies
    read(1,*) trash16
    do s = 1, Nshells+1
        read(1,*) trash2
        if (all(trash2 .eq. -1137.0)) exit
        K(s) = trash2(2)
    end do
    close(1)
end subroutine read_atomic_parameters

subroutine read_relaxation_data(relaxation_type, Z, Nshells, Eij, eta, cscpershell)
    implicit none
    character(*), intent(in) :: relaxation_type
    integer, intent(in) :: Z
    integer, intent(in) :: Nshells ! I can (and should) remove this as an input, because it can be obtained while reading for a given Z

    real, dimension(:), allocatable, intent(inout) :: Eij ! MeV
    real, dimension(:), allocatable, intent(inout) :: eta
    integer, dimension(:), allocatable, intent(inout) :: cscpershell

    integer :: shell, cascade
    integer :: Ncascades
    character(50) :: Zstr
    real, dimension(2) :: trash2
    integer :: sectionctr
    integer :: cscctr
    real, dimension(:,:), allocatable :: augdata
    real, dimension(:,:), allocatable :: B_augdata

    write(Zstr,*) Z

    sectionctr = 0
    allocate(cscpershell(Nshells))
    cscpershell = 0
    allocate(augdata(2,0))
    if (relaxation_type .eq. "auger") then
        open(1, file = "Physics data/Auger/"//trim(adjustl(Zstr))//".dat")
    else
        open(1, file = "Physics data/Fluorescence/"//trim(adjustl(Zstr))//".dat")
    end if
    do shell = 1, Nshells
        cscctr = 0
        read(1,*) trash2
        if (all(trash2 .eq. -1)) exit
        do
            read(1,*) trash2
            if (all(trash2 .eq. -1137)) exit
            allocate(B_augdata(2,size(augdata,2)+1))
            B_augdata(1:2,1:size(augdata,2)) = augdata
            B_augdata(1:2,size(B_augdata,2)) = trash2
            call move_alloc(B_augdata, augdata)
            cscctr = cscctr + 1
        end do
        cscpershell(shell) = cscctr
    end do
    close(1)

    Ncascades = sum(cscpershell)
    allocate(Eij(Ncascades))
    allocate(eta(Ncascades))
    Eij = augdata(2,1:Ncascades)
    eta = augdata(1,1:Ncascades)
end subroutine read_relaxation_data

subroutine read_photoelectric_effect(Z, Emin, Emax, peWvals, Espec, pexs)
    implicit none
    integer, intent(in) :: Z
    real, intent(in) :: Emin
    real, intent(in) :: Emax

    real, dimension(:), allocatable, intent(inout) :: peWvals ! MeV
    real, dimension(:), allocatable, intent(inout) :: Espec ! mc^2
    real, dimension(:,:), allocatable, intent(inout) :: pexs ! (E,shell), cm^2

    integer :: shell, i
    integer :: max_NEspec = 420
    integer :: Nshells
    integer :: NEspec
    real :: u = 1.001
    character(50) :: Zstr
    integer, dimension(2) :: reader2
    real, dimension(:,:), allocatable :: Sigma_array
    real, dimension(:,:), allocatable :: B_Sigma_array
    integer :: eof
    integer :: startind
    integer :: finind
    real, dimension(:), allocatable :: l_Espec

    write(Zstr, *) Z

    open(1, file = "Physics data/Photoelectric/Nshells.txt")
    do i = 1, Z
        read(1,*) Nshells
    end do
    close(1)

    open(1, file = "Physics data/Photoelectric/"//trim(adjustl(Zstr))//".txt")

    allocate(Sigma_array(Nshells+2,max_NEspec))
    Sigma_array = 0
    do i = 1, max_NEspec
        read(1,*, iostat = eof) Sigma_array(:,i)
        if (eof .lt. 0) exit
    end do
    close(1)

    NEspec = findloc(Sigma_array(1,:) .eq. 0, .true., dim = 1) - 1
    allocate(l_Espec(NEspec))
    l_Espec = Sigma_array(1,1:NEspec)*10**(-6.0)

    allocate(peWvals(Nshells))
    do shell = 1, Nshells
        peWvals(shell) = l_Espec(findloc(Sigma_array(shell+2,:) .gt. 0, .true., dim = 1))
    end do
    peWvals = peWvals*e_mass_E ! To convert it back to MeV

    startind = findloc(l_Espec .ge. Emin, .true., dim = 1)
    if (l_Espec(startind) .gt. Emin) startind = startind - 1
    finind = findloc(l_Espec .ge. Emax, .true., dim = 1)

    allocate(Espec, source = l_Espec(startind:finind))

    allocate(B_Sigma_array(startind:finind,Nshells+1))
    B_Sigma_array(startind:finind,1:Nshells) = transpose(Sigma_array(3:Nshells+2,startind:finind))
    B_Sigma_array(startind:finind,Nshells+1) = Sigma_array(2,startind:finind) ! TOTAL XS IS LAST ENTRY
    allocate(pexs, source = B_Sigma_array/u)
end subroutine read_photoelectric_effect

subroutine read_pair_production(Z, Emax, ENist, ppxs)
    implicit none
    integer, intent(in) :: Z
    real, intent(in) :: Emax

    real, dimension(:), allocatable, intent(inout) :: ENist ! In electron rest masses
    real, dimension(:), allocatable, intent(inout) :: ppxs

    integer :: i
    integer :: NISTmax = 200
    integer :: NNIST
    integer :: startind
    character(50) :: Zstr
    character(50), dimension(3) :: trash3
    real :: trash1
    real, dimension(:,:), allocatable :: reader
    integer :: eof

    allocate(reader(3,NISTmax))

    write(Zstr,*) Z

    open(1, file = "Physics data/Pair production/red_screening_rad.txt")
    do i = 1, XSZ - 1
        read(1,*) trash1
    end do
    read(1,*) XSR
    close(1)

    open(1, file = "Physics data/Pair production/eta_infty.txt")
    do i = 1, XSZ - 1
        read(1,*) trash1
    end do
    read(1,*) XSeta
    close(1)

    open(1, file = "Physics data/Pair production/"//trim(adjustl(Zstr))//".txt" )
    do i = 1, 2
        read(1,*) trash3
    end do

    reader = 0
    do i = 1, NISTmax
        read(1,*, iostat = eof) reader(:,i)
        if (eof .lt. 0) exit
    end do
    close(1)

    startind = findloc(reader(1,:), 1.022, dim = 1)
    if (startind .eq. 0) then
        print *, "MGXS pair production: NIST cross section appears not to have E' = 1.022. Program ending."
        stop
    end if
    NNIST = findloc(reader(1,:) .gt. Emax, .true., dim = 1) - &
            startind + 1

    allocate(ENist(NNIST))
    allocate(ppxs(NNIST))

    ENist = reader(1,startind:NNIST+startind-1)
    ppxs = reader(2,startind:NNIST+startind-1) + reader(3,startind:NNIST+startind-1)

    ENist = ENist/e_mass_E
end subroutine read_pair_production

! ----- Other
real function MGXS_K_MMS(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: y
    real :: l_gamma
    real :: l_delta

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_gamma = XSc
    l_delta = max(min(XSd,x),XSc)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    MGXS_K_MMS = 0.25*(XSb-XSa)*(l_delta-l_gamma)*max(sign(1.0,x-y),0.0)/(2*x+1)
end function MGXS_K_MMS

real function FEXS_K_MMS(u,v)
    implicit none
    real, intent(in) :: u
    real, intent(in) :: v

    real :: x
    real :: y
    real :: l_gamma
    real :: l_delta

    x = 0.5*(XSb-XSa)*u + 0.5*(XSb+XSa)
    l_gamma = XSc
    l_delta = max(min(XSd,x),XSc)
    y = 0.5*(l_delta-l_gamma)*v + 0.5*(l_delta+l_gamma)

    FEXS_K_MMS = 0.25*(XSEgp-x)*(XSEg-y)*(XSb-XSa)*&
    (l_delta-l_gamma)*max(sign(1.0,x-y),0.0)/(2*x+1)
end function FEXS_K_MMS

real function K_MMS_XS(x,y)
    implicit none
    real, intent(in) :: x ! INCOMING ENERGY
    real, intent(in) :: y ! OUTGOING ENERGY

    K_MMS_XS = max(sign(1.0,x-y),0.0)/(2*x+1)
end function K_MMS_XS

real function K_MMS_vector(x)
    implicit none
    real, intent(in) :: x

    K_MMS_vector = 0.5*max(1-abs(2*x-1),0.0)
end function K_MMS_vector

real function FEXS_K_MMS_kernel(x)
    implicit none
    real, intent(in) :: x

    FEXS_K_MMS_kernel = K_MMS_XS(x,XSa)*K_MMS_vector(x)
end function FEXS_K_MMS_kernel

real function MGXS_K_MMS_kernel(x,y)
    implicit none
    real, intent(in) :: x
    real, intent(in) :: y

    MGXS_K_MMS_kernel = K_MMS_XS(x,y)*K_MMS_vector(x)
end function MGXS_K_MMS_kernel

end module physics