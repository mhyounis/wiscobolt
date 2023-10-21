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

module MGXS
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

! WIP LIST:
! Bremsstrahlung
! Electron -> Photon fluorescence

contains
subroutine MGXS_read_custom_XS(particle, S, St, EDEP, CDEP)
    implicit none
    integer, intent(in) :: particle

    real, dimension(:,:,:,:), allocatable, intent(inout) :: S
    real, dimension(:,:), allocatable, intent(inout) :: St
    real, dimension(:,:), allocatable, intent(inout) :: EDEP
    real, dimension(:,:), allocatable, intent(inout) :: CDEP

    integer :: l
    character(100) :: fname

    if (particle .eq. 1) then
        allocate(S(Gp,Gp,NL+1,Nmats))
        allocate(St(Gp,Nmats))
        allocate(EDEP(Gp,Nmats))
        allocate(CDEP(Gp,Nmats))

        fname = "Physics data/User-defined XS/photon_scattering.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) S
        close(1)
        fname = "Physics data/User-defined XS/photon_total.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) St
        close(1)
        fname = "Physics data/User-defined XS/photon_energy_dep.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) EDEP
        close(1)
        fname = "Physics data/User-defined XS/photon_charge_dep.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) CDEP
        close(1)

        do l = 0, NL
            S(:,:,l+1,:) = (2*l+1)*S(:,:,l+1,:)/fourpi
        end do
    else if (particle .eq. 12) then
        allocate(S(Gp,Ge,NL+1,Nmats))

        fname = "Physics data/User-defined XS/photon_to_electron_scattering.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) S
        close(1)

        do l = 0, NL
            S(:,:,l+1,:) = (2*l+1)*S(:,:,l+1,:)/fourpi
        end do
    else if (particle .eq. 2) then
        if (inelastic_ETC .or. exact_RCSDA_angular) then
            allocate(S(Ge,Ge,NL+2,Nmats))
        else
            allocate(S(Ge,Ge,NL+1,Nmats))
        end if
        allocate(St(Ge,Nmats))
        allocate(EDEP(Ge,Nmats))
        allocate(CDEP(Ge,Nmats))

        fname = "Physics data/User-defined XS/electron_scattering.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) S
        close(1)
        fname = "Physics data/User-defined XS/electron_total.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) St
        close(1)
        fname = "Physics data/User-defined XS/electron_energy_dep.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) EDEP
        close(1)
        fname = "Physics data/User-defined XS/electron_charge_dep.dat"
        open(1, file = fname, form = "unformatted", action = "read")
        read(1) CDEP
        close(1)

        do l = 0, NL
            S(:,:,l+1,:) = (2*l+1)*S(:,:,l+1,:)/fourpi
        end do
    end if
end subroutine MGXS_read_custom_XS

! SCATTERING
!--------------------- Photon-Photon --------------------------------------------------------------
subroutine MGXS_photon_photon(Ep, dEp, adaptive, tol, Nint, Sigmapp, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigmapp
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: mat
    real, dimension(:,:,:), allocatable :: Sigmappc
    real, dimension(:,:,:), allocatable :: SigmappF

    allocate(Sigmapp(Gp,Gp,NL+1,Nmats))
    allocate(EDEP(Gp,Nmats))
    Sigmapp = 0
    EDEP = 0

    call MGXS_photon_photon_Compton(Ep, dEp, adaptive, tol, Nint, Sigmappc, EDEP)
    call MGXS_photon_photon_fluorescence(Ep, dEp, SigmappF, EDEP)

    Sigmapp(1:Gp,1:Gp,1,1:Nmats) = SigmappF

    do mat = 1, Nmats
        Sigmapp(1:Gp,1:Gp,1:NL+1,mat) = Sigmapp(1:Gp,1:Gp,1:NL+1,mat) + mat_rho_e(mat)*Sigmappc
    end do
end subroutine MGXS_photon_photon

subroutine MGXS_photon_photon_Compton(Ep, dEp, adaptive, tol, Nint, Sigmappc, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:), allocatable, intent(inout) :: Sigmappc
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: g, gpr, l, lp
    real, dimension(:,:), allocatable :: coeff
    real, dimension(:,:,:), allocatable :: integral
    real :: l_integral
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:), allocatable :: portion

    allocate(Sigmappc(Gp,Gp,NL+1))
    allocate(coeff(0:NL,0:NL))
    allocate(integral(Gp,Gp,NL+1))
    allocate(portion(0:NL,Gp))

    coeff = Legendre_polynomial_coeffs(NL)

    integral = 0

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_legendre_quadrature(Nint, xq1, wq1)
    end if

    do gpr = 1, Gp
        XSa = Ep(gpr+1)/e_mass_E
        XSb = Ep(gpr)/e_mass_E
        do g = gpr, Gp
            XSc = Ep(g+1)/e_mass_E
            XSd = Ep(g)/e_mass_E
            do lp = 0, NL
                XSlp = lp
                if (adaptive) then
                    call Gauss_Legendre_adaptive_quadrature_2D &
                    (MGXS_photon_photon_Compton_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, l_integral)
                else
                    call Gauss_Legendre_quadrature_2D &
                    (MGXS_photon_photon_Compton_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, l_integral)
                end if
                integral(gpr,g,lp+1) = l_integral
            end do
        end do
    end do

    Sigmappc = 0
    do g = 1, Gp
        do gpr = 1, g
            do l = 0, NL
                do lp = 0, l
                    Sigmappc(gpr,g,l+1) = Sigmappc(gpr,g,l+1) + &
                    (e_radius**2)*(2*l+1)*e_mass_E*coeff(l,lp)*integral(gpr,g,lp+1)/(4*dEp(gpr))
                end do
            end do
        end do
    end do

    !outerloop: do gpr = 1, Gp
    !    do l = 0, NL
    !        portion(l,gpr) = sum(Sigmappc(gpr,:,l+1))
    !        if (portion(l,gpr)/portion(0,gpr) .lt. 0.01) then
    !            Sigmappc(gpr,:,l+1:NL+1) = 0.0
    !            cycle outerloop
    !        end if
    !    end do
    !end do outerloop

    do g = 1, Gp
        do gpr = g, Gp
            EDEP(g,:) = EDEP(g,:) - fourpi*mat_rho_e*Sigmappc(g,gpr,1)*0.5*(Ep(gpr)+Ep(gpr+1))
        end do
    end do
end subroutine MGXS_photon_photon_Compton

subroutine MGXS_photon_photon_fluorescence(Ep, dEp, SigmappF, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp

    real, dimension(:,:,:), allocatable, intent(inout) :: SigmappF
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: cascade, ele, g, gpr, i, j, mat, shell
    integer :: Nshells
    integer :: Ncascades
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: Wvals ! MeV
    real, dimension(:), allocatable :: trash_kinetic ! MeV
    integer, dimension(:), allocatable :: cscpershell
    integer, dimension(:), allocatable :: B_cscpershell
    integer, dimension(:), allocatable :: C_cscpershell
    real, dimension(:), allocatable :: Eij ! MeV
    real, dimension(:), allocatable :: B_Eij
    real, dimension(:), allocatable :: C_Eij
    real, dimension(:), allocatable :: eta
    real, dimension(:), allocatable :: B_eta
    real, dimension(:), allocatable :: C_eta
    real, dimension(:), allocatable :: peWvals ! MeV
    real, dimension(:), allocatable :: Espec
    real, dimension(:,:), allocatable :: pexs
    integer, dimension(:), allocatable :: matches
    integer, dimension(:,:), allocatable :: deltas
    integer :: globalcascade
    real :: l_alpha
    real :: l_beta
    real :: integral
    real, dimension(:), allocatable :: integrals

    allocate(SigmappF(Gp,Gp,Nmats))
    allocate(integrals(Gp))

    SigmappF = 0

    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            if (eles_list(mat,ele) .lt. 6) cycle

            ! READ PARAMETERS
            call read_atomic_parameters(eles_list(mat,ele), degen, Wvals, trash_kinetic)
            deallocate(trash_kinetic)

            ! READ AUGER DATA
            call read_relaxation_data("fluorescence", eles_list(mat,ele), size(Wvals), Eij, eta, cscpershell)

            ! READ PHOTOELECTRIC EFFECT DATA
            call read_photoelectric_effect(eles_list(mat,ele), Epmin, Epmax, peWvals, Espec, pexs)
            Nshells = size(peWvals)

            ! Determine which shells have photoelectric effect data, and which do not
            allocate(matches(Nshells))
            matches = 0
            do i = 1, Nshells
                do j = 1, size(Wvals)
                    if (abs(peWvals(i) - Wvals(j)) .le. 1.0E-6) then
                        matches(i) = j
                    end if
                end do
            end do

            if (any(matches .eq. 0)) then
                matches = 0
                do i = 1, Nshells
                    matches(i) = i
                end do
            end if

            do i = 1, Nshells
                peWvals(i) = Wvals(matches(i))
            end do
            call move_alloc(peWvals, Wvals)

            Ncascades = 0
            allocate(B_cscpershell(Nshells))
            do i = 1, Nshells
                Ncascades = Ncascades + cscpershell(matches(i))
                B_cscpershell(i) = cscpershell(matches(i))
            end do
            call move_alloc(B_cscpershell, cscpershell)

            allocate(B_Eij(0))
            allocate(B_eta(0))
            do i = 1, Nshells
                cascade = 1
                do j = 1, matches(i)-1
                    cascade = cascade + cscpershell(matches(j))
                end do
                allocate(C_Eij(size(B_Eij)+cscpershell(matches(i))))
                C_Eij(1:size(B_Eij)) = B_Eij
                C_Eij(size(B_Eij)+1:size(C_Eij)) = Eij(cascade:cascade+cscpershell(matches(i))-1)
                call move_alloc(C_Eij, B_Eij)
                allocate(C_eta(size(B_eta)+cscpershell(matches(i))))
                C_eta(1:size(B_eta)) = B_eta
                C_eta(size(B_eta)+1:size(C_eta)) = eta(cascade:cascade+cscpershell(matches(i))-1)
                call move_alloc(C_eta, B_eta)
            end do
            call move_alloc(B_Eij, Eij)
            call move_alloc(B_eta, eta)

            if (size(Eij) .ne. Ncascades .or. size(eta) .ne. Ncascades) then
                print *, "STOP!"
                print *, "MODULE MGXS, SUBROUTINE MGXS_photon_photon_fluorescence:"
                print *, "Number of cascades available for photoionization doesn't match"
                print *, "the number of energies or efficiencies after truncating from EADL data."
                print *, "PROGRAM ENDING."
                stop
            end if

            if (all(Eij .lt. Eemin)) go to 1

            allocate(deltas(Gp,Ncascades))
            do cascade = 1, Ncascades
                do g = 1, Gp
                    if (Ep(g+1) .le. Eij(cascade) .and. Eij(cascade) .le. Ep(g)) then
                        deltas(g,cascade) = 1
                    else
                        deltas(g,cascade) = 0
                    end if
                end do
            end do

            globalcascade = 0
            do shell = 1, Nshells
                ! INTEGRATE PEXS
                integrals = 0
                do gpr = 1, Gp
                    if (photon_exc_shells) then
                        if (Ep(gpr+1) + 1.0E-6 .lt. Wvals(shell)) cycle
                        l_alpha = Ep(gpr+1)/e_mass_E
                        l_beta = Ep(gpr)/e_mass_E
                    else
                        l_alpha = min(max(Ep(gpr+1),Wvals(shell)),Ep(gpr))/e_mass_E
                        l_beta = Ep(gpr)/e_mass_E
                    end if
                    call log_log_integration_1D &
                        (0, pexs(:,shell), Espec, l_alpha, l_beta, integrals(gpr))
                end do

                do cascade = 1, cscpershell(shell)
                    globalcascade = globalcascade + 1
                    do g = 1, Gp
                        do gpr = 1, Gp
                            SigmappF(gpr,g,mat) = SigmappF(gpr,g,mat) + &
                            mat_rho_a(mat,ele)*deltas(g,globalcascade)*&
                            eta(globalcascade)*e_mass_E*integrals(gpr)/(dEp(gpr)*fourpi)
                        end do
                    end do
                end do
            end do

            deallocate(deltas)
    1       deallocate(degen)
            deallocate(Wvals)
            deallocate(pexs)
            deallocate(Espec)
            deallocate(cscpershell)
            deallocate(Eij)
            deallocate(eta)
            deallocate(matches)
        end do
    end do

    do g = 1, Gp
        do gpr = 1, Gp
            EDEP(g,:) = EDEP(g,:) - fourpi*SigmappF(g,gpr,:)*0.5*(Ep(gpr)+Ep(gpr+1))
        end do
    end do
end subroutine MGXS_photon_photon_fluorescence

!--------------------- Photon-Electron ------------------------------------------------------------
subroutine MGXS_photon_electron(Ep, dEp, Ee, dEe, Eemid, adaptive, tol, Nint, Sigmape, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigmape
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: g, mat
    real, dimension(:,:,:), allocatable :: Sigmapec
    real, dimension(:,:,:,:), allocatable :: Sigmapepp
    real, dimension(:,:,:,:), allocatable :: Sigmapepe
    real, dimension(:,:,:), allocatable :: SigmapeA

    call MGXS_photon_electron_Compton(Ep, dEp, Ee, dEe, adaptive, tol, Nint, Sigmapec, EDEP)
    if (Epmax .lt. 1.022) then
        allocate(Sigmapepp(Gp,Ge,NL+1,Nmats))
        Sigmapepp = 0
    else
        call MGXS_photon_electron_pair_production(Ep, dEp, Ee, dEe, Eemid, Nint, Sigmapepp, EDEP)
    end if
    call MGXS_photon_electron_photoelectric_effect(Ep, dEp, Ee, dEe, Eemid, Sigmapepe, EDEP)
    call MGXS_photon_electron_Auger(Ep, dEp, Ee, dEe, SigmapeA, EDEP)

    allocate(Sigmape(Gp,Ge,NL+1,Nmats))
    Sigmape = 0
    Sigmape(1:Gp,1:Ge,1,1:Nmats) = SigmapeA

    do mat = 1, Nmats
        Sigmape(1:Gp,1:Ge,1:NL+1,mat) = Sigmape(1:Gp,1:Ge,1:NL+1,mat) + mat_rho_e(mat)*Sigmapec
    end do

    Sigmape = Sigmape + Sigmapepp + Sigmapepe
end subroutine MGXS_photon_electron

subroutine MGXS_photon_electron_Compton(Ep, dEp, Ee, dEe, adaptive, tol, Nint, Sigmapec, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:), allocatable, intent(inout) :: Sigmapec
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: g, gpr, l, lp
    real, dimension(:,:), allocatable :: coeff
    real, dimension(:,:,:), allocatable :: integral
    real :: l_integral
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:), allocatable :: portion

    allocate(Sigmapec(Gp,Ge,NL+1))
    allocate(coeff(0:NL,0:NL))
    allocate(integral(Gp,Ge,NL+1))
    allocate(portion(0:NL,Gp))

    coeff = Legendre_polynomial_coeffs(NL)

    integral = 0
    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)

        do gpr = 1, Gp
            XSa = Ep(gpr+1)/e_mass_E
            XSb = Ep(gpr)/e_mass_E
            do g = 1, Ge
                if (Ee(g+1) .ge. Emaxec(Ep(gpr))) cycle
                XSc = Ee(g+1)/e_mass_E
                XSd = Ee(g)/e_mass_E
                do lp = 0, NL
                    XSlp = lp
                    call Gauss_Legendre_adaptive_quadrature_2D &
                    (MGXS_photon_electron_Compton_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, l_integral)
                    integral(gpr,g,lp+1) = l_integral
                end do
            end do
        end do
    else
        call populate_Gauss_Legendre_quadrature(Nint, xq1, wq1)

        do gpr = 1, Gp
            XSa = Ep(gpr+1)/e_mass_E
            XSb = Ep(gpr)/e_mass_E
            do g = 1, Ge
                if (Ee(g+1) .ge. Emaxec(Ep(gpr))) cycle
                XSc = Ee(g+1)/e_mass_E
                XSd = Ee(g)/e_mass_E
                do lp = 0, NL
                    XSlp = lp
                    call Gauss_Legendre_quadrature_2D &
                    (MGXS_photon_electron_Compton_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, l_integral)
                    integral(gpr,g,lp+1) = l_integral
                end do
            end do
        end do
    end if
    Sigmapec = 0
    do g = 1, Ge
        do gpr = 1, Gp
            do l = 0, NL
                do lp = 0, l
                    Sigmapec(gpr,g,l+1) = Sigmapec(gpr,g,l+1) + &
                    (e_radius**2)*(2*l+1)*e_mass_E*coeff(l,lp)*integral(gpr,g,lp+1)/(4*dEp(gpr))
                end do
            end do
        end do
    end do

    !outerloop: do gpr = 1, Gp
    !    do l = 0, NL
    !        portion(l,gpr) = sum(Sigmapec(gpr,1:Ge,l+1))
    !        if (portion(l,gpr)/portion(0,gpr) .lt. 0.01) then
    !            Sigmapec(gpr,1:Ge,l+1:NL+1) = 0
    !            cycle outerloop
    !        end if
    !    end do
    !end do outerloop

    do g = 1, Gp
        do gpr = 1, Ge
            EDEP(g,:) = EDEP(g,:) - fourpi*mat_rho_e*Sigmapec(g,gpr,1)*0.5*(Ee(gpr)+Ee(gpr+1))
        end do
    end do
end subroutine MGXS_photon_electron_Compton

subroutine MGXS_photon_electron_pair_production(Ep, dEp, Ee, dEe, Eemid, Nint, Sigmapepp, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigmapepp
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: ele, mat, g, gpr, i, j, l
    integer :: NNIST
    integer :: Gpp
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:), allocatable :: ENist
    real, dimension(:), allocatable :: ppxs
    real, dimension(:,:), allocatable :: angular
    real, dimension(:,:,:), allocatable :: Fints
    real, dimension(:), allocatable :: norms
    real :: l_eps
    real :: l_delta
    real :: l_gamma
    real :: l_beta
    real :: l_E
    real, dimension(:,:,:), allocatable :: yvals
    real :: l_integral
    real, dimension(:,:,:), allocatable :: integrals

    Gpp = findloc(Ep, 1.022, dim = 1)
    if (Gpp .eq. 0) then
        print *, "MGXS pair production: Photon energy grouping does not contain 1.022 MeV (2*electron mass energy)."
        stop
    end if
    Gpp = Gpp - 1 ! Last group for which PP is valid

    allocate(Sigmapepp(Gp,Ge,NL+1,Nmats))
    allocate(angular(NL+1,Nint))
    allocate(integrals(Gpp,Ge,NL+1))

    Sigmapepp = 0

    call populate_Gauss_Legendre_quadrature(Nint, xq, wq)

    do mat = 1, Nmats
        integrals = 0
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            XSZ = eles_list(mat,ele)
            call read_pair_production(eles_list(mat,ele), Epmax, ENIST, ppxs)
            NNIST = size(ENIST)

            allocate(Fints(NNIST,NL+1,Ge))
            allocate(yvals(NNIST,NL+1,Ge))
            allocate(norms(NNIST))

            norms = 0
            do i = 1, NNIST
                do j = 1, Nint
                    l_eps = 0.5*(0.5 - 1/ENIST(i))*xq(j) + 0.5*(0.5 + 1/ENIST(i))
                    norms(i) = norms(i) + (0.5-1/ENIST(i))*wq(j)*Bethe_Heitler_DED(ENIST(i),l_eps)
                end do
            end do

            Fints = 0
            do g = 1, Ge
                l_gamma = Ee(g+1)/e_mass_E
                do i = 1, NNIST
                    l_delta = max(min(Ee(g)/e_mass_E,ENIST(i)-2),Ee(g+1)/e_mass_E)
                    do j = 1, Nint
                        l_eps = 0.5*((l_delta + 1)/ENIST(i) - (l_gamma + 1)/ENIST(i))*xq(j) &
                              + 0.5*((l_delta + 1)/ENIST(i) + (l_gamma + 1)/ENIST(i))
                        l_E = l_eps*ENIST(i) - 1
                        l_beta = unitless_e_beta(l_E)
                        angular(1,j) = 1
                        angular(2,j) = &
                            (2*l_beta + (1-l_beta**2)*log((1-l_beta)/(1+l_beta)))/(2*l_beta**2)
                        do l = 2, NL
                            angular(l+1,j) = &
                            ((2*l-1)*angular(l,j) - l*l_beta*angular(l-1,j))/((l-1)*l_beta)
                            if (l .ge. 6 .and. l_beta .lt. 0.5) angular(l+1,j) = 0
                        end do
                    end do

                    do l = 0, NL
                        do j = 1, Nint
                            l_eps = 0.5*((l_delta + 1)/ENIST(i) - (l_gamma + 1)/ENIST(i))*xq(j) &
                                  + 0.5*((l_delta + 1)/ENIST(i) + (l_gamma + 1)/ENIST(i))
                            Fints(i,l+1,g) = Fints(i,l+1,g) + &
                            0.5*((l_delta + 1)/ENIST(i) - (l_gamma + 1)/ENIST(i))*&
                            wq(j)*angular(l+1,j)*Bethe_Heitler_DED(ENIST(i),l_eps)/norms(i)
                        end do
                    end do
                end do
            end do

            do i = 1, NNIST
                yvals(i,1:NL+1,1:Ge) = Fints(i,1:NL+1,1:Ge)*ppxs(i)
            end do

            do g = 1, Ge
                do gpr = 1, Gpp
                    do l = 0, NL
                        !call trapezoidal_integration_1D &
                        !(0, yvals(:,l+1,g), ENIST, Ep(gpr+1)/e_mass_E, Ep(gpr)/e_mass_E, l_integral)
                        call log_log_integration_1D &
                        (0, yvals(:,l+1,g), ENIST, Ep(gpr+1)/e_mass_E, Ep(gpr)/e_mass_E, l_integral)
                        integrals(gpr,g,l+1) = integrals(gpr,g,l+1) + &
                            mat_rho(mat)*mass_frac(mat,ele)*e_mass_E*(2*l+1)*l_integral/(twopi*dEp(gpr))
                    end do
                end do
            end do
            deallocate(ENIST)
            deallocate(ppxs)
            deallocate(Fints)
            deallocate(yvals)
            deallocate(norms)
        end do

        Sigmapepp(1:Gpp,1:Ge,1:NL+1,mat) = integrals
    end do

    do g = 1, Gp
        do gpr = 1, Ge
            EDEP(g,:) = EDEP(g,:) - fourpi*Sigmapepp(g,gpr,1,:)*0.5*(Ee(gpr)+Ee(gpr+1))
        end do
    end do
end subroutine MGXS_photon_electron_pair_production

subroutine MGXS_photon_electron_photoelectric_effect(Ep, dEp, Ee, dEe, Eemid, Sigmapepe, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigmapepe
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: ele, g, gpr, i, ip, j, l, mat, shell
    integer :: Nint = 64
    integer :: NEspec
    integer :: Nshells
    real, dimension(:), allocatable :: E
    real, dimension(:), allocatable :: peWvals
    real, dimension(:,:), allocatable :: pexs
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:,:), allocatable :: Legendre
    real, dimension(:,:), allocatable :: angular
    real, dimension(:,:), allocatable :: apexs
    real :: l_alpha
    real :: l_beta
    real :: integral
    integer :: eof

    allocate(Sigmapepe(Gp,Ge,NL+1,Nmats))
    allocate(Legendre(Nint,NL+1))

    call populate_Gauss_Legendre_quadrature(Nint, xq, wq)

    Legendre(:,1:NL+1) = Legendre_polynomials(NL,xq)

    Sigmapepe = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            call read_photoelectric_effect(eles_list(mat,ele), Epmin, Epmax, peWvals, E, pexs)
            Nshells = size(peWvals)
            NEspec = size(E)

            ! COMPARE peWvals TO EADL Wvals

            allocate(angular(NEspec,NL+1))
            allocate(apexs(NEspec,NL+1))

            do shell = 1, Nshells
                angular = 0
                angular(:,1) = 1
                do l = 1, NL
                    do i = 1, NEspec
                        if (E(i) .le. peWvals(shell)/e_mass_E) cycle
                        do ip = 1, Nint
                            angular(i,l+1) = angular(i,l+1) + &
                                wq(ip)*Legendre(ip,l+1)*Sauter(xq(ip),E(i)-peWvals(shell)/e_mass_E)
                        end do
                    end do
                end do

                do l = 0, NL
                    apexs(:,l+1) = angular(:,l+1)*pexs(:,shell)
                end do

                ! Am I supposed to use total at all? Am I supposed to use probabilities of shell ionization?

                do gpr = 1, Gp
                    do g = 1, Ge
                        l_alpha = min(max(Ep(gpr+1),Ee(g+1)+peWvals(shell)),Ep(gpr))/e_mass_E
                        l_beta = max(min(Ep(gpr),Ee(g)+peWvals(shell)),Ep(gpr+1))/e_mass_E
                        if (l_alpha .eq. l_beta) cycle
                        do l = 0, NL
                            call trapezoidal_integration_1D &
                            (0, apexs(:,l+1), E, l_alpha, l_beta, integral)
                            Sigmapepe(gpr,g,l+1,mat) = Sigmapepe(gpr,g,l+1,mat) + &
                                (2*l+1)*mat_rho_a(mat,ele)*e_mass_E*integral/&
                                (fourpi*dEp(gpr))
                        end do
                    end do
                end do
            end do

            deallocate(peWvals)
            deallocate(E)
            deallocate(pexs)
            deallocate(angular)
            deallocate(apexs)
        end do
    end do

    do g = 1, Gp
        do gpr = 1, Ge
            EDEP(g,:) = EDEP(g,:) - fourpi*Sigmapepe(g,gpr,1,:)*0.5*(Ee(gpr)+Ee(gpr+1))
        end do
    end do
end subroutine MGXS_photon_electron_photoelectric_effect

subroutine MGXS_photon_electron_Auger(Ep, dEp, Ee, dEe, SigmapeA, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe

    real, dimension(:,:,:), allocatable, intent(inout) :: SigmapeA
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: cascade, ele, g, gpr, i, j, mat, shell
    integer :: Nshells
    integer :: Ncascades
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: Wvals ! MeV
    real, dimension(:), allocatable :: trash_kinetic ! MeV
    integer, dimension(:), allocatable :: cscpershell
    integer, dimension(:), allocatable :: B_cscpershell
    integer, dimension(:), allocatable :: C_cscpershell
    real, dimension(:), allocatable :: Eij ! MeV
    real, dimension(:), allocatable :: B_Eij
    real, dimension(:), allocatable :: C_Eij
    real, dimension(:), allocatable :: eta
    real, dimension(:), allocatable :: B_eta
    real, dimension(:), allocatable :: C_eta
    real, dimension(:), allocatable :: peWvals ! MeV
    real, dimension(:), allocatable :: Espec
    real, dimension(:,:), allocatable :: pexs
    integer, dimension(:), allocatable :: matches
    integer, dimension(:,:), allocatable :: deltas
    integer :: globalcascade
    real :: l_alpha
    real :: l_beta
    real :: integral
    real, dimension(:), allocatable :: integrals

    allocate(SigmapeA(Gp,Ge,Nmats))
    allocate(integrals(Gp))

    SigmapeA = 0

    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            if (eles_list(mat,ele) .lt. 6) cycle

            ! READ PARAMETERS
            call read_atomic_parameters(eles_list(mat,ele), degen, Wvals, trash_kinetic)
            deallocate(trash_kinetic)

            ! READ AUGER DATA
            call read_relaxation_data("auger", eles_list(mat,ele), size(Wvals), Eij, eta, cscpershell)

            ! READ PHOTOELECTRIC EFFECT DATA
            call read_photoelectric_effect(eles_list(mat,ele), Epmin, Epmax, peWvals, Espec, pexs)
            Nshells = size(peWvals)

            ! Determine which shells have photoelectric effect data, and which do not
            allocate(matches(Nshells))
            matches = 0
            do i = 1, Nshells
                do j = 1, size(Wvals)
                    if (abs(peWvals(i) - Wvals(j)) .le. 1.0E-6) then
                        matches(i) = j
                    end if
                end do
            end do

            if (any(matches .eq. 0)) then
                matches = 0
                do i = 1, Nshells
                    matches(i) = i
                end do
            end if

            do i = 1, Nshells
                peWvals(i) = Wvals(matches(i))
            end do
            call move_alloc(peWvals, Wvals)

            Ncascades = 0
            allocate(B_cscpershell(Nshells))
            do i = 1, Nshells
                Ncascades = Ncascades + cscpershell(matches(i))
                B_cscpershell(i) = cscpershell(matches(i))
            end do
            call move_alloc(B_cscpershell, cscpershell)

            allocate(B_Eij(0))
            allocate(B_eta(0))
            do i = 1, Nshells
                cascade = 1
                do j = 1, matches(i)-1
                    cascade = cascade + cscpershell(matches(j))
                end do
                allocate(C_Eij(size(B_Eij)+cscpershell(matches(i))))
                C_Eij(1:size(B_Eij)) = B_Eij
                C_Eij(size(B_Eij)+1:size(C_Eij)) = Eij(cascade:cascade+cscpershell(matches(i))-1)
                call move_alloc(C_Eij, B_Eij)
                allocate(C_eta(size(B_eta)+cscpershell(matches(i))))
                C_eta(1:size(B_eta)) = B_eta
                C_eta(size(B_eta)+1:size(C_eta)) = eta(cascade:cascade+cscpershell(matches(i))-1)
                call move_alloc(C_eta, B_eta)
            end do
            call move_alloc(B_Eij, Eij)
            call move_alloc(B_eta, eta)

            if (size(Eij) .ne. Ncascades .or. size(eta) .ne. Ncascades) then
                print *, "STOP!"
                print *, "MODULE MGXS, SUBROUTINE MGXS_photon_electron_Auger:"
                print *, "Number of cascades available for photoionization doesn't match"
                print *, "the number of energies or efficiencies after truncating from EADL data."
                print *, "PROGRAM ENDING."
                stop
            end if

            if (all(Eij .lt. Eemin)) go to 1

            allocate(deltas(Ge,Ncascades))
            do cascade = 1, Ncascades
                do g = 1, Ge
                    if (Ee(g+1) .le. Eij(cascade) .and. Eij(cascade) .le. Ee(g)) then
                        deltas(g,cascade) = 1
                    else
                        deltas(g,cascade) = 0
                    end if
                end do
            end do

            globalcascade = 0
            do shell = 1, Nshells
                ! INTEGRATE PEXS
                integrals = 0
                do gpr = 1, Gp
                    if (photon_exc_shells) then
                        if (Ep(gpr+1) + 1.0E-6 .lt. Wvals(shell)) cycle
                        l_alpha = Ep(gpr+1)/e_mass_E
                        l_beta = Ep(gpr)/e_mass_E
                    else
                        l_alpha = min(max(Ep(gpr+1),Wvals(shell)),Ep(gpr))/e_mass_E
                        l_beta = Ep(gpr)/e_mass_E
                    end if
                    call log_log_integration_1D &
                        (0, pexs(:,shell), Espec, l_alpha, l_beta, integrals(gpr))
                end do

                do cascade = 1, cscpershell(shell) ! Could double check this does what I want it to, ~4/22/23
                    globalcascade = globalcascade + 1
                    do g = 1, Ge
                        do gpr = 1, Gp
                            SigmapeA(gpr,g,mat) = SigmapeA(gpr,g,mat) + &
                            mat_rho_a(mat,ele)*deltas(g,globalcascade)*&
                            eta(globalcascade)*e_mass_E*integrals(gpr)/(fourpi*dEp(gpr))
                        end do
                    end do
                end do
            end do

            deallocate(deltas)
    1       deallocate(degen)
            deallocate(Wvals)
            deallocate(cscpershell)
            deallocate(Eij)
            deallocate(eta)
            deallocate(matches)
        end do
    end do

    do g = 1, Gp
        do gpr = 1, Ge
            EDEP(g,:) = EDEP(g,:) - fourpi*SigmapeA(g,gpr,:)*0.5*(Ee(gpr)+Ee(gpr+1))
        end do
    end do
end subroutine MGXS_photon_electron_Auger

!--------------------- Electron-Electron ----------------------------------------------------------
subroutine MGXS_electron_electron &
    (Ee, dEe, Eemid, adaptive, tol, N, S, St, Sa, S2, EDEP, CDEP)
    ! Change name to MGXS_electron? Gotta figure out what to specialize
    ! for case with no photons and with photons. Like fluorescence.
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: N

    real, dimension(:,:,:,:), allocatable, intent(inout) :: S
    real, dimension(:,:), allocatable, intent(inout) :: St
    real, dimension(:,:), allocatable, intent(inout) :: Sa
    real, dimension(:,:), allocatable, intent(inout) :: S2
    real, dimension(:,:), allocatable, intent(inout) :: EDEP
    real, dimension(:,:), allocatable, intent(inout) :: CDEP

    integer :: g, gpr, l, mat
    integer :: l_NL
    ! Inelastic
    real, dimension(:,:,:,:), allocatable :: S_inel
    real, dimension(:,:), allocatable :: St_inel
    real, dimension(:,:), allocatable :: S2_inel
    ! Bremsstrahlung
    real, dimension(:,:,:), allocatable :: S_B
    real, dimension(:,:), allocatable :: St_B
    ! RCSDA
    real, dimension(:,:,:), allocatable :: R
    real, dimension(:,:), allocatable :: Rt
    real, dimension(:,:), allocatable :: Ra
    real, dimension(:,:,:,:), allocatable :: RTEST
    ! Elastic
    real, dimension(:,:,:), allocatable :: S_el
    real, dimension(:,:), allocatable :: St_el
    ! Relaxation
    real, dimension(:,:,:), allocatable :: S_A
    real, dimension(:,:,:), allocatable :: S_F

    XSlogic = .false.

    l_NL = NL
    if (exact_RCSDA_angular .or. inelastic_etc) l_NL = NL + 1

    allocate(S(Ge,Ge,l_NL+1,Nmats))
    allocate(Sa(Ge,Nmats))
    allocate(St(Ge,Nmats))
    allocate(S2(Ge,Nmats))
    allocate(EDEP(Ge,Nmats))
    allocate(CDEP(Ge,Nmats))
    S = 0
    Sa = 0
    St = 0
    S2 = 0
    EDEP = 0
    CDEP = 0

    ! --------------
    ! ### SCATTERING
    ! --------------

    ! Inelastic
    call MGXS_electron_electron_inelastic &
        (Ee, dEe, adaptive, tol, N, S_inel, St_inel, S2_inel)

    !! Bremsstrahlung
    !call MGXS_electron_electron_Bremsstrahlung(Ee, dEe, adaptive, tol, N, S_B, St_B)

    ! Elastic
    if (electron_elastic .eq. "fitted") then
        call MGXS_electron_electron_elastic_fitted(Ee, dEe, adaptive, tol, N, S_el, St_el)
    else if (electron_elastic .eq. "Wentzel-Moliere") then
        call MGXS_electron_electron_elastic_Wentzel_Moliere &
            (Ee, dEe, adaptive, tol, N, S_el, St_el)
    else if (electron_elastic .eq. "ELSEPA") then
        call MGXS_electron_electron_elastic_ELSEPA(Ee, dEe, S_el, St_el)
    end if

    ! RCSDA
    if (electron_RCSDA .eq. "first order") then
        call MGXS_electron_electron_FO_RCSDA(Ee, dEe, Eemid, R, Rt)
    else if (electron_RCSDA .eq. "second order") then
        call MGXS_electron_electron_SO_RCSDA(Ee, dEe, Eemid, R, Rt, Ra)
    else
        call MGXS_electron_electron_RCSDA(Ee, dEe, Eemid, adaptive, tol, N, R, Rt)
    end if

    !! Relaxation
    !call MGXS_electron_electron_Auger(Ee, dEe, adaptive, tol, N, S_A)
    !!call MGXS_electron_photon_fluorescence ! NEED TO DO THIS. IT CONTRIBUTES TO EDEP.

    ! ----------------------------------------
    ! ### ABSORPTION, SECONDARY, DEPOSITION XS
    ! ----------------------------------------

    ! Inelastic (T, A, S, E, C)
    St = St + St_inel
    S2 = S2 + S2_inel
    do g = 1, Ge
        Sa(g,:) = Sa(g,:) + St_inel(g,:) - S2_inel(g,:)
        EDEP(g,:) = EDEP(g,:) + St_inel(g,:)*Eemid(g)
        CDEP(g,:) = CDEP(g,:) + St_inel(g,:)
        do gpr = g+2, Ge
            Sa(g,:) = Sa(g,:) - fourpi*S_inel(g,gpr,1,:)
            EDEP(g,:) = EDEP(g,:) - fourpi*S_inel(g,gpr,1,:)*Eemid(gpr)
            CDEP(g,:) = CDEP(g,:) - fourpi*S_inel(g,gpr,1,:)
        end do
    end do

    !! Bremsstrahlung (T, A, S, E, C) ! Double check absorption, secondary, etc. No 4pi because assumed no angular deflection.
    !St = St + St_B
    !do g = 1, Ge
    !    Sa(g,:) = Sa(g,:) + St_B(g,:)
    !    EDEP(g,:) = EDEP(g,:) + St_B(g,:)*Eemid(g)
    !    CDEP(g,:) = CDEP(g,:) + St_B(g,:)
    !    do gpr = g+2, Ge
    !        Sa(g,:) = Sa(g,:) - S_B(g,gpr,:)
    !        EDEP(g,:) = EDEP(g,:) - S_B(g,gpr,:)*Eemid(gpr)
    !        CDEP(g,:) = CDEP(g,:) - S_B(g,gpr,:)
    !    end do
    !end do

    ! RCSDA (T, A, E, C)
    St = St + Rt
    if (electron_RCSDA .eq. "first order") then
        Sa(Ge,:) = Sa(Ge,:) + Rt(Ge,:)
        do g = 1, Ge-1
            EDEP(g,:) = EDEP(g,:) + R(g,g+1,:)*(Eemid(g) - Eemid(g+1))
        end do
        EDEP(Ge,:) = EDEP(Ge,:) + Rt(Ge,:)*Eemid(Ge)
        CDEP(Ge,:) = CDEP(Ge,:) + Rt(Ge,:)
    else if (electron_RCSDA .eq. "second order") then
        Sa = Sa + Ra
        do g = 1, Ge
            EDEP(g,:) = EDEP(g,:) + Rt(g,:)*Eemid(g)
            do gpr = g+1, Ge
                EDEP(g,:) = EDEP(g,:) - R(g,gpr,:)*Eemid(gpr)
            end do
        end do
        CDEP = CDEP + Ra
    else
        do g = 1, Ge
            Sa(g,:) = Sa(g,:) + Rt(g,:)
            EDEP(g,:) = EDEP(g,:) + Rt(g,:)*Eemid(g)
            CDEP(g,:) = CDEP(g,:) + Rt(g,:)
            if (g+1 .gt. Ge) cycle
            do gpr = g+1, merge(Ge, g+2, g+2 .gt. Ge)
                Sa(g,:) = Sa(g,:) - R(g,gpr,:)
                EDEP(g,:) = EDEP(g,:) - R(g,gpr,:)*Eemid(gpr)
                CDEP(g,:) = CDEP(g,:) - R(g,gpr,:)
            end do
        end do
    end if

    !! Relaxation (S, E, C)
    !do g = 1, Ge
    !    do gpr = g, Ge
    !        S2(g,:) = S2(g,:) + fourpi*S_A(g,gpr,:) ! Should use S_F
    !        EDEP(g,:) = EDEP(g,:) - fourpi*S_A(g,gpr,:)*Eemid(gpr) ! Should use S_F
    !        CDEP(g,:) = CDEP(g,:) - fourpi*S_A(g,gpr,:) ! Doesn't use S_F (it's the only one).
    !    end do
    !end do

    ! Elastic (T)
    St = St + St_el

    !S(1:Ge,1:Ge,1,1:Nmats) = S_A
    if (inelastic_etc) then
        do l = 0, NL
            S(:,:,l+1,:) = S_inel(:,:,l+1,:) - (2*l+1)*S_inel(:,:,NL+2,:)/(2*NL+3)
        end do

        S(:,:,NL+2,:) = fourpi*S_inel(:,:,NL+2,:)/(2*NL+3)

        if (exact_RCSDA_angular) then
            S(:,:,NL+2,:) = S(:,:,NL+2,:) + R
        else
            do l = 0, NL
                S(:,:,l+1,:) = S(:,:,l+1,:) + (2*l+1)*R/fourpi
            end do
        end if

        do g = 1, Ge
            S(g,g,1:NL+1,:) = S(g,g,1:NL+1,:) + S_el(g,:,:)
        end do
    else
        if (exact_RCSDA_angular) then
            S(:,:,NL+2,:) = S(:,:,NL+2,:) + R
        else
            do l = 0, NL
                S(:,:,l+1,:) = S(:,:,l+1,:) + (2*l+1)*R/fourpi
            end do
        end if

        S(:,:,1:NL+1,:) = S(:,:,1:NL+1,:) + S_inel

        do g = 1, Ge
            S(g,g,1:NL+1,:) = S(g,g,1:NL+1,:) + S_el(g,1:NL+1,:)
        end do
    end if
end subroutine MGXS_electron_electron

subroutine MGXS_electron_electron_inelastic(Ee, dEe, adaptive, tol, N, S_inel, St_inel, S2_inel)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: N

    real, dimension(:,:,:,:), allocatable, intent(inout) :: S_inel
    real, dimension(:,:), allocatable, intent(inout) :: St_inel
    real, dimension(:,:), allocatable, intent(inout) :: S2_inel

    integer :: ele, g, gpr, l, lp, mat, s
    integer :: l_NL
    integer :: Nshells
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:), allocatable :: integrals
    real :: integral
    real :: integral1
    real :: integral2
    real :: secondary
    real, dimension(:,:), allocatable :: coeff

    l_NL = NL
    if (inelastic_etc) l_NL = NL + 1

    allocate(S_inel(Ge,Ge,l_NL+1,Nmats))
    allocate(St_inel(Ge,Nmats))
    allocate(S2_inel(Ge,Nmats))
    allocate(integrals(l_NL+1))
    allocate(coeff(0:l_NL,0:l_NL))

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(N, xq1, wq1)
    end if

    coeff = Legendre_polynomial_coeffs(l_NL)

    S_inel = 0
    St_inel = 0
    S2_inel = 0
    if (electron_inelastic .eq. "Moller") then
        ! ccut can not be zero.
        XSw = 0
        do gpr = 1, Ge-1
            XSa = Ee(gpr+1)/e_mass_E
            XSb = Ee(gpr)/e_mass_E
            XScut = Ee(gpr+2)/e_mass_E
            do g = gpr+2, Ge
                XSc = Ee(g+1)/e_mass_E
                XSd = Ee(g)/e_mass_E
                do lp = 0, l_NL
                    XSlp = lp
                    if (adaptive) then
                        call Gauss_Legendre_adaptive_quadrature_2D &
                        (MGXS_electron_electron_Moller_primary_integrand, &
                        -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral1)
                        !call Gauss_Legendre_adaptive_quadrature_2D &
                        !(MGXS_electron_electron_Moller_secondary_integrand, &
                        !-1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral2)
                    else
                        call Gauss_Legendre_quadrature_2D &
                        (MGXS_electron_electron_Moller_primary_integrand, &
                        -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral1)
                        !call Gauss_Legendre_quadrature_2D &
                        !(MGXS_electron_electron_Moller_secondary_integrand, &
                        !-1.0, 1.0, -1.0, 1.0, xq1, wq1, integral2)
                    end if
                    if (adaptive) then
                        !call Gauss_Legendre_adaptive_quadrature_2D &
                        !(MOLLER_TEST_PRIMARY, &
                        !-1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral1)
                        call Gauss_Legendre_adaptive_quadrature_2D &
                        (MOLLER_TEST_SECONDARY, &
                        -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral2)
                    else
                        !call Gauss_Legendre_quadrature_2D &
                        !(MOLLER_TEST_PRIMARY, &
                        !-1.0, 1.0, -1.0, 1.0, xq1, wq1, integral1)
                        call Gauss_Legendre_quadrature_2D &
                        (MOLLER_TEST_SECONDARY, &
                        -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral2)
                    end if
                    !if (abs(integral1-integral2) .gt. 1.0E-10) then
                    !    print *, gpr, g, integral1-integral2
                    !end if

                    !if (gpr .eq. 31 .and. g .eq. 41) then
                    !    print *, "---------------"
                    !    print *, gpr, g
                    !    print *, XSlp
                    !    print *, XSa, XSb
                    !    print *, XSc, XSd
                    !    print *, XScut
                    !    print *, integral1
                    !    print *, "---------------"
                    !end if

                    integrals(lp+1) = integral1 + integral2
                    if (lp .eq. 0) secondary = integral2
                end do

                do l = 0, l_NL
                    do lp = 0, l
                        do mat = 1, Nmats
                            S_inel(gpr,g,l+1,mat) = S_inel(gpr,g,l+1,mat) + &
                                mat_rho_e(mat)*(e_radius**2)*coeff(l,lp)*&
                                (2*l+1)*e_mass_E*integrals(lp+1)/(2*dEe(gpr))
                        end do
                    end do
                end do

                S2_inel(gpr,:) = S2_inel(gpr,:) + &
                    fourpi*mat_rho_e*(e_radius**2)*e_mass_E*integral2/(2*dEe(gpr))
            end do

            !if (adaptive) then
            !    call Gauss_Legendre_adaptive_quadrature_1D &
            !    (MGXS_electron_Moller_R_total_integrand, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
            !else
            !    call Gauss_Legendre_quadrature_1D &
            !    (MGXS_electron_Moller_R_total_integrand, -1.0, 1.0, xq1, wq1, integral)
            !end if

            call Gauss_Legendre_quadrature_2D &
                (MOLLER_TEST_TOTAL, &
                -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral)

            !if (abs(integral - integral2) .gt. 1.0E-10) then
            !    print *, gpr, integral, integral2
            !end if

            ! The last group has no catastrophic attn., so this is fine to be in this loop
            St_inel(gpr,:) = mat_rho_e*twopi*(e_radius**2)*e_mass_E*integral/dEe(gpr)
        end do
    else if (electron_inelastic .eq. "RBED") then
        do mat = 1, Nmats
            do ele = 1, Nmaxeles
                if (eles_list(mat,ele) .eq. 0) cycle
                call read_atomic_parameters(eles_list(mat,ele), degen, W, U)
                Nshells = size(degen)

                do s = 1, Nshells
                    XSw = W(s)/e_mass_E
                    XSbeta_U2 = beta(U(s))**2
                    XSbeta_W2 = beta(W(s))**2
                    do gpr = 1, Ge-1
                        XSa = Ee(gpr+1)/e_mass_E
                        XSb = Ee(gpr)/e_mass_E
                        XScut = Ee(gpr+2)/e_mass_E
                        do g = gpr+2, Ge
                            XSc = Ee(g+1)/e_mass_E
                            XSd = Ee(g)/e_mass_E
                            do lp = 0, l_NL
                                XSlp = lp
                                if (adaptive) then
                                    call Gauss_Legendre_adaptive_quadrature_2D &
                                    (MGXS_electron_electron_RBED_primary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral1)
                                    call Gauss_Legendre_adaptive_quadrature_2D &
                                    (MGXS_electron_electron_RBED_secondary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral2)
                                else
                                    call Gauss_Legendre_quadrature_2D &
                                    (MGXS_electron_electron_RBED_primary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral1)
                                    call Gauss_Legendre_quadrature_2D &
                                    (MGXS_electron_electron_RBED_secondary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral2)
                                end if
                                integrals(lp+1) = integral1 + integral2
                                if (lp .eq. 0) secondary = integral2
                            end do

                            do l = 0, l_NL
                                do lp = 0, l
                                    S_inel(gpr,g,l+1,mat) = S_inel(gpr,g,l+1,mat) + &
                                        degen(s)*mat_rho_a(mat,ele)*(e_radius**2)*coeff(l,lp)*&
                                        (2*l+1)*e_mass_E*integrals(lp+1)/(2*dEe(gpr))
                                end do
                            end do

                            S2_inel(gpr,mat) = S2_inel(gpr,mat) + &
                                degen(s)*fourpi*mat_rho_a(mat,ele)*&
                                (e_radius**2)*e_mass_E*secondary/(2*dEe(g))
                        end do

                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_1D &
                            (MGXS_electron_RBED_total_integrand, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                        else
                            call Gauss_Legendre_quadrature_1D &
                            (MGXS_electron_RBED_total_integrand, -1.0, 1.0, xq1, wq1, integral)
                        end if

                        ! The last group has no catastrophic attn., so this is fine to be in this loop
                        St_inel(gpr,mat) = St_inel(gpr,mat) + &
                            degen(s)*mat_rho_a(mat,ele)*twopi*(e_radius**2)*e_mass_E*integral/dEe(gpr)
                    end do
                end do

                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do
        end do
    end if
end subroutine MGXS_electron_electron_inelastic

subroutine MGXS_electron_electron_FO_RCSDA(Ee, dEe, Eemid, R, Rt)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid

    real, dimension(:,:,:), allocatable, intent(inout) :: R
    real, dimension(:,:), allocatable, intent(inout) :: Rt

    integer :: ele, g, gpr, mat, s
    integer :: Nshells
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    real :: Ibar
    real :: x
    real, dimension(:), allocatable :: l_LET
    real :: E_u
    real :: temp

    E_u = Ee(Ge+1) - dEe(Ge) + 0.5*dEe(Ge)

    allocate(R(Ge,Ge,Nmats))
    allocate(Rt(Ge,Nmats))
    allocate(l_LET(Ge))

    if (electron_inelastic .eq. "MMS") then
        deallocate(R, Rt, l_LET)
        allocate(R(size(Ee)-1,size(Ee)-1,1))
        allocate(l_LET(size(Ee)-1))
        allocate(Rt(size(Ee)-1,Nmats))

        R = 0

        do g = 1, size(Ee)-1
            l_LET(g) = RCSDA_MMS_LET(Eemid(g))
        end do

        do g = 2, size(Ee)-1
            R(g-1,g,1) = l_LET(g-1)/(Eemid(g-1)-Eemid(g))
        end do

        !Rt(1:size(Ee)-1-1,1) = Rt(1:size(Ee)-1-1,1) + &
        !    sum(R(1:size(Ee)-1-1,1:size(Ee)-1,1),2)
        do g = 1, size(Ee)-2
            Rt(g,1) = Rt(g,1) + R(g,g+1,1)
        end do
        !Rt(size(Ee)-1,1) = Rt(size(Ee)-1,1) + &
        !    l_LET(size(Ee)-1)/(Eemid(Rt(Ee)-1)-E_u)
        Rt(size(Ee)-1,1) = Rt(size(Ee)-1,1) + &
            interp1d(Eemid(size(Ee)-1-1:1:-1), Rt(size(Ee)-1-1:1:-1,1), Eemid(size(Ee)-1))

        return
    end if

    R = 0
    Rt = 0
    do mat = 1, Nmats
        call mean_excitation_energy(eles_list(mat,:), Natoms(mat,:), Ibar)
        call set_density_effect_parameters &
            (mat_rho(mat), total_GAM(mat), eles_list(mat,:), Natoms(mat,:), Ibar, .false.)

        l_LET = 0
        if (electron_inelastic .eq. "Moller") then
            XSw = 0
            XSI = Ibar*1.0E-6/e_mass_E

            do g = 1, Ge
                XScut = merge(0.0, Ee(g+2)/e_mass_E, g .eq. Ge)
                x = Eemid(g)/e_mass_E

                if (g .eq. Ge) then
                    l_LET(g) = mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*&
                        unitless_Bethe(x)
                    exit
                end if
                l_LET(g) = mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*&
                    unitless_Moller_RLET(x)
            end do
        else if (electron_inelastic .eq. "RBED") then
            do ele = 1, Nmaxeles
                if (eles_list(mat,ele) .eq. 0) cycle
                call read_atomic_parameters(eles_list(mat,ele), degen, W, U)
                Nshells = size(degen)

                do g = 1, Ge-1
                    XScut = Ee(g+2)/e_mass_E
                    x = Eemid(g)/e_mass_E

                    do s = 1, Nshells
                        XSw = W(s)/e_mass_E
                        XSbeta_W2 = beta(W(s))**2
                        XSbeta_U2 = beta(U(s))**2

                        l_LET(g) = l_LET(g) + &
                            mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius)**2*e_mass_E*&
                            unitless_RBED_RLET(x)
                    end do
                end do

                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do
        end if

        do g = 1, Ge-1
            R(g,g+1,mat) = l_LET(g)/(Eemid(g)-Eemid(g+1))
            Rt(g,mat) = R(g,g+1,mat)
        end do

        Rt(Ge,mat) = l_LET(Ge)/(Eemid(Ge)-E_u)
    end do
end subroutine MGXS_electron_electron_FO_RCSDA

subroutine MGXS_electron_electron_SO_RCSDA(Ee, dEe, Eemid, R, Rt, Ra)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid

    real, dimension(:,:,:), allocatable, intent(inout) :: R
    real, dimension(:,:), allocatable, intent(inout) :: Rt
    real, dimension(:,:), allocatable, intent(inout) :: Ra

    integer :: g, gpr, mat
    real :: Ibar
    real :: x
    !real, dimension(:), allocatable :: xq1
    !real, dimension(:), allocatable :: wq1
    !real, dimension(:), allocatable :: xq2
    !real, dimension(:), allocatable :: wq2
    real, dimension(:), allocatable :: Lupper
    real, dimension(:), allocatable :: Llower
    real, dimension(:), allocatable :: Llast
    real, dimension(:,:), allocatable :: LDelta

    allocate(R(Ge,Ge,Nmats))
    allocate(Rt(Ge,Nmats))
    allocate(Ra(Ge,Nmats))
    allocate(Lupper(Ge))
    allocate(Llower(Ge))
    allocate(Llast(Nmats))
    allocate(LDelta(Ge,Nmats))

    R = 0
    Rt = 0
    Ra = 0

    if (electron_inelastic .eq. "MMS") then
        deallocate(R, Rt, Lupper, Llower, Llast, LDelta)
        allocate(R(size(Ee)-1,size(Ee)-1,1))
        allocate(Rt(size(Ee)-1,1))
        allocate(Lupper(size(Ee)-1))
        allocate(Llower(size(Ee)-1))
        allocate(Llast(1))
        allocate(LDelta(size(Ee)-1,1))

        do g = 2, size(Ee)-1-1
            Lupper(g) = RCSDA_MMS_LET(Ee(g))
        end do
        do g = 2, size(Ee)-1
            Llower(g) = RCSDA_MMS_LET(Ee(g))
        end do

        Llast(1) = RCSDA_MMS_LET(Ee(size(Ee)))

        do g = 1, size(Ee)-1
            LDelta(g,1) = 0.5*(Lupper(g)+Llower(g))
        end do

        do gpr = 1, size(Ee)-1-1
            do g = gpr+1, size(Ee)-1-1
                R(gpr,g,:) = (-1)**(g-gpr+1)*2*(LDelta(g,:) + LDelta(g+1,:))/dEe(gpr)
            end do
            R(gpr,size(Ee)-1,:) = (-1)**(size(Ee)-1-gpr+1)*2*(LDelta(size(Ee)-1,:) + Llast)/dEe(gpr)
        end do

        do g = 1, size(Ee)-1
            Rt(g,:) = Rt(g,:) + (-1)**(size(Ee)-1-g)*2*Llast/dEe(g)

            do gpr = g+1, size(Ee)-1
                Rt(g,:) = Rt(g,:) + R(g,gpr,:)
            end do
        end do

        return
    end if

    if (electron_inelastic .eq. "Moller") then
        XSw = 0

        do mat = 1, Nmats
            call mean_excitation_energy(eles_list(mat,:), Natoms(mat,:), Ibar)
            call set_density_effect_parameters &
                (mat_rho(mat), total_GAM(mat), eles_list(mat,:), Natoms(mat,:), Ibar, .false.)

            Lupper = 0
            Llower = 0
            XSI = Ibar*1.0E-6

            do g = 2, Ge-1
                XScut = Ee(g+2)/e_mass_E
                x = Ee(g)/e_mass_E

                Lupper(g) = mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*&
                    unitless_Moller_RLET(x)
            end do
            do g = 2, Ge
                XScut = Ee(g+1)/e_mass_E
                x = Ee(g)/e_mass_E

                Llower(g) = mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*&
                    unitless_Moller_RLET(x)
            end do

            x = Ee(Ge+1)/e_mass_E
            Llast(mat) = mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*unitless_Bethe(x)

            do g = 2, Ge
                LDelta(g,mat) = 0.5*(Lupper(g)+Llower(g))
            end do
        end do
    else if (electron_inelastic .eq. "RBED") then
        print *, "RBED SO RCSDA WIP"
        stop
    end if

    !!! NOW DO BREMSSTRAHLUNG. JUST ADDS TO THE L'S.

    do gpr = 1, Ge-1
        do g = gpr+1, Ge-1
            R(gpr,g,:) = (-1)**(g-gpr+1)*2*(LDelta(g,:) + LDelta(g+1,:))/dEe(gpr)
        end do
        R(gpr,Ge,:) = (-1)**(Ge-gpr+1)*2*(LDelta(Ge,:) + Llast)/dEe(gpr)
    end do

    do g = 1, Ge
        Ra(g,:) = (-1)**(Ge-g)*2*Llast/dEe(g)

        Rt(g,:) = Ra(g,:) + sum(R(g,g:Ge,:))
    end do
end subroutine MGXS_electron_electron_SO_RCSDA

subroutine MGXS_electron_electron_RCSDA(Ee, dEe, Eemid, adaptive, tol, N, R, Rt)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: N

    real, dimension(:,:,:), allocatable, intent(inout) :: R
    real, dimension(:,:), allocatable, intent(inout) :: Rt

    integer :: ele, g, gpr, mat, s
    integer :: Nshells
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    real :: Ibar
    real :: plasma_freq
    real :: Cparam
    real :: Bparam
    real :: x1param
    real :: x0param
    real :: xparam
    real :: Dcorr
    real :: x
    real, dimension(:,:), allocatable :: l_LET
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real :: l_alpha
    real :: l_beta
    real :: integral

    allocate(R(Ge,Ge,Nmats))
    allocate(Rt(Ge,Nmats))
    allocate(l_LET(Ge,Nmats))

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(N, xq1, wq1)
    end if

    if (electron_inelastic .eq. "MMS") then
        deallocate(R, Rt, l_LET)
        allocate(R(size(Ee)-1,size(Ee)-1,1))
        allocate(Rt(size(Ee)-1,1))
        allocate(l_LET(size(Ee)-1,1))

        R = 0
        l_LET = 0

        do g = 1, size(Ee)-1
            if (adaptive) then
                call Gauss_Legendre_adaptive_quadrature_1D &
                (RCSDA_MMS_LET, Ee(g+1), Ee(g), xq1, wq1, xq2, wq2, tol, integral)
            else
                call Gauss_Legendre_quadrature_1D &
                (RCSDA_MMS_LET, Ee(g+1), Ee(g), xq1, wq1, integral)
            end if

            l_LET(g,1) = integral/dEe(g)
        end do

        do g = 1, size(Ee)-1
            Rt(g,:) = Rt(g,:) + 3*l_LET(g,:)/(2*dEe(g))
        end do

        R(1,2,:) = 2*l_LET(1,:)/dEe(1)
        do g = 3, size(Ee)-1
            R(g-1,g,:) = 2*l_LET(g-1,:)/dEe(g-1)
            R(g-2,g,:) = -0.5*l_LET(g-2,:)/dEe(g-2)
        end do

        return
    end if

    R = 0
    Rt = 0
    l_LET = 0
    if (electron_inelastic .eq. "Moller") then
        do mat = 1, Nmats
            call mean_excitation_energy(eles_list(mat,:), Natoms(mat,:), Ibar)
            call set_density_effect_parameters &
                (mat_rho(mat), total_GAM(mat), eles_list(mat,:), Natoms(mat,:), Ibar, .false.)

            XSw = 0
            XSI = Ibar*1.0E-6/e_mass_E
            do g = 1, Ge
                XScut = merge(0.0, Ee(g+2)/e_mass_E, g .eq. Ge)
                l_alpha = Ee(g+1)/e_mass_E
                l_beta = Ee(g)/e_mass_E

                !if (adaptive) then
                !    call Gauss_Legendre_adaptive_quadrature_1D &
                !    (unitless_Bethe, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                !else
                !    call Gauss_Legendre_quadrature_1D &
                !    (unitless_Bethe, l_alpha, l_beta, xq1, wq1, integral)
                !end if

                if (adaptive) then
                    call Gauss_Legendre_adaptive_quadrature_1D &
                    (unitless_Moller_RLET, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                else
                    call Gauss_Legendre_quadrature_1D &
                    (unitless_Moller_RLET, l_alpha, l_beta, xq1, wq1, integral)
                end if

                l_LET(g,mat) = mat_rho_e(mat)*twopi*(e_radius**2)*(e_mass_E**2)*&
                    integral/dEe(g)
            end do
        end do
    else if (electron_inelastic .eq. "RBED") then
        do mat = 1, Nmats
            do ele = 1, Nmaxeles
                if (eles_list(mat,ele) .eq. 0) cycle
                call read_atomic_parameters(eles_list(mat,ele), degen, W, U)
                Nshells = size(degen)

                do s = 1, Nshells
                    XSw = W(s)/e_mass_E
                    XSbeta_W2 = beta(W(s))**2
                    XSbeta_U2 = beta(U(s))**2
                    do g = 1, Ge-1
                        XScut = Ee(g+2)/e_mass_E
                        l_alpha = Ee(g+1)/e_mass_E
                        l_beta = Ee(g)/e_mass_E

                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_1D &
                            (unitless_RBED_RLET, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                        else
                            call Gauss_Legendre_quadrature_1D &
                            (unitless_RBED_RLET, l_alpha, l_beta, xq1, wq1, integral)
                        end if

                        l_LET(g,mat) = l_LET(g,mat) + &
                            mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius)**2*(e_mass_E**2)*&
                            integral/dEe(g)
                    end do
                    print *, "RBED CEPTRE RCSDA: do expression for g = Ge"
                    stop
                    !l_LET(Ge,mat) = l_LET(Ge,mat) +
                end do
            end do
        end do
    end if

    do g = 1, Ge
        Rt(g,:) = 3*l_LET(g,:)/(2*dEe(g))
    end do

    R(1,2,:) = 2*l_LET(1,:)/dEe(1)
    do g = 3, Ge
        R(g-1,g,:) = 2*l_LET(g-1,:)/dEe(g-1)
        R(g-2,g,:) = -0.5*l_LET(g-2,:)/dEe(g-2)
    end do
end subroutine MGXS_electron_electron_RCSDA

subroutine MGXS_electron_electron_elastic_fitted(Ee, dEe, adaptive, tol, N, S_el, St_el)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: N

    real, dimension(:,:,:), allocatable, intent(inout) :: S_el
    real, dimension(:,:), allocatable, intent(inout) :: St_el

    integer :: ele, g, i, j, jp, k, l, lp, mat
    integer :: l_NL
    integer :: nNint = 64
    character(50) :: Zstr
    character(150) :: fname
    real, dimension(:,:,:), allocatable :: integral
    real :: l_integral
    real, dimension(:,:,:,:), allocatable :: test_integrals
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:), allocatable :: Pil
    real*4, dimension(:,:), allocatable :: bcoeff ! k,j
    real :: mucut = 0.99
    real :: wt
    real :: nd

    integer :: NMT = 606
    integer :: NEs
    real, dimension(:,:), allocatable :: xs
    real, dimension(:), allocatable :: xsell0
    real, dimension(:), allocatable :: mu
    real, dimension(:), allocatable :: E
    real, dimension(5) :: term
    real :: l_beta
    real :: l_A
    real, dimension(:,:,:), allocatable :: B_S_el

    if (Eemin .lt. 1.0E-3 - 1.0E-6) then
        print *, "STOP!"
        print *, "MODULE MGXS, SUBROUTINE MGXS_electron_electron_elastic_fitted:"
        print *, "The electron cutoff energy is below 1 keV, which is the lowest energy"
        print *, "for which the fitting method is intended. You must either use"
        print *, "a larger cutoff energy, or use ELSEPA."
        print *, "PROGRAM ENDING."
    end if

    l_NL = NL
    if (elastic_etc) l_NL = NL + 1

    allocate(S_el(Ge,l_NL+1,Nmats))
    allocate(St_el(Ge,Nmats))
    allocate(bcoeff(6,5))
    allocate(integral(Ge,2*nNint,6))
    allocate(test_integrals(l_NL+1,Ge,5,6))
    allocate(Pil(2*nNint,0:l_NL))
    if (elastic_etc) allocate(B_S_el(Ge,NL+1,Nmats))

    call populate_Gauss_Legendre_quadrature(nNint, xq, wq)

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(N, xq1, wq1)
    end if

    Pil(1:nNint,0:l_NL) = Legendre_polynomials(l_NL,0.5*(mucut+1)*xq+0.5*(mucut-1))
    Pil(nNint+1:2*nNint,0:l_NL) = Legendre_polynomials(l_NL,0.5*(1-mucut)*xq+0.5*(1+mucut))

    S_el = 0
    St_el = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            write(Zstr,*) eles_list(mat,ele)

            fname = trim(adjustl(Zstr)) //".txt"

            open(1, file = "Physics data/Electron elastic/Fitted/"//fname)
            read(1,*) bcoeff
            close(1)

            XSZ = eles_list(mat,ele)

            do k = 1, 6
                XSk = k
                do i = 1, 2*nNint
                    if (i .le. nNint) then
                        XSmu = 0.5*(mucut+1)*xq(i) + 0.5*(mucut-1)
                    else
                        XSmu = 0.5*(1-mucut)*xq(i-nNint) + 0.5*(1+mucut)
                    end if
                    do g = 1, Ge
                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_1D &
                                (MGXS_electron_electron_elastic_fitted_integrand, &
                                beta(Ee(g+1)), beta(Ee(g)), xq1, wq1, xq2, wq2, tol, integral(g,i,k))
                        else
                            call Gauss_Legendre_quadrature_1D &
                                (MGXS_electron_electron_elastic_fitted_integrand, &
                                beta(Ee(g+1)), beta(Ee(g)), xq1, wq1, integral(g,i,k))
                        end if
                    end do
                end do
            end do

            do l = 0, l_NL
                do g = 1, Ge
                    do i = 1, 2*nNint
                        if (i .le. nNint) then
                            wt = 0.5*(mucut+1)*wq(i)
                            nd = 0.5*(mucut+1)*xq(i)+0.5*(mucut-1)
                        else
                            wt = 0.5*(1-mucut)*wq(i-nNint)
                            nd = 0.5*(1-mucut)*xq(i-nNint)+0.5*(1-mucut)
                        end if
                        do j = 0, 4
                            do k = 1, 6
                                S_el(g,l+1,mat) = S_el(g,l+1,mat) + &
                                    wt*mat_rho_a(mat,ele)*eles_list(mat,ele)**2*&
                                    e_mass_E*(e_radius**2)*0.5*(2*l+1)*((1-nd)**(0.5*j))*&
                                    Pil(i,l)*bcoeff(k,j+1)*integral(g,i,k)/dEe(g)
                                !S_el(g,l+1,mat) = S_el(g,l+1,mat) + &
                                !    wt*eles_list(mat,ele)**2*&
                                !    e_mass_E*(e_radius**2)*0.5*(2*l+1)*((1-nd)**(0.5*j))*&
                                !    Pil(i,l)*bcoeff(k,j+1)*integral(g,i,k)/dEe(g)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    if (elastic_etc) then
        do l = 0, NL
            B_S_el(:,l+1,:) = &
                S_el(:,l+1,:) - (2*l+1)*S_el(:,NL+2,:)/(2*NL+3)
        end do

        call move_alloc(B_S_el, S_el)
    end if

    St_el = fourpi*S_el(1:Ge,1,1:Nmats)

    !!allocate(E, source = [(0.5*(Ee(g)+Ee(g+1)),g=1,Ge)])
    !allocate(E, source = [0.001])
    !NEs = size(E)
    !allocate(xs(NMT,NEs))
    !allocate(mu(NMT))
    !allocate(xsell0(NEs))
    !
    !!mu = [(cos((j-1)*pi/(NMT-1)),j=1,NMT)]
    !mu = [(10**(-(i-1)/50.0),i=1,NMT)]
    !mu = 1-2*mu
    !xs = 0
    !xsell0 = 0

    !! log-log or trapz
    !do i = 1, NEs
    !    l_beta = beta(E(i))
    !    l_A = ((1-l_beta**2)/(l_beta**2))*&
    !    (1.13+3.76*(fsa*eles_list(1,1)/l_beta)**2)*(fsa*(eles_list(1,1))**(1.0/3)/1.77068)**2
    !    do j = 1, NMT
    !        do k = 1, 6
    !            do jp = 0, 4
    !                xs(j,i) = xs(j,i) + &
    !                ((e_radius*eles_list(1,1))**2)*((1-l_beta**2)/(l_beta**4))*&
    !                bcoeff(k,jp+1)*(1-mu(j))**(0.5*jp)*(l_beta-0.7181287)**(k-1)/&
    !                (2*l_A+1-mu(j))**2
    !            end do
    !        end do
    !    end do
    !end do
    !
    !do i = 1, NEs
    !    call trapezoidal_integration_1D &
    !    (0, xs(:,i), mu, -1.0, 1.0, xsell0(i))
    !    xsell0(i) = xsell0(i)/2
    !end do
    !
    !mu = (1 - mu)/2
    !
    !do i = 1, NEs
    !    call log_log_integration_1D &
    !    (0, xs(:,i), mu, 0.0, 1.0, xsell0(i))
    !    ! No factor of 2, because of integral transformation.
    !end do

    !! explicit
    !do i = 1, NEs
    !    l_beta = beta(E(i))
    !    l_A = ((1-l_beta**2)/(l_beta**2))*&
    !    (1.13+3.76*(fsa*eles_list(1,1)/l_beta)**2)*(fsa*(eles_list(1,1))**(1.0/3)/1.77068)**2
    !    term(1) = 1.0/(2*l_A*(l_A+1))
    !    term(2) = (l_A*atan(1/sqrt(l_A))-sqrt(l_A)+atan(1/sqrt(l_A)))/(sqrt(2*l_A)*(l_A+1))
    !    term(3) = log((2*l_A+2)/(2*l_A)) + l_A/(l_A+1) - 1
    !    term(4) = sqrt(2.0)*(-3*(l_A**(3.0/2))*atan(1/sqrt(l_A)) + 3*l_A - 3*sqrt(l_A)*atan(1/sqrt(l_A)) + 2)/&
    !        (l_A+1)
        !    term(5) = -(4*(l_A**2+l_A)*log((2*l_A+2)/(2*l_A)) - 4*l_A - 2)/(l_A+1)
    !
    !    do k = 1, 6
    !        do jp = 0, 4
    !            xsell0(i) = xsell0(i) + &
    !            0.5*((e_radius*eles_list(1,1))**2)*((1-l_beta**2)/(l_beta**4))*&
    !            bcoeff(k,jp+1)*(l_beta-0.7181287)**(k-1)*term(jp+1)
    !        end do
    !    end do
    !end do

    !! quadrature
    !
    !call populate_Gauss_Legendre_quadrature(64, xq1, wq1)
    !
    !do i = 1, NEs
    !    l_beta = beta(E(i))
    !    l_A = ((1-l_beta**2)/(l_beta**2))*&
    !    (1.13+3.76*(fsa*eles_list(1,1)/l_beta)**2)*(fsa*(eles_list(1,1))**(1.0/3)/1.77068)**2
    !    do j = 1, 64
    !        do k = 1, 6
    !            do jp = 0, 4
    !                xsell0(i) = xsell0(i) + &
    !                0.5*wq1(j)*((e_radius*eles_list(1,1))**2)*((1-l_beta**2)/(l_beta**4))*&
    !                bcoeff(k,jp+1)*(1-xq1(j))**(0.5*jp)*(l_beta-0.7181287)**(k-1)/&
    !                (2*l_A+1-xq1(j))**2
    !            end do
    !        end do
    !    end do
    !end do

    !!!!!!!!!!

    !do i = 1, NEs
    !    print *, xsell0(i)
    !end do

    !do j = 1, NMT
    !    print *, (1-mu(j))/2
    !end do
    !print *, "BREAK"
    !do j = 1, NMT
    !    print *, xs(j,1)
    !end do

    !stop

    !do l = 0, NL
    !    print *, "STARTING PRINTOUT:", l
    !    do g = 1, Ge
    !        print *, fourpi*Sigmaeeel(g,l+1,1)
    !    end do
    !end do

    !print *, "STARTING PRINTOUT:", 15
    !do g = 1, Ge
    !    print *, fourpi*Sigmaeeel(g,15+1,1)
    !end do
end subroutine MGXS_electron_electron_elastic_fitted

subroutine MGXS_electron_electron_elastic_Wentzel_Moliere(Ee, dEe, adaptive, tol, N, S_el, St_el)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: N

    real, dimension(:,:,:), allocatable, intent(inout) :: S_el
    real, dimension(:,:), allocatable, intent(inout) :: St_el

    integer :: ele, g, l, lp, mat
    integer :: l_NL
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:), allocatable :: coeff
    real :: integral
    real :: integral2
    real, dimension(:,:), allocatable :: integrals
    real, dimension(:,:,:), allocatable :: B_S_el

    l_NL = NL
    if (elastic_etc) l_NL = NL+1

    allocate(S_el(Ge,l_NL+1,Nmats))
    allocate(St_el(Ge,Nmats))
    allocate(coeff(0:l_NL,0:l_NL))
    allocate(integrals(Ge,l_NL+1))

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(N, xq1, wq1)
    end if

    coeff = Legendre_polynomial_coeffs(l_NL)

    S_el = 0
    St_el = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            XSZ = eles_list(mat,ele)

            do g = 1, Ge
                XSa = Ee(g+1)/e_mass_E
                XSb = Ee(g)/e_mass_E
                do l = 0, l_NL
                    XSlp = l
                    if (adaptive) then
                        call Gauss_Legendre_adaptive_quadrature_1D &
                            (MGXS_electron_electron_elastic_Wentzel_Moliere_integrand, &
                            -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                    else
                        call Gauss_Legendre_quadrature_1D &
                            (MGXS_electron_electron_elastic_Wentzel_Moliere_integrand, &
                            -1.0, 1.0, xq1, wq1, integral)
                    end if
                    integrals(g,l+1) = integral
                end do
            end do

            do g = 1, Ge
                do l = 0, l_NL
                    do lp = 0, l
                        S_el(g,l+1,mat) = S_el(g,l+1,mat) + &
                            mat_rho_a(mat,ele)*(e_radius**2)*(2*l+1)*&
                            coeff(l,lp)*e_mass_E*integrals(g,lp+1)/(2*dEe(g))
                    end do
                end do
            end do
        end do
    end do

    if (elastic_etc) then
        allocate(B_S_el(Ge,NL+1,Nmats))
        do l = 0, NL
            B_S_el(:,l+1,:) = &
                S_el(:,l+1,:) - (2*l+1)*S_el(:,NL+2,:)/(2*NL+3)
        end do

        call move_alloc(B_S_el, S_el)
    end if

    St_el = fourpi*S_el(:,1,:)
end subroutine MGXS_electron_electron_elastic_Wentzel_Moliere

subroutine MGXS_electron_electron_elastic_ELSEPA(Ee, dEe, S_el, St_el)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe

    real, dimension(:,:,:), allocatable, intent(inout) :: S_el
    real, dimension(:,:), allocatable, intent(inout) :: St_el

    integer :: ele, g, i, j, l, lp, mat, n
    integer :: l_NL
    integer :: NEs = 73
    integer :: NMT = 606
    real, dimension(:), allocatable :: E
    real, dimension(:), allocatable :: B_E
    integer :: start
    integer :: finish
    integer :: power
    integer :: second
    integer :: first
    character(50) :: trash1
    character(50) :: powerstr
    character(50) :: secondstr
    character(50) :: firststr
    character(50) :: Zstr
    character(150) :: fname
    real, dimension(5) :: trash5
    integer :: eof
    real, dimension(:), allocatable :: u
    real, dimension(:), allocatable :: xs
    real, dimension(:,:), allocatable :: coeff
    real, dimension(:,:), allocatable :: section
    real :: Mval
    real :: Bval
    real :: integral
    real :: x
    !real, dimension(:), allocatable :: xs_ell
    real, dimension(:,:), allocatable :: xs_ell
    real, dimension(:,:), allocatable :: xs_elli
    real, dimension(:,:,:), allocatable :: B_S_el

    l_NL = NL
    if (elastic_etc) l_NL = NL + 1

    allocate(S_el(Ge,l_NL+1,Nmats))
    allocate(St_el(Ge,Nmats))
    allocate(E(NEs))
    open(1, file = "Physics data/Electron elastic/ELSEPA/Evals.txt")
    do i = 1, NEs
        read(1,*) E(i)
    end do
    close(1)

    E = E/1E6

    start = findloc(Eemin .gt. E(1:NEs-1) .and. Eemin .le. E(2:NEs), .true., dim = 1)
    if (start .eq. 0) start = start + 1 ! Should I have to do this?
    finish = findloc(Eemax .gt. E(1:NEs-1) .and. Eemax .le. E(2:NEs), .true., dim = 1) + 1

    NEs = finish - start + 1

    allocate(B_E, source = E(start:finish))
    call move_alloc(B_E, E)
    allocate(B_E, source = E/e_mass_E)

    allocate(u(NMT))
    allocate(xs(NMT))
    allocate(coeff(0:l_NL,0:l_NL))
    allocate(section(l_NL+1,NMT-1))
    !allocate(xs_ell(l_NL+1))
    allocate(xs_ell(NMT,l_NL+1))
    allocate(xs_elli(NEs,l_NL+1))

    coeff = Legendre_polynomial_coeffs(l_NL)

    ! NOTE: If ever doing molecules, don't use mat_rho_a, use mat_rho_molecule or something like that.

    S_el = 0
    St_el = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            write(Zstr,*) eles_list(mat,ele)

            do i = 1, NEs
                power = floor(log10(E(i)))
                second = nint(mod(E(i)*10.0**(-power),1.0)*1000)
                if (second .eq. 1000) then ! SHOULDN'T NEED TO DO THIS, FIX!!
                    second = 0
                end if
                first = nint(E(i)*10.0**(-power) - real(second)/1000)

                write(powerstr,*) power + 6
                write(firststr,*) first

                if (second .eq. 0) then ! CHANGE THESE
                    fname = trim(adjustl(Zstr)) //"/dcs_"&
                    // trim(adjustl(firststr)) //"p000e0" &
                    // trim(adjustl(powerstr)) //".dat"
                else
                    write(secondstr,*) second

                    fname = trim(adjustl(Zstr)) //"/dcs_"&
                    // trim(adjustl(firststr)) //"p"// trim(adjustl(secondstr)) // &
                    "e0" // trim(adjustl(powerstr)) //".dat"
                end if

                open(1, file = "Physics data/Electron elastic/ELSEPA/"//fname)
                trashloop: do
                    read(1,*) trash1
                    if (trash1 .ne. "#") then
                        exit trashloop
                    end if
                end do trashloop
                xs = 0
                do j = 1, NMT
                    read(1,*, iostat = eof) trash5
                    u(j) = trash5(2)
                    xs(j) = trash5(3)
                    !print *, xs(j)
                    if (eof .lt. 0) exit
                end do
                close(1)

                !if (abs(E(i) - 0.001) .le. 0.000001) then
                !    do j = 1, NMT
                !        print *, u(j)
                !    end do
                !    print *, "BREAK"
                !    do j = 1, NMT
                !        print *, xs(j)
                !    end do
                !    print *, "BREAK"
                !    print *, "ENERGY:", E(i)
                !    stop
                !end if

                call log_log_integration_1D &
                    (0, xs, u, 0.0, 1.0, xs_elli(i,1))

                do l = 1, l_NL
                    integral = 0
                    Mval = (xs(2)-xs(1))/(u(2)-u(1))
                    Bval = xs(1) - Mval*u(1)
                    do lp = 0, l
                        do n = 0, lp
                            integral = integral + binomial_int(lp,n)*coeff(l,lp)*(-1)**(n)*2**(n)*&
                                (Mval*(u(2)**(n+2)-u(1)**(n+2))/(n+2) + &
                                Bval*(u(2)**(n+1)-u(1)**(n+1))/(n+1))
                        end do
                    end do
                    do j = 2, NMT-1
                        Mval = log10(xs(j+1)/xs(j))/log10(u(j+1)/u(j))
                        Bval = log10(xs(j)/(u(j)**Mval))
                        do lp = 0, l
                            do n = 0, lp
                                integral = integral + 10**(Bval)*binomial_int(lp,n)*coeff(l,lp)*(-1)**(n)*2**(n)*&
                                    (u(j+1)**(Mval+n+1)-u(j)**(Mval+n+1))/(Mval+n+1)
                            end do
                        end do
                    end do

                    xs_elli(i,l+1) = (2*l+1)*integral
                end do

                !if (i .eq. 21) then
                !    print *, E(i), mat_rho_a(mat,ele)*5.14403E-17, fourpi*mat_rho_a(mat,ele)*xs_elli(i,1)
                !    stop
                !end if

                !print *, E(i), fourpi*mat_rho_a(mat,ele)*xs_elli(i,1)

                !print *, &
                ![E(i),fourpi*xs_elli(i,1), fourpi*(xs_elli(i,1)-xs_elli(i,2)/3),fourpi*(xs_elli(i,1)-xs_elli(i,3)/5)]
            end do
            !print *, "BREAK"

            !do i = 1, NEs
            !    print *, E(i)
            !end do
            !print *, "BREAK"
            !do i = 1, NEs
            !    print *, fourpi*xs_elli(i,1)
            !end do
            !print *, "BREAK"
            !do i = 1, NEs
            !    print *, fourpi*(xs_elli(i,1)-xs_elli(i,2)/3)
            !end do
            !print *, "BREAK"
            !do i = 1, NEs
            !    print *, fourpi*(xs_elli(i,1)-xs_elli(i,3)/5)
            !end do
            !stop

            !do i = 1, size(xs_elli,1)
            !    print *, E(i)
            !end do
            !print *, "BREAK"
            !do i = 1, size(xs_elli,1)
            !    print *, xs_elli(i,1)
            !end do

            do l = 0, l_NL
                if (any(xs_elli(1:NEs,l+1) .le. 0)) then
                    do g = 1, Ge
                        call trapezoidal_integration_1D & ! INVESTIGATE NEGATIVES????
                            (0, xs_elli(1:NEs,l+1), B_E, Ee(g+1)/e_mass_E, Ee(g)/e_mass_E, integral)
                        S_el(g,l+1,mat) = S_el(g,l+1,mat) + &
                            mat_rho_a(mat,ele)*e_mass_E*integral/dEe(g)
                        !S_el(g,l+1,mat) = S_el(g,l+1,mat) + &
                        !    e_mass_E*integral/dEe(g)
                    end do
                else
                    do g = 1, Ge
                        call log_log_integration_1D &
                            (0, xs_elli(1:NEs,l+1), B_E, Ee(g+1)/e_mass_E, Ee(g)/e_mass_E, integral)
                        S_el(g,l+1,mat) = S_el(g,l+1,mat) + &
                            mat_rho_a(mat,ele)*e_mass_E*integral/dEe(g)
                        !S_el(g,l+1,mat) = S_el(g,l+1,mat) + &
                        !    e_mass_E*integral/dEe(g)
                    end do
                end if
            end do
        end do
    end do

    if (elastic_etc) then
        allocate(B_S_el(Ge,NL+1,Nmats))
        do l = 0, NL
            B_S_el(:,l+1,:) = &
                S_el(:,l+1,:) - (2*l+1)*S_el(:,NL+2,:)/(2*NL+3)
        end do

        call move_alloc(B_S_el, S_el)
    end if

    St_el = fourpi*S_el(1:Ge,1,1:Nmats)
end subroutine MGXS_electron_electron_elastic_ELSEPA

subroutine MGXS_electron_electron_Auger(Ee, dEe, adaptive, tol, N, S_A)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: N

    real, dimension(:,:,:), allocatable, intent(inout) :: S_A

    integer :: cascade, ele, g, gpr, i, j, mat, s
    integer :: Nshells
    integer :: Ncascades
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W ! MeV
    real, dimension(:), allocatable :: U ! MeV
    integer, dimension(:), allocatable :: cscpershell
    real, dimension(:), allocatable :: Eij ! MeV
    real, dimension(:), allocatable :: eta
    integer, dimension(:,:), allocatable :: deltas
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    integer :: globalcascade
    real :: l_alpha
    real :: l_beta
    real :: integral
    real, dimension(:), allocatable :: integrals

    allocate(S_A(Ge,Ge,Nmats))
    allocate(integrals(Ge))

    S_A = 0
    integrals = 0

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(N, xq1, wq1)
    end if

    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            if (eles_list(mat,ele) .lt. 6) cycle

            ! READ PARAMETERS
            call read_atomic_parameters(eles_list(mat,ele), degen, W, U)
            Nshells = size(degen)

            ! READ AUGER DATA
            call read_relaxation_data("auger", eles_list(mat,ele), Nshells, Eij, eta, cscpershell)
            Ncascades = size(Eij)

            if (all(Eij .lt. Eemin)) go to 1

            allocate(deltas(Ge,Ncascades))
            do cascade = 1, Ncascades
                do g = 1, Ge
                    if (Ee(g+1) .le. Eij(cascade) .and. Eij(cascade) .le. Ee(g)) then
                        deltas(g,cascade) = 1
                    else
                        deltas(g,cascade) = 0
                    end if
                end do
            end do

            globalcascade = 0
            do s = 1, Nshells
                XSbeta_U2 = beta(U(s))**2
                XSbeta_W2 = beta(W(s))**2
                XSw = W(s)/e_mass_E
                integrals = 0

                do gpr = 1, Ge
                    if (electron_exc_shells) then
                        if (Ee(gpr+1) + 1.0E-6 .lt. W(s)) cycle
                        l_alpha = Ee(gpr+1)/e_mass_E
                        l_beta = Ee(gpr)/e_mass_E
                    else
                        l_alpha = min(max(Ee(gpr+1),W(s)),Ee(gpr))/e_mass_E
                        l_beta = Ee(gpr)/e_mass_E
                    end if
                    if (electron_inelastic .eq. "Moller") then
                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_1D &
                                (unitless_Gryzinski, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                        else
                            call Gauss_Legendre_quadrature_1D &
                                (unitless_Gryzinski, l_alpha, l_beta, xq1, wq1, integral)
                        end if
                    else if (electron_inelastic .eq. "RBED") then
                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_1D &
                                (unitless_RBEB, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                        else
                            call Gauss_Legendre_quadrature_1D &
                                (unitless_RBEB, l_alpha, l_beta, xq1, wq1, integral)
                        end if
                    end if
                    integrals(gpr) = twopi*(e_radius**2)*degen(s)*e_mass_E*&
                        integral
                end do

                do cascade = 1, cscpershell(s)
                    globalcascade = globalcascade + 1
                    do g = 1, Ge
                        do gpr = 1, Ge
                            S_A(gpr,g,mat) = S_A(gpr,g,mat) + &
                                mat_rho_a(mat,ele)*deltas(g,globalcascade)*eta(globalcascade)*&
                                integrals(gpr)/(fourpi*dEe(gpr))
                        end do
                    end do
                end do
            end do

            deallocate(deltas)
            1 continue
            deallocate(degen)
            deallocate(W)
            deallocate(U)
            deallocate(cscpershell)
            deallocate(Eij)
            deallocate(eta)
        end do
    end do
end subroutine MGXS_electron_electron_Auger

subroutine MGXS_electron_electron_Bremsstrahlung(Ee, dEe, adaptive, tol, N, S_B, St_B)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: N

    real, dimension(:,:,:), allocatable, intent(inout) :: S_B ! Treated like RCSDA, in that it's just the non angular part of the DCS
    real, dimension(:,:), allocatable, intent(inout) :: St_B

    integer :: ele, g, gpr, mat
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    integer :: Z
    real, dimension(:), allocatable :: trash1
    real, dimension(:), allocatable :: trash2
    real :: integral

    allocate(S_B(Ge,Ge,Nmats))
    allocate(St_B(Ge,Nmats))

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(N, xq1, wq1)
    end if

    S_B = 0
    St_B = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            Z = eles_list(mat,ele)
            call read_pair_production(Z, 2.0, trash1, trash2) ! Sets XSR and XSeta. Also, 2.0 is meaningless here.

            do gpr = 1, Ge-1
                XSa = Ee(gpr+1)/e_mass_E
                XSb = Ee(gpr)/e_mass_E
                XScut = Ee(gpr+2)/e_mass_E
                do g = gpr+2, Ge
                    XSc = Ee(g+1)/e_mass_E
                    XSd = Ee(g)/e_mass_E
                    if (adaptive) then
                        call Gauss_Legendre_adaptive_quadrature_2D &
                        (MGXS_electron_electron_Bremsstrahlung_integrand, &
                        -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                    else
                        call Gauss_Legendre_quadrature_2D &
                        (MGXS_electron_electron_Bremsstrahlung_integrand, &
                        -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral)
                    end if

                    S_B(gpr,g,mat) = S_B(gpr,g,mat) + &
                        mat_rho_a(mat,ele)*e_mass_E*twopi*(e_radius**2)*fsa*Z*(Z+XSeta)*&
                        integral/dEe(gpr)
                end do

                if (adaptive) then
                    call Gauss_Legendre_adaptive_quadrature_2D &
                    (MGXS_electron_electron_Bremsstrahlung_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                else
                    call Gauss_Legendre_quadrature_2D &
                    (MGXS_electron_electron_Bremsstrahlung_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral)
                end if

                St_B(gpr,mat) = St_B(gpr,mat) + &
                    mat_rho_a(mat,ele)*e_mass_E*twopi*(e_radius**2)*fsa*Z*(Z+XSeta)*&
                    integral/dEe(gpr)
            end do

            deallocate(trash1)
            deallocate(trash2)
        end do
    end do
end subroutine MGXS_electron_electron_Bremsstrahlung

!--------------------- Electron-Photon ------------------------------------------------------------



! TOTAL ATTENUATION
!--------------------- Photon ---------------------------------------------------------------------
subroutine MGXS_photon_attn(Ep, dEp, adaptive, tol, Nint, Sigmapt, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:), allocatable, intent(inout) :: Sigmapt
    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: ele, g, gpr, mat
    integer :: Gpp
    integer :: NNIST
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:), allocatable :: Sigmapct
    real, dimension(:,:), allocatable :: Sigmaptpp
    real, dimension(:,:), allocatable :: Sigmaptpe
    real :: integral
    real, dimension(:), allocatable :: ENist
    real, dimension(:), allocatable :: ppxs
    real, dimension(:), allocatable :: peWvals
    real, dimension(:), allocatable :: Espec
    real, dimension(:,:), allocatable :: pexs
    real, dimension(:), allocatable :: pexst
    integer :: NEspec

    allocate(Sigmapt(Gp,Nmats))
    allocate(Sigmapct(Gp))
    allocate(Sigmaptpp(Gp,Nmats))
    allocate(Sigmaptpe(Gp,Nmats))

    if (Epmax .ge. 1.022) then
        Gpp = findloc(Ep, 1.022, dim = 1)
        if (Gpp .eq. 0) then
            print *, "subroutine ----: Photon energy grouping does not contain 1.022 MeV."
            stop
        end if
        Gpp = Gpp - 1
    else
        Gpp = 0
    end if

    Sigmapct = 0
    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
        do g = 1, Gp
            call Gauss_Legendre_adaptive_quadrature_1D &
            (unitless_sigmapct, Ep(g+1)/e_mass_E, Ep(g)/e_mass_E, xq1, wq1, xq2, wq2, tol, integral) ! Breaks from convention because bounds are simple.
            Sigmapct(g) = twopi*(e_radius**2)*e_mass_E*integral/dEp(g)
        end do
    else
        print *, "MGXS_photon_attn: non adaptive is WIP"
    end if

    Sigmaptpp = 0
    Sigmaptpe = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            if (Epmax .ge. 1.022) then
                call read_pair_production(eles_list(mat,ele), Epmax, ENist, ppxs)
                NNist = size(ENist)

                do g = 1, Gpp
                    call log_log_integration_1D &
                        (0, ppxs, ENist, Ep(g+1)/e_mass_E, Ep(g)/e_mass_E, integral)
                    Sigmaptpp(g,mat) = Sigmaptpp(g,mat) + mat_rho(mat)*mass_frac(mat,ele)*e_mass_E*&
                        integral/dEp(g)
                end do
            end if

            call read_photoelectric_effect(eles_list(mat,ele), Epmin, Epmax, peWvals, Espec, pexs)
            NEspec = size(Espec)

            allocate(pexst(NEspec))

            pexst = pexs(:,size(peWvals)+1)

            do g = 1, Gp
                call log_log_integration_1D &
                    (0, pexst, Espec, Ep(g+1)/e_mass_E, Ep(g)/e_mass_E, integral)
                Sigmaptpe(g,mat) = Sigmaptpe(g,mat) + mat_rho_a(mat,ele)*e_mass_E*integral/dEp(g)
            end do

            if (Epmax .ge. 1.022) then
                deallocate(ENist)
                deallocate(ppxs)
            end if
            deallocate(peWvals)
            deallocate(Espec)
            deallocate(pexs)
            deallocate(pexst)
        end do
    end do

    do mat = 1, Nmats
        Sigmapt(1:Gp,mat) = mat_rho_e(mat)*Sigmapct + Sigmaptpp(1:Gp,mat) + Sigmaptpe(1:Gp,mat)
    end do

    do g = 1, Gp
        EDEP(g,:) = EDEP(g,:) + (mat_rho_e*Sigmapct(g) + Sigmaptpe(g,:))*0.5*(Ep(g)+Ep(g+1))
    end do

    do g = 1, Gp
        EDEP(g,:) = EDEP(g,:) + Sigmaptpp(g,:)*0.5*(Ep(g)+Ep(g+1))
    end do
end subroutine MGXS_photon_attn

!--------------------- Electron -------------------------------------------------------------------
subroutine MGXS_electron_attn(Ee, dEe, p_Sigmaet, adaptive, tol, Nint, Sigmaet, Sigmaea, EDEP, CDEP)
    ! NOT USED ANYMORE. REMOVE??
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:,:), intent(in) :: p_Sigmaet
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:), allocatable, intent(inout) :: Sigmaet
    real, dimension(:,:), allocatable, intent(inout) :: Sigmaea
    real, dimension(:,:), allocatable, intent(inout) :: EDEP
    real, dimension(:,:), allocatable, intent(inout) :: CDEP

    integer :: ele, g, mat, s
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:), allocatable :: Sigmaemt
    real, dimension(:,:), allocatable :: SigmaeBt
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    integer :: Z
    real, dimension(:), allocatable :: trash1
    real, dimension(:), allocatable :: trash2
    integer :: Nshells
    real :: integral
    real :: l_alpha
    real :: l_beta
    real :: x

    if (ccut .eq. 0) XScut = 0

    allocate(Sigmaet(Ge,Nmats))
    allocate(Sigmaemt(Ge,Nmats))
    allocate(SigmaeBt(Ge,Nmats))

    Sigmaemt = 0
    SigmaeBt = 0

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(Nint, xq1, wq1)
    end if

    if (electron_inelastic .eq. "Moller") then
        XSw = 0

        do g = 1, Ge-1
            XSa = Ee(g+1)/e_mass_E
            XSb = Ee(g)/e_mass_E
            XScut = Ee(g+2)/e_mass_E
            if (adaptive) then
                call Gauss_Legendre_adaptive_quadrature_1D &
                (MGXS_electron_Moller_R_total_integrand, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
            else
                call Gauss_Legendre_quadrature_1D &
                (MGXS_electron_Moller_R_total_integrand, -1.0, 1.0, xq1, wq1, integral)
            end if

            do mat = 1, Nmats
                Sigmaemt(g,mat) = mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*integral/dEe(g)
            end do
        end do
    else if (electron_inelastic .eq. "RBED") then
        do mat = 1, Nmats
            do ele = 1, Nmaxeles
                if (eles_list(mat,ele) .eq. 0) cycle
                call read_atomic_parameters(eles_list(mat,ele), degen, W, U)
                Nshells = size(degen)

                integral = 0
                do s = 1, Nshells
                    XSw = W(s)/e_mass_E
                    XSbeta_U2 = beta(U(s))**2
                    XSbeta_W2 = beta(W(s))**2
                    do g = 1, merge(Ge-1, Ge, ccut .gt. 0)
                        XSa = Ee(g+1)/e_mass_E
                        XSb = Ee(g)/e_mass_E
                        if (ccut .gt. 0) XScut = Ee(g+2)/e_mass_E
                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_1D &
                            (MGXS_electron_RBED_total_integrand, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                        else
                            call Gauss_Legendre_quadrature_1D &
                            (MGXS_electron_RBED_total_integrand, -1.0, 1.0, xq1, wq1, integral)
                        end if

                        Sigmaemt(g,mat) = Sigmaemt(g,mat) + &
                            mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius**2)*e_mass_E*integral/dEe(g)
                    end do
                end do
                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do
        end do
    end if

    SigmaeBt = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            Z = eles_list(mat,ele)
            call read_pair_production(Z, 2.0, trash1, trash2) ! Sets XSR and XSeta

            do g = 1, Ge-1
                XSa = Ee(g+1)/e_mass_E
                XSb = Ee(g)/e_mass_E
                XScut = Ee(g+2)/e_mass_E
                if (adaptive) then
                    call Gauss_Legendre_adaptive_quadrature_2D &
                    (MGXS_electron_electron_Bremsstrahlung_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                else
                    call Gauss_Legendre_quadrature_2D &
                    (MGXS_electron_electron_Bremsstrahlung_integrand, &
                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral)
                end if

                SigmaeBt(g,mat) = SigmaeBt(g,mat) + &
                    mat_rho_a(mat,ele)*e_mass_E*twopi*(e_radius**2)*fsa*Z*(Z+XSeta)*&
                    integral/dEe(g)
            end do

            deallocate(trash1)
            deallocate(trash2)
        end do
    end do

    Sigmaet = Sigmaemt + SigmaeBt + p_Sigmaet

    ! DEPOSITION CROSS SECTIONS
    do g = 1, Ge
        Sigmaea(g,:) = Sigmaea(g,:) + (Sigmaemt(g,:) + SigmaeBt(g,:))
        EDEP(g,:) = EDEP(g,:) + (Sigmaemt(g,:) + SigmaeBt(g,:))*0.5*(Ee(g)+Ee(g+1))
        CDEP(g,:) = CDEP(g,:) + (Sigmaemt(g,:) + SigmaeBt(g,:))
    end do
end subroutine MGXS_electron_attn

! STOPPING POWER
subroutine MGXS_stopping_power(Ee, dEe, adaptive, tol, Nint, EDEP)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:), allocatable, intent(inout) :: EDEP

    integer :: ele, g, mat, s
    real :: Ibar
    real :: integral
    real :: l_alpha
    real :: l_beta
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    integer :: Nshells
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U

    allocate(EDEP(Ge,Nmats))

    EDEP = 0

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(Nint, xq1, wq1)
    end if

    if (electron_inelastic .eq. "Moller") then
        XSw = 0

        do mat = 1, Nmats
            call mean_excitation_energy(eles_list(mat,:), Natoms(mat,:), Ibar)
            call set_density_effect_parameters &
                (mat_rho(mat), total_GAM(mat), eles_list(mat,:), Natoms(mat,:), Ibar, .false.)

            XSI = Ibar*1.0E-6/e_mass_E
            do g = 1, Ge
                l_alpha = Ee(g+1)/e_mass_E
                l_beta = Ee(g)/e_mass_E
                if (adaptive) then
                    call Gauss_Legendre_adaptive_quadrature_1D &
                        (unitless_Bethe, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                else
                    call Gauss_Legendre_quadrature_1D &
                        (unitless_Bethe, l_alpha, l_beta, xq1, wq1, integral)
                end if
                EDEP(g,mat) = &
                    mat_rho_e(mat)*twopi*(e_radius**2)*(e_mass_E**2)*integral/(mat_rho(mat)*dEe(g))
            end do
        end do

        !do g = 1, Ge
        !    print *, 0.5*(Ee(g)+Ee(g+1))
        !end do
        !print *, "BREAK"
        !do g = 1, Ge
        !    print *, EDEP(g,1)
        !end do
        !stop
    else if (electron_inelastic .eq. "RBED") then
        XSlogic = .true.
        XScut = 0 ! Set because we want full LET, not any kind of restricted LET. Doesn't
                  ! need to be set for Moller because unitless_Bethe does not reference XScut at all.

        do mat = 1, Nmats
            do ele = 1, Nmaxeles
                if (eles_list(mat,ele) .eq. 0) cycle
                call read_atomic_parameters(eles_list(mat,ele), degen, W, U)

                Nshells = size(degen)
                do s = 1, Nshells
                    XSw = W(s)/e_mass_E
                    XSbeta_U2 = beta(U(s))**2
                    XSbeta_W2 = beta(W(s))**2
                    do g = 1, Ge
                        l_alpha = Ee(g+1)/e_mass_E
                        l_beta = Ee(g)/e_mass_E
                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_1D &
                                (unitless_RBED_LET, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                        else
                            call Gauss_Legendre_quadrature_1D &
                                (unitless_RBED_LET, l_alpha, l_beta, xq1, wq1, integral)
                        end if
                        EDEP(g,mat) = EDEP(g,mat) + &
                            mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius**2)*(e_mass_E**2)*&
                            integral/(mat_rho(mat)*dEe(g))
                    end do
                end do

                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do
        end do
        XSlogic = .false.
    end if
end subroutine MGXS_stopping_power

end module MGXS