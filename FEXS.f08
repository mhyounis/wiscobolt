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

module FEXS
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

contains
subroutine FEXS_construct_cross_sections(particle)
    implicit none
    integer, intent(in) :: particle

    !if (particle .eq. 1) then
    !    call FEXS_photon_photon()
    !else if (particle .eq. 12) then
    !    call FEXS_photon_electron()
    !else if (particle .eq. 2) then
    !    call FEXS_electron_electron()
    !end if
end subroutine FEXS_construct_cross_sections

! SCATTERING
!--------------------- Photon-Photon --------------------------------------------------------------
subroutine FEXS_photon_photon(Ep, dEp, adaptive, tol, Nint, Sigmapp)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Sigmapp

    integer :: g, gpr, l, mat, n, np
    real, dimension(:,:,:,:,:), allocatable :: Sigmappc
    real, dimension(:,:,:,:,:), allocatable :: SigmappF
    real, dimension(:,:,:,:,:,:), allocatable :: B_Sigma

    call FEXS_photon_photon_Compton(Ep, dEp, adaptive, tol, Nint, Sigmappc)
    !call FEXS_photon_photon_fluorescence(Ep, dEp, SigmappF)
    !call FEXS_photon_photon_Rayleigh

    allocate(Sigmapp(2,2,Gp,Gp,NL+1,Nmats))
    Sigmapp = 0
    !Sigmapp(1:2,1:2,1:Gp,1:Gp,1,1:Nmats) = SigmappF
    do mat = 1, Nmats
        Sigmapp(1:2,1:2,1:Gp,1:Gp,1:NL+1,mat) = Sigmapp(1:2,1:2,1:Gp,1:Gp,1:NL+1,mat) + &
            mat_rho_e(mat)*Sigmappc
    end do

    allocate(B_Sigma(2,2,Gp,Gp,NL+1,Nmats))
    do g = 1, Gp
        B_Sigma(1,1,:,g,:,:) = Sigmapp(1,1,:,g,:,:)*4/dEp(g) - Sigmapp(1,2,:,g,:,:)*2/dEp(g)
        B_Sigma(2,1,:,g,:,:) = Sigmapp(2,1,:,g,:,:)*4/dEp(g) - Sigmapp(2,2,:,g,:,:)*2/dEp(g)
        B_Sigma(1,2,:,g,:,:) = -Sigmapp(1,1,:,g,:,:)*2/dEp(g) + Sigmapp(1,2,:,g,:,:)*4/dEp(g)
        B_Sigma(2,2,:,g,:,:) = -Sigmapp(2,1,:,g,:,:)*2/dEp(g) + Sigmapp(2,2,:,g,:,:)*4/dEp(g)
    end do

    Sigmapp = B_Sigma
end subroutine FEXS_photon_photon

subroutine FEXS_photon_photon_Compton(Ep, dEp, adaptive, tol, Nint, Sigmappc)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: Sigmappc

    integer :: g, gpr, l, lp, n, np
    real, dimension(:,:), allocatable :: coeff
    real :: l_integral
    real, dimension(:,:,:,:,:), allocatable :: integral
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:), allocatable :: portion
    real, dimension(:,:,:,:,:), allocatable :: B_Sigmappc

    allocate(Sigmappc(2,2,Gp,Gp,NL+1))
    allocate(coeff(0:NL,0:NL))
    allocate(integral(2,2,Gp,Gp,NL+1))
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
            do np = 1, 2
                XSEgp = Ep(gpr-np+2)/e_mass_E
                do n = 1, 2
                    XSEg = Ep(g-n+2)/e_mass_E
                    do lp = 0, NL
                        XSlp = lp
                        if (adaptive) then
                            call Gauss_Legendre_adaptive_quadrature_2D &
                                (FEXS_photon_photon_Compton_integrand, &
                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, l_integral)
                        else
                            call Gauss_Legendre_quadrature_2D &
                                (FEXS_photon_photon_Compton_integrand, &
                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, l_integral)
                        end if
                        integral(np,n,gpr,g,lp+1) = l_integral
                    end do
                end do
            end do
        end do
    end do

    Sigmappc = 0
    do g = 1, Gp
        do gpr = 1, g
            do n = 1, 2
                do np = 1, 2
                    do l = 0, NL
                        do lp = 0, l
                            Sigmappc(np,n,gpr,g,l+1) = Sigmappc(np,n,gpr,g,l+1) + &
                                ((-1)**(n+np))*(e_radius**2)*(2*l+1)*(e_mass_E**3)*&
                                coeff(l,lp)*integral(np,n,gpr,g,lp+1)/(4*dEp(gpr)*dEp(g))
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine FEXS_photon_photon_Compton

subroutine FEXS_photon_photon_Rayleigh()
    implicit none
end subroutine FEXS_photon_photon_Rayleigh

subroutine FEXS_photon_photon_fluorescence(Ep, dEp, SigmappF)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: SigmappF

    integer :: cascade, ele, g, gpr, i, j, mat, n, np, shell
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
    integer, dimension(:,:,:), allocatable :: Lambdas
    integer :: globalcascade
    real :: l_alpha
    real :: l_beta
    real :: integral
    real, dimension(:,:), allocatable :: integrals

    allocate(SigmappF(2,2,Gp,Gp,Nmats))
    allocate(integrals(2,Gp))

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
                print *, "MODULE FEXS, SUBROUTINE FEXS_photon_photon_Fluorescence:"
                print *, "Number of cascades available for photoionization doesn't match"
                print *, "the number of energies or efficiencies after truncating from EADL data."
                print *, "PROGRAM ENDING."
                stop
            end if

            if (all(Eij .lt. Eemin)) go to 1

            allocate(deltas(Gp,Ncascades))
            allocate(Lambdas(2,Gp,Ncascades))
            do cascade = 1, Ncascades
                do g = 1, Gp
                    if (Ep(g+1) .le. Eij(cascade) .and. Eij(cascade) .le. Ep(g)) then
                        deltas(g,cascade) = 1
                    else
                        deltas(g,cascade) = 0
                    end if
                    do n = 1, 2
                        Lambdas(n,g,cascade) = ((-1)**n)*(Ep(g-n+2)-Eij(cascade))/dEp(g)
                    end do
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
                    do np = 1, 2
                        MATHA = Ep(gpr-np+2)/e_mass_E
                        call log_log_integration_1D &
                            (2, pexs(:,shell), Espec, l_alpha, l_beta, integrals(np,gpr))
                    end do
                end do

                do cascade = 1, cscpershell(shell) ! Could double check this does what I want it to, ~4/22/23
                    globalcascade = globalcascade + 1
                    do g = 1, Gp
                        do gpr = 1, Gp
                            do np = 1, 2
                                do n = 1, 2
                                    SigmappF(np,n,gpr,g,mat) = SigmappF(np,n,gpr,g,mat) + ((-1)**np)*&
                                        mat_rho_a(mat,ele)*deltas(g,globalcascade)*&
                                        Lambdas(n,g,globalcascade)*eta(globalcascade)*&
                                        (e_mass_E**2)*integrals(np,gpr)/(fourpi*dEp(gpr))
                                end do
                            end do
                        end do
                    end do
                end do
            end do

            deallocate(deltas)
    1       deallocate(degen)
            deallocate(Wvals)
            deallocate(Espec)
            deallocate(pexs)
            deallocate(Lambdas)
            deallocate(cscpershell)
            deallocate(Eij)
            deallocate(eta)
            deallocate(matches)
        end do
    end do
end subroutine FEXS_photon_photon_fluorescence

!--------------------- Photon-Electron ------------------------------------------------------------
subroutine FEXS_photon_electron(Ep, dEp, Ee, dEe, Eemid, adaptive, tol, Nint, Sigmape)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Sigmape

    integer :: g, mat
    real, dimension(:,:,:,:,:), allocatable :: Sigmapec
    real, dimension(:,:,:,:,:,:), allocatable :: Sigmapepp
    real, dimension(:,:,:,:,:,:), allocatable :: Sigmapepe
    real, dimension(:,:,:,:,:), allocatable :: SigmapeA
    real, dimension(:,:,:,:,:,:), allocatable :: B_Sigma

    call FEXS_photon_electron_Compton(Ep, dEp, Ee, dEe, adaptive, tol, Nint, Sigmapec)
    if (Epmax .le. 1.022) then
        allocate(Sigmapepp(2,2,Gp,Ge,NL+1,Nmats))
        Sigmapepp = 0
    else
        call FEXS_photon_electron_pair_production(Ep, dEp, Ee, dEe, Eemid, Nint, Sigmapepp)
    end if
    call FEXS_photon_electron_photoelectric_effect(Ep, dEp, Ee, dEe, Eemid, Sigmapepe)
    !call FEXS_photon_electron_Auger(Ep, dEp, Ee, dEe, SigmapeA)

    allocate(Sigmape(2,2,Gp,Ge,NL+1,Nmats))
    Sigmape = 0
    !Sigmape(1:2,1:2,1:Gp,1:Ge,1,1:Nmats) = SigmapeA
    do mat = 1, Nmats
        Sigmape(:,:,:,:,:,mat) = Sigmape(:,:,:,:,:,mat) + mat_rho_e(mat)*Sigmapec
    end do
    Sigmape = Sigmape + Sigmapepp + Sigmapepe

    allocate(B_Sigma(2,2,Gp,Ge,NL+1,Nmats))
    do g = 1, Ge
        B_Sigma(1,1,:,g,:,:) = Sigmape(1,1,:,g,:,:)*4/dEe(g) - Sigmape(1,2,:,g,:,:)*2/dEe(g)
        B_Sigma(2,1,:,g,:,:) = Sigmape(2,1,:,g,:,:)*4/dEe(g) - Sigmape(2,2,:,g,:,:)*2/dEe(g)
        B_Sigma(1,2,:,g,:,:) = -Sigmape(1,1,:,g,:,:)*2/dEe(g) + Sigmape(1,2,:,g,:,:)*4/dEe(g)
        B_Sigma(2,2,:,g,:,:) = -Sigmape(2,1,:,g,:,:)*2/dEe(g) + Sigmape(2,2,:,g,:,:)*4/dEe(g)
    end do

    Sigmape = B_Sigma
end subroutine FEXS_photon_electron

subroutine FEXS_photon_electron_Compton(Ep, dEp, Ee, dEe, adaptive, tol, Nint, Sigmapec)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: Sigmapec

    integer :: g, gpr, l, lp, n, np
    real, dimension(:,:), allocatable :: coeff
    real, dimension(:,:,:,:,:), allocatable :: integral
    real :: l_integral
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:), allocatable :: portion

    allocate(Sigmapec(2,2,Gp,Ge,NL+1))
    allocate(coeff(0:NL,0:NL))
    allocate(integral(2,2,Gp,Ge,NL+1))
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
                do np = 1, 2
                    XSEgp = Ep(gpr-np+2)/e_mass_E
                    do n = 1, 2
                        XSEg = Ee(g-n+2)/e_mass_E
                        do lp = 0, NL
                            XSlp = lp
                            call Gauss_Legendre_adaptive_quadrature_2D &
                                (FEXS_photon_electron_Compton_integrand, &
                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, l_integral)
                            integral(np,n,gpr,g,lp+1) = l_integral
                        end do
                    end do
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
                do np = 1, 2
                    XSEgp = Ep(gpr-np+2)/e_mass_E
                    do n = 1, 2
                        XSEg = Ee(g-n+2)/e_mass_E
                        do lp = 0, NL
                            XSlp = lp
                            call Gauss_Legendre_quadrature_2D &
                                (FEXS_photon_electron_Compton_integrand, &
                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, l_integral)
                            integral(np,n,gpr,g,lp+1) = l_integral
                        end do
                    end do
                end do
            end do
        end do
    end if

    Sigmapec = 0
    do g = 1, Ge
        do gpr = 1, Gp
            do n = 1, 2
                do np = 1, 2
                    do l = 0, NL
                        do lp = 0, l
                            Sigmapec(np,n,gpr,g,l+1) = Sigmapec(np,n,gpr,g,l+1) + &
                                ((-1)**(n+np))*(e_radius**2)*(2*l+1)*(e_mass_E**3)*coeff(l,lp)*&
                                integral(np,n,gpr,g,lp+1)/(4*dEp(gpr)*dEe(g))
                            !Sigmapec(np,n,gpr,g,l+1) = Sigmapec(np,n,gpr,g,l+1) + &
                            !    (e_radius**2)*(2*l+1)*e_mass_E*coeff(l,lp)*&
                            !    integral(np,n,gpr,g,lp+1)/(4*dEp(gpr))
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine FEXS_photon_electron_Compton

subroutine FEXS_photon_electron_pair_production(Ep, dEp, Ee, dEe, Eemid, Nint, Sigmapepp)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Sigmapepp

    integer :: ele, mat, g, gpr, i, j, l, n, np
    integer :: NNIST
    integer :: Gpp
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:), allocatable :: ENist
    real, dimension(:), allocatable :: ppxs
    real, dimension(:,:), allocatable :: angular
    real, dimension(:,:,:,:), allocatable :: Fints
    real, dimension(:), allocatable :: norms
    real :: l_eps
    real :: l_delta
    real :: l_gamma
    real :: l_beta
    real :: l_E
    real, dimension(:,:,:,:), allocatable :: yvals
    real :: l_integral
    real, dimension(:,:,:,:,:), allocatable :: integrals

        !real, dimension(:), allocatable :: BETAS
        !
        !allocate(angular(NL+1,Ge))
        !allocate(BETAS(Ge))
        !
        !do i = 1, Ge
        !    BETAS(i) = beta(0.5*(Ee(i)+Ee(i+1)))
        !    !print *, BETAS(i)
        !end do
        !
        !do i = 1, Ge
        !    angular(1,i) = 1
        !    angular(2,i) = &
        !        (2*BETAS(i) + (1-BETAS(i)**2)*log((1-BETAS(i))/(1+BETAS(i))))/(2*BETAS(i)**2)
        !    do l = 2, NL
        !        angular(l+1,i) = &
        !        ((2*l-1)*angular(l,i) - l*BETAS(i)*angular(l-1,i))/((l-1)*BETAS(i))
        !        if (l .ge. 6 .and. BETAS(i) .lt. 0.5) angular(l+1,i) = 0
        !    end do
        !end do
        !
        !!where (angular .le. 1.0E-8) angular = 0
        !
        !print *, "BREAK"
        !do l = 5, 5
        !    print *, "STARTING PRINTOUT FOR:", l
        !    do i = 1, Ge
        !        print *, angular(l+1,i)
        !    end do
        !end do
        !
        !stop

    Gpp = findloc(Ep, 1.022, dim = 1)
    if (Gpp .eq. 0) then
        print *, "MGXS pair production: Photon energy grouping does not contain 1.022 MeV (2*electron mass energy)."
        stop
    end if
    Gpp = Gpp - 1 ! Last group for which PP is valid

    if (findloc(Ep, 2.044, dim = 1) .eq. 0) then
        print *, "MGXS pair production: Photon energy grouping does not contain 2.044 MeV (4*electron mass energy)."
        stop
    end if

    allocate(Sigmapepp(2,2,Gp,Ge,NL+1,Nmats))
    allocate(angular(NL+1,Nint))
    allocate(integrals(2,2,Gpp,Ge,NL+1))

    Sigmapepp = 0

    call populate_Gauss_Legendre_quadrature(Nint, xq, wq)

    do mat = 1, Nmats
        integrals = 0
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            XSZ = eles_list(mat,ele)
            call read_pair_production(eles_list(mat,ele), Epmax, ENIST, ppxs)
            NNIST = size(ENIST)

            allocate(Fints(2,NNIST,NL+1,Ge))
            allocate(yvals(2,NNIST,NL+1,Ge))
            allocate(norms(NNIST))

            ! Form normalization constant of Bethe-Heitler differential energy distribution
            ! for energies in ENIST
            norms = 0
            do i = 1, NNIST
                do j = 1, Nint
                    l_eps = 0.5*(0.5 - 1/ENIST(i))*xq(j) + 0.5*(0.5 + 1/ENIST(i))
                    norms(i) = norms(i) + (0.5-1/ENIST(i))*wq(j)*Bethe_Heitler_DED(ENIST(i),l_eps)
                end do
            end do

            Fints = 0
            ! Do epsilon integrals for a given group g
            do g = 1, Ge
                l_gamma = Ee(g+1)/e_mass_E
                ! Do epsilon integrals for a given energy ENIST(i)
                do i = 1, NNIST
                    l_delta = max(min(Ee(g)/e_mass_E,ENIST(i)-2),Ee(g+1)/e_mass_E)

                    ! PREPARE ANGULAR
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

                    !where (angular .lt. 1.0E-8) angular = 1.0E-8

                    ! Integrate with quadrature. Use epsilon values given by epsilon(ENIST(i),l_delta)
                    ! as upper limit, and epsilon(ENIST(i),l_gamma) as lower limit. Form quadrature node
                    ! with these limits.
                    do l = 0, NL
                        do j = 1, Nint
                            l_eps = 0.5*((l_delta + 1)/ENIST(i) - (l_gamma + 1)/ENIST(i))*xq(j) &
                                  + 0.5*((l_delta + 1)/ENIST(i) + (l_gamma + 1)/ENIST(i))
                            do n = 1, 2
                                Fints(n,i,l+1,g) = Fints(n,i,l+1,g) + &
                                    0.5*(Ee(g-n+2)/e_mass_E - (l_eps*ENIST(i)-1))*&
                                    ((l_delta + 1)/ENIST(i) - (l_gamma + 1)/ENIST(i))*&
                                    wq(j)*angular(l+1,j)*Bethe_Heitler_DED(ENIST(i),l_eps)/norms(i)
                                !Fints(n,i,l+1,g) = Fints(n,i,l+1,g) + &
                                !0.5*((l_delta + 1)/ENIST(i) - (l_gamma + 1)/ENIST(i))*&
                                !wq(j)*angular(l+1,j)*Bethe_Heitler_DED(ENIST(i),l_eps)/norms(i)
                            end do
                        end do
                    end do
                end do
            end do

            do i = 1, NNIST
                do n = 1, 2
                    yvals(n,i,1:NL+1,1:Ge) = Fints(n,i,1:NL+1,1:Ge)*ppxs(i)
                end do
            end do

            do g = 1, Ge
                do gpr = 1, Gpp
                    do l = 0, NL
                        do n = 1, 2
                            do np = 1, 2
                                MATHA = Ep(gpr-np+2)/e_mass_E
                                call trapezoidal_integration_1D &
                                    (2, yvals(n,:,l+1,g), ENIST, Ep(gpr+1)/e_mass_E, Ep(gpr)/e_mass_E, l_integral)
                                integrals(np,n,gpr,g,l+1) = integrals(np,n,gpr,g,l+1) + &
                                    ((-1)**(n+np))*mat_rho(mat)*mass_frac(mat,ele)*(e_mass_E**3)*(2*l+1)*&
                                    l_integral/(twopi*dEp(gpr)*dEe(g))
                                !call trapezoidal_integration_1D &
                                !(0, yvals(n,:,l+1,g), ENIST, Ep(gpr+1)/e_mass_E, Ep(gpr)/e_mass_E, l_integral)
                                !integrals(np,n,gpr,g,l+1) = integrals(np,n,gpr,g,l+1) + &
                                !mat_rho(mat)*mass_frac(mat,ele)*e_mass_E*(2*l+1)*&
                                !l_integral/(twopi*dEp(gpr))
                            end do
                        end do
                    end do
                end do
            end do

            deallocate(ENIST)
            deallocate(ppxs)
            deallocate(Fints)
            deallocate(yvals)
            deallocate(norms)
        end do

        Sigmapepp(1:2,1:2,1:Gpp,1:Ge,1:NL+1,mat) = integrals
    end do
end subroutine FEXS_photon_electron_pair_production

subroutine FEXS_photon_electron_photoelectric_effect(Ep, dEp, Ee, dEe, Eemid, Sigmapepe)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid

    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Sigmapepe

    integer :: ele, g, gpr, i, ip, j, l, mat, n, np, shell
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

    allocate(Sigmapepe(2,2,Gp,Ge,NL+1,Nmats))
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

                do gpr = 1, Gp
                    do g = 1, Ge
                        l_alpha = min(max(Ep(gpr+1),Ee(g+1)+peWvals(shell)),Ep(gpr))/e_mass_E
                        l_beta = max(min(Ep(gpr),Ee(g)+peWvals(shell)),Ep(gpr+1))/e_mass_E
                        if (l_alpha .eq. l_beta) cycle
                        do np = 1, 2
                            MATHA = Ep(gpr-np+2)/e_mass_E
                            do n = 1, 2
                                MATHB = (Ee(g-n+2) + peWvals(shell))/e_mass_E
                                do l = 0, NL
                                    call trapezoidal_integration_1D &
                                    (1, apexs(:,l+1), E, l_alpha, l_beta, integral)
                                    Sigmapepe(np,n,gpr,g,l+1,mat) = Sigmapepe(np,n,gpr,g,l+1,mat) + &
                                        ((-1)**(np+n))*mat_rho_a(mat,ele)*(2*l+1)*&
                                        (e_mass_E**3)*integral/(fourpi*dEp(gpr)*dEe(g))
                                    !call trapezoidal_integration_1D &
                                    !(0, apexs(:,l+1), E, l_alpha, l_beta, integral)
                                    !Sigmapepe(np,n,gpr,g,l+1,mat) = Sigmapepe(np,n,gpr,g,l+1,mat) + &
                                    !    mat_rho_a(mat,ele)*(2*l+1)*&
                                    !    e_mass_E*integral/(fourpi*dEp(gpr))
                                end do
                            end do
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
end subroutine FEXS_photon_electron_photoelectric_effect

subroutine FEXS_photon_electron_Auger(Ep, dEp, Ee, dEe, SigmapeA)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: SigmapeA

    integer :: cascade, ele, g, gpr, i, j, mat, n, np, shell
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
    real, dimension(:,:,:), allocatable :: Lambdas
    integer :: globalcascade
    real :: l_alpha
    real :: l_beta
    real :: integral
    real, dimension(:,:), allocatable :: integrals

    allocate(SigmapeA(2,2,Gp,Ge,Nmats))
    allocate(integrals(2,Gp))

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
            allocate(Lambdas(2,Ge,Ncascades))
            do cascade = 1, Ncascades
                do g = 1, Ge
                    if (Ee(g+1) .le. Eij(cascade) .and. Eij(cascade) .le. Ee(g)) then
                        deltas(g,cascade) = 1
                    else
                        deltas(g,cascade) = 0
                    end if
                    !deltas(g,cascade) = merge(1,0,Ee(g+1) .le. Eij(cascade) .and. Eij(cascade) .le. Ee(g))
                    do n = 1, 2
                        Lambdas(n,g,cascade) = ((-1)**n)*(Ee(g-n+2)-Eij(cascade))/dEe(g)
                    end do
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
                    do np = 1, 2
                        MATHA = Ep(gpr-np+2)/e_mass_E
                        call log_log_integration_1D &
                            (2, pexs(:,shell), Espec, l_alpha, l_beta, integrals(np,gpr))
                    end do
                end do

                do cascade = 1, cscpershell(shell)
                    globalcascade = globalcascade + 1
                    do g = 1, Ge
                        do gpr = 1, Gp
                            do n = 1, 2
                                do np = 1, 2
                                    SigmapeA(np,n,gpr,g,mat) = SigmapeA(np,n,gpr,g,mat) + ((-1)**np)*&
                                    mat_rho_a(mat,ele)*deltas(g,globalcascade)*&
                                    Lambdas(n,g,globalcascade)*eta(globalcascade)*&
                                    (e_mass_E**2)*integrals(np,gpr)/(fourpi*dEp(gpr))
                                end do
                            end do
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
end subroutine FEXS_photon_electron_Auger

!--------------------- Electron-Electron ----------------------------------------------------------
subroutine FEXS_electron_electron(Ee, dEe, Eemid, adaptive, tol, Nint, Sigmaee, p_Sigmaet, pi_Sigmaet)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:), intent(in) :: Eemid
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Sigmaee
    real, dimension(:,:,:,:), allocatable, intent(inout) :: p_Sigmaet
    real, dimension(:,:), allocatable, intent(inout) :: pi_Sigmaet

    integer :: g, gpr, l, mat
    integer :: l_NL
    real, dimension(:,:,:,:,:,:), allocatable :: Sigmaeeinel
    real, dimension(:,:,:,:,:), allocatable :: RCSDA
    real, dimension(:,:,:,:,:), allocatable :: Sigmaeeel
    real, dimension(:,:,:,:,:,:), allocatable :: SigmaeeA
    real, dimension(:,:,:,:,:,:), allocatable :: B_Sigma

    XSlogic = .true.

    l_NL = NL
    if (exact_RCSDA_angular .and. electron_angular .eq. "SN") l_NL = NL + 1

    allocate(Sigmaee(2,2,Ge,Ge,l_NL+1,Nmats))
    allocate(p_Sigmaet(2,2,Ge,Nmats))
    allocate(pi_Sigmaet(Ge+1,Nmats))

    Sigmaee = 0
    p_Sigmaet = 0
    pi_Sigmaet = 0

    call FEXS_electron_electron_inelastic(Ee, dEe, adaptive, tol, Nint, Sigmaeeinel)
    if (ccut .gt. 0) then
        call FEXS_electron_electron_RCSDA(Ee, dEe, adaptive, tol, Nint, RCSDA, p_Sigmaet)
    end if
    call FEXS_electron_electron_elastic_ELSEPA(Ee, dEe, Sigmaeeel, p_Sigmaet, pi_Sigmaet)

    if (inelastic_etc) then
        do l = 0, NL
            Sigmaee(:,:,:,:,l+1,:) = Sigmaeeinel(:,:,:,:,l+1,:) - (2*l+1)*Sigmaeeinel(:,:,:,:,NL+2,:)/(2*NL+3)
        end do

        Sigmaee(:,:,:,:,NL+2,:) = fourpi*Sigmaeeinel(:,:,:,:,NL+2,:)/(2*NL+3)

        if (ccut .gt. 0) then
            if (exact_RCSDA_angular) then
                Sigmaee(:,:,:,:,NL+2,:) = Sigmaee(:,:,:,:,NL+2,:) + RCSDA
            else
                do l = 0, NL
                    Sigmaee(:,:,:,:,l+1,:) = Sigmaee(:,:,:,:,l+1,:) + (2*l+1)*RCSDA/fourpi
                end do
            end if
        end if

        do g = 1, Ge
            Sigmaee(:,:,g,g,1:NL+1,:) = Sigmaee(:,:,g,g,1:NL+1,:) + Sigmaeeel(:,:,g,:,:)
        end do
    else
        if (ccut .gt. 0) then
            if (exact_RCSDA_angular) then
                Sigmaee(:,:,:,:,NL+2,:) = Sigmaee(:,:,:,:,NL+2,:) + RCSDA
            else
                do l = 0, NL
                    Sigmaee(:,:,:,:,l+1,:) = Sigmaee(:,:,:,:,l+1,:) + (2*l+1)*RCSDA/fourpi
                end do
            end if
        end if

        Sigmaee(:,:,:,:,1:NL+1,:) = Sigmaee(:,:,:,:,1:NL+1,:) + Sigmaeeinel

        do g = 1, Ge
            Sigmaee(:,:,g,g,1:NL+1,:) = Sigmaee(:,:,g,g,1:NL+1,:) + Sigmaeeel(:,:,g,1:NL+1,:)
        end do
    end if

    allocate(B_Sigma(2,2,Ge,Ge,l_NL+1,Nmats))
    do g = 1, Ge
        B_Sigma(1,1,:,g,:,:) = Sigmaee(1,1,:,g,:,:)*4/dEe(g) - Sigmaee(1,2,:,g,:,:)*2/dEe(g)
        B_Sigma(2,1,:,g,:,:) = Sigmaee(2,1,:,g,:,:)*4/dEe(g) - Sigmaee(2,2,:,g,:,:)*2/dEe(g)
        B_Sigma(1,2,:,g,:,:) = -Sigmaee(1,1,:,g,:,:)*2/dEe(g) + Sigmaee(1,2,:,g,:,:)*4/dEe(g)
        B_Sigma(2,2,:,g,:,:) = -Sigmaee(2,1,:,g,:,:)*2/dEe(g) + Sigmaee(2,2,:,g,:,:)*4/dEe(g)
    end do
    Sigmaee = B_Sigma
end subroutine FEXS_electron_electron

subroutine FEXS_electron_electron_inelastic(Ee, dEe, adaptive, tol, Nint, Sigmaeeinel)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Sigmaeeinel

    integer :: ele, g, gpr, l, lp, mat, n, np, s
    integer :: l_NL
    integer :: Nshells
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real :: integral1
    real :: integral2
    real, dimension(:,:,:,:,:), allocatable :: integrals
    real, dimension(:,:), allocatable :: coeff
    real, dimension(:,:,:,:,:), allocatable :: time

    XScut = ccut/e_mass_E

    l_NL = NL
    if (inelastic_etc) l_NL = NL + 1

    allocate(Sigmaeeinel(2,2,Ge,Ge,l_NL+1,Nmats))
    allocate(integrals(2,2,Ge,Ge,l_NL+1))
    allocate(coeff(0:l_NL,0:l_NL))

    !allocate(time(2,2,Ge,Ge,NL+1))

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(Nint, xq1, wq1)
    end if

    coeff = Legendre_polynomial_coeffs(l_NL)

    Sigmaeeinel = 0
    if (electron_inelastic .eq. "Moller") then
        XSw = 0
        do gpr = 1, Ge
            XSa = Ee(gpr+1)/e_mass_E
            XSb = Ee(gpr)/e_mass_E
            do g = gpr, Ge
                XSc = Ee(g+1)/e_mass_E
                XSd = Ee(g)/e_mass_E
                do np = 1, 2
                    XSEgp = Ee(gpr-np+2)/e_mass_E
                    do n = 1, 2
                        XSEg = Ee(g-n+2)/e_mass_E
                        do lp = 0, l_NL
                            XSlp = lp
                            if (adaptive) then
                                call Gauss_Legendre_adaptive_quadrature_2D &
                                    (FEXS_electron_electron_Moller_primary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral1)
                                call Gauss_Legendre_adaptive_quadrature_2D &
                                    (FEXS_electron_electron_Moller_secondary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral2)
                            else
                                call Gauss_Legendre_quadrature_2D &
                                    (FEXS_electron_electron_Moller_primary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral1)
                                call Gauss_Legendre_quadrature_2D &
                                    (FEXS_electron_electron_Moller_secondary_integrand, &
                                    -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral2)
                            end if

                            do l = lp, l_NL
                                do mat = 1, Nmats
                                    Sigmaeeinel(np,n,gpr,g,l+1,1) = Sigmaeeinel(np,n,gpr,g,l+1,1) + &
                                        ((-1)**(n+np))*mat_rho_e(mat)*(e_radius**2)*coeff(l,lp)*&
                                        (2*l+1)*(e_mass_E**3)*(integral1+integral2)/(2*dEe(gpr)*dEe(g))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    else if (electron_inelastic .eq. "RBED") then
        do mat = 1, Nmats
            do ele = 1, Nmaxeles
                if (eles_list(mat,ele) .eq. 0) cycle
                call read_atomic_parameters(eles_list(mat,ele), degen, W, U)
                Nshells = size(degen)

                integrals = 0
                do s = 1, Nshells
                    XSw = W(s)/e_mass_E
                    XSbeta_U2 = beta(U(s))**2
                    XSbeta_W2 = beta(W(s))**2
                    do gpr = 1, Ge
                        XSa = Ee(gpr+1)/e_mass_E
                        XSb = Ee(gpr)/e_mass_E
                        do g = gpr, Ge
                            XSc = Ee(g+1)/e_mass_E
                            XSd = Ee(g)/e_mass_E
                            do np = 1, 2
                                XSEgp = Ee(gpr-np+2)/e_mass_E
                                do n = 1, 2
                                    XSEg = Ee(g-n+2)/e_mass_E
                                    do lp = 0, l_NL
                                        XSlp = lp
                                        !time(np,n,gpr,g,lp+1) = omp_get_wtime()
                                        if (adaptive) then
                                            call Gauss_Legendre_adaptive_quadrature_2D &
                                                (FEXS_electron_electron_RBED_primary_integrand, &
                                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral1)
                                            call Gauss_Legendre_adaptive_quadrature_2D &
                                                (FEXS_electron_electron_RBED_secondary_integrand, &
                                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral2)
                                        else
                                            call Gauss_Legendre_quadrature_2D &
                                                (FEXS_electron_electron_RBED_primary_integrand, &
                                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral1)
                                            call Gauss_Legendre_quadrature_2D &
                                                (FEXS_electron_electron_RBED_secondary_integrand, &
                                                -1.0, 1.0, -1.0, 1.0, xq1, wq1, integral2)
                                        end if
                                        integrals(np,n,gpr,g,lp+1) = integrals(np,n,gpr,g,lp+1) + &
                                            degen(s)*(integral1 + integral2)
                                        !time(np,n,gpr,g,lp+1) = omp_get_wtime() - time(np,n,gpr,g,lp+1)
                                    end do
                                    !if (gpr .ge. 26) print *, XSa, XSb, XSc, XSd, XSEgp, XSEg
                                end do
                            end do
                        end do
                    end do
                end do

                do gpr = 1, Ge
                    do g = gpr, Ge
                        do np = 1, 2
                            do n = 1, 2
                                do l = 0, l_NL
                                    do lp = 0, l
                                        Sigmaeeinel(np,n,gpr,g,l+1,mat) = Sigmaeeinel(np,n,gpr,g,l+1,mat) + &
                                            ((-1)**(np+n))*mat_rho_a(mat,ele)*(e_radius**2)*coeff(l,lp)*&
                                            (2*l+1)*(e_mass_E**3)*integrals(np,n,gpr,g,lp+1)/(2*dEe(gpr)*dEe(g))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do

                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do
        end do
    end if
end subroutine FEXS_electron_electron_inelastic

subroutine FEXS_electron_electron_RCSDA(Ee, dEe, adaptive, tol, Nint, RCSDA, p_Sigmaet)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: RCSDA
    real, dimension(:,:,:,:), intent(inout) :: p_Sigmaet

    integer :: ele, g, gpr, mat, n, np, s
    real :: Ibar
    integer :: Nshells
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real :: integral
    real :: x
    real, dimension(:,:,:), allocatable :: Lints
    real, dimension(:), allocatable :: LDeltas

    XScut = ccut/e_mass_E

    allocate(RCSDA(2,2,Ge,Ge,Nmats))
    allocate(Lints(2,2,Ge))
    allocate(LDeltas(Ge+1))

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(Nint, xq1, wq1)
    end if

    ! NOTE: Boundary condition is assumed that psi = 0. Should I not use boundary condition or what?
    ! The effect this has is on the g = 1 term. See MMS for normal (not assuming any bdy condition, so no sweep)
    ! See others for psi = 0 condition (first line of g = 1 term is assumed 0).

    if (electron_inelastic .eq. "MMS") then
        deallocate(RCSDA, Lints, LDeltas)
        allocate(RCSDA(2,2,size(Ee)-1,size(Ee)-1,1))
        allocate(Lints(2,2,size(Ee)-1))
        allocate(LDeltas(size(Ee)))

        do g = 1, size(Ee)-1
            XSa = Ee(g+1)
            XSb = Ee(g)
            do np = 1, 2
                XSEgp = Ee(g-np+2)
                call Gauss_Legendre_adaptive_quadrature_1D &
                    (RCSDA_MMS_integrand, -1.0, 1.0, xq1, wq1, &
                    xq2, wq2, tol, integral)
                do n = 1, 2
                    Lints(np,n,g) = (-1)**(n+1+np)*integral/(dEe(g)**2)
                end do
            end do
        end do

        !do np = 1, 2
        !    do n = 1, 2
        !        print *, [np,n], Lints(np,n,1)*2
        !    end do
        !end do
        !stop

        do g = 1, size(Ee)
            LDeltas(g) = RCSDA_MMS_LET(Ee(g))
        end do

        do gpr = 1, size(Ee)-1
            do g = gpr, size(Ee)-1
                do np = 1, 2
                    do n = 1, 2
                        if (g .eq. 1) then
                            RCSDA(np,n,gpr,1,1) = & ! NOTE: MMS like this doesn't assume zero bdy cond. Should do that though... Or use a bdycond.
                                Kron_delta(n,2)*Kron_delta(np,1)*Kron_delta(gpr,1)*LDeltas(gpr+np-1) &
                                - Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,1)*LDeltas(gpr+np-1) &
                                - Kron_delta(gpr,1)*Lints(np,n,1)
                        else
                            RCSDA(np,n,gpr,g,1) = &
                                Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,g-1)*LDeltas(gpr+np-1) &
                                - Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,g)*LDeltas(gpr+np-1) &
                                - Kron_delta(gpr,g)*Lints(np,n,g)
                        end if
                    end do
                end do
            end do
        end do
        print *, "GOTTA LOOK AT SIGMAT FOR RCSDA FEXS"
        stop
        return
    end if

    RCSDA = 0
    do mat = 1, Nmats
        Lints = 0
        LDeltas = 0

        if (electron_inelastic .eq. "Moller") then
            XSw = 0
            call mean_excitation_energy(eles_list(mat,:), Natoms(mat,:), Ibar)
            call set_density_effect_parameters &
                (mat_rho(mat), total_GAM(mat), eles_list(mat,:), Natoms(mat,:), Ibar, .false.)

            XSI = Ibar*1.0E-6/e_mass_E

            do g = 1, Ge
                XSa = Ee(g+1)/e_mass_E
                XSb = Ee(g)/e_mass_E
                do np = 1, 2
                    XSEgp = Ee(g-np+2)/e_mass_E

                    if (adaptive) then
                        call Gauss_Legendre_adaptive_quadrature_1D &
                            (FEXS_electron_electron_Moller_RCSDA_integrand, -1.0, 1.0, xq1, wq1, &
                            xq2, wq2, tol, integral)
                    else
                        call Gauss_Legendre_quadrature_1D &
                            (FEXS_electron_electron_Moller_RCSDA_integrand, -1.0, 1.0, xq1, wq1, &
                            integral)
                    end if

                    do n = 1, 2
                        Lints(np,n,g) = &
                            (-1)**(n+1+np)*mat_rho_e(mat)*twopi*(e_radius**2)*&
                            (e_mass_E**3)*integral/(dEe(g)**2)
                    end do
                end do
            end do

            do g = 1, Ge+1
                x = Ee(g)/e_mass_E
                LDeltas(g) = mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*unitless_Moller_RLET(x)
            end do

            do g = 1, Ge
                do np = 1, 2
                    do n = 1, 2
                        p_Sigmaet(np,n,g,mat) = p_Sigmaet(np,n,g,mat) + &
                            Kron_delta(n,2)*Kron_delta(np,2)*LDeltas(g+np-1) + &
                            Lints(np,n,g)
                        if (g .eq. 1) cycle
                        RCSDA(np,n,g-1,g,mat) = Kron_delta(n,1)*Kron_delta(np,2)*LDeltas(g-1+np-1)
                    end do
                end do
            end do

            !do gpr = 1, Ge
            !    do g = gpr, Ge
            !        do np = 1, 2
            !            do n = 1, 2
            !                if (g .eq. 1) then
            !                    RCSDA(np,n,gpr,1,mat) = &
            !                        - Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,1)*LDeltas(gpr+np-1) &
            !                        - Kron_delta(gpr,1)*Lints(np,n,1)
            !                else
            !                    RCSDA(np,n,gpr,g,mat) = &
            !                        Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,g-1)*LDeltas(gpr+np-1) &
            !                        - Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,g)*LDeltas(gpr+np-1) &
            !                        - Kron_delta(gpr,g)*Lints(np,n,g)
            !                end if
            !            end do
            !        end do
            !    end do
            !end do
        else if (electron_inelastic .eq. "RBED") then
            do ele = 1, Nmaxeles
                if (eles_list(mat,ele) .eq. 0) cycle
                call read_atomic_parameters(eles_list(mat,ele), degen, W, U)
                Nshells = size(degen)

                do s = 1, Nshells
                    XSw = W(s)/e_mass_E
                    XSbeta_U2 = beta(U(s))**2
                    XSbeta_W2 = beta(W(s))**2
                    do g = 1, Ge
                        XSa = Ee(g+1)/e_mass_E
                        XSb = Ee(g)/e_mass_E
                        do np = 1, 2
                            XSEgp = Ee(g-np+2)/e_mass_E

                            if (adaptive) then
                                call Gauss_Legendre_adaptive_quadrature_1D &
                                    (FEXS_electron_electron_RBED_RCSDA_integrand, -1.0, 1.0, &
                                    xq1, wq1, xq2, wq2, tol, integral)
                            else
                                call Gauss_Legendre_quadrature_1D &
                                    (FEXS_electron_electron_RBED_RCSDA_integrand, -1.0, 1.0, &
                                    xq1, wq1, integral)
                            end if

                            do n = 1, 2
                                Lints(np,n,g) = Lints(np,n,g) + &
                                    (-1)**(n+1+np)*mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius**2)*&
                                    (e_mass_E**3)*integral/(dEe(g)**2)
                            end do
                        end do
                    end do

                    do g = 1, Ge+1
                        x = Ee(g)/e_mass_E
                        LDeltas(g) = LDeltas(g) + &
                            mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius**2)*e_mass_E*&
                            unitless_RBED_RLET(x)
                    end do
                end do

                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do

            do g = 1, Ge
                do np = 1, 2
                    do n = 1, 2
                        p_Sigmaet(np,n,g,mat) = p_Sigmaet(np,n,g,mat) + &
                            Kron_delta(n,2)*Kron_delta(np,2)*LDeltas(g+np-1) + &
                            Lints(np,n,g)
                        if (g .eq. 1) cycle
                        RCSDA(np,n,g-1,g,mat) = Kron_delta(n,1)*Kron_delta(np,2)*LDeltas(g-1+np-1)
                    end do
                end do
            end do

            !do gpr = 1, Ge
            !    do g = gpr, Ge
            !        do np = 1, 2
            !            do n = 1, 2
            !                if (g .eq. 1) then
            !                    RCSDA(np,n,gpr,1,mat) = &
            !                        - Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,1)*LDeltas(gpr+np-1) &
            !                        - Kron_delta(gpr,1)*Lints(np,n,1)
            !                else
            !                    RCSDA(np,n,gpr,g,mat) = &
            !                        Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,g-1)*LDeltas(gpr+np-1) &
            !                        - Kron_delta(n,2)*Kron_delta(np,2)*Kron_delta(gpr,g)*LDeltas(gpr+np-1) &
            !                        - Kron_delta(gpr,g)*Lints(np,n,g)
            !                end if
            !            end do
            !        end do
            !    end do
            !end do
        end if
    end do
end subroutine FEXS_electron_electron_RCSDA

subroutine FEXS_electron_electron_elastic_ELSEPA(Ee, dEe, Sigmaeeel, p_Sigmaet, pi_Sigmaet)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: Sigmaeeel
    real, dimension(:,:,:,:), intent(inout) :: p_Sigmaet
    real, dimension(:,:), intent(inout) :: pi_Sigmaet

    integer :: ele, g, i, j, l, lp, mat, n, np
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
    real :: Mval
    real :: Bval
    real :: integral
    real, dimension(:,:), allocatable :: xs_ell
    real, dimension(:,:), allocatable :: xs_elli
    real, dimension(:,:,:,:,:), allocatable :: B_Sigmaeeel

    l_NL = NL
    if (elastic_etc) l_NL = NL + 1

    allocate(Sigmaeeel(2,2,Ge,l_NL+1,Nmats))
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
    allocate(xs_ell(NMT,l_NL+1))
    allocate(xs_elli(NEs,l_NL+1))
    if (elastic_etc) allocate(B_Sigmaeeel(2,2,Ge,NL+1,Nmats))

    coeff = Legendre_polynomial_coeffs(l_NL)

    Sigmaeeel = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            write(Zstr, *) eles_list(mat,ele)

            do i = 1, NEs
                power = floor(log10(E(i)))
                second = nint(mod(E(i)*10.0**(-power),1.0)*1000)
                if (second .eq. 1000) then ! SHOULDN'T NEED TO DO THIS, FIX!!
                    second = 0
                end if
                first = nint(E(i)*10.0**(-power) - real(second)/1000)

                write(powerstr, *) power + 6
                write(firststr, *) first

                if (second .eq. 0) then ! CHANGE THESE
                    fname = trim(adjustl(Zstr)) //"/dcs_"&
                    // trim(adjustl(firststr)) //"p000e0" &
                    // trim(adjustl(powerstr)) //".dat"
                else
                    write(secondstr, *) second

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
            end do

            do l = 0, l_NL
                if (any(xs_elli(1:NEs,l+1) .le. 0)) then
                    do g = 1, Ge
                        do np = 1, 2
                            MATHA = Ee(g-np+2)/e_mass_E
                            do n = 1, 2
                                MATHB = Ee(g-n+2)/e_mass_E
                                call trapezoidal_integration_1D &
                                    (1, xs_elli(1:NEs,l+1), B_E, Ee(g+1)/e_mass_E, Ee(g)/e_mass_E, integral)
                                Sigmaeeel(np,n,g,l+1,mat) = Sigmaeeel(np,n,g,l+1,mat) + &
                                    ((-1)**(n+np))*mat_rho_a(mat,ele)*(e_mass_E**3)*integral/(dEe(g)**2)
                            end do
                        end do
                    end do
                else
                    do g = 1, Ge
                        do np = 1, 2
                            MATHA = Ee(g-np+2)/e_mass_E
                            do n = 1, 2
                                MATHB = Ee(g-n+2)/e_mass_E
                                call log_log_integration_1D &
                                    (1, xs_elli(1:NEs,l+1), B_E, Ee(g+1)/e_mass_E, Ee(g)/e_mass_E, integral)
                                Sigmaeeel(np,n,g,l+1,mat) = Sigmaeeel(np,n,g,l+1,mat) + &
                                    ((-1)**(n+np))*mat_rho_a(mat,ele)*(e_mass_E**3)*integral/(dEe(g)**2)
                            end do
                        end do

                        !call log_log_integration_1D &
                        !    (0, xs_elli(1:NEs,l+1), B_E, Ee(g+1)/e_mass_E, Ee(g)/e_mass_E, integral)
                        !Sigmaeeel(1,1,g,l+1,mat) = Sigmaeeel(1,1,g,l+1,mat) + &
                        !    mat_rho_a(mat,ele)*e_mass_E*integral/dEe(g)
                    end do
                end if
            end do

            if (elastic_etc) then
                do l = 0, NL
                    xs_elli(1:NEs,l+1) = xs_elli(1:NEs,l+1) - (2*l+1)*xs_elli(1:NEs,NL+2)/(2*NL+3)
                end do
            end if

            do g = 1, Ge+1
                pi_Sigmaet(g,mat) = pi_Sigmaet(g,mat) + &
                    fourpi*mat_rho_a(mat,ele)*log_interp1d(B_E, xs_elli(1:NEs,1), Ee(g)/e_mass_E)
            end do
        end do

        if (elastic_etc) then
            do l = 0, NL
                B_Sigmaeeel(:,:,:,l+1,mat) = &
                    Sigmaeeel(:,:,:,l+1,mat) - (2*l+1)*Sigmaeeel(:,:,:,NL+2,mat)/(2*NL+3)
            end do
        end if
    end do

    call move_alloc(B_Sigmaeeel, Sigmaeeel)

    p_Sigmaet = p_Sigmaet + fourpi*Sigmaeeel(1:2,1:2,1:Ge,1,1:Nmats)
end subroutine FEXS_electron_electron_elastic_ELSEPA

subroutine FEXS_electron_electron_Auger(Ee, dEe, adaptive, tol, Nint, SigmaeeA) ! WIP (~6/11/23)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:), allocatable, intent(inout) :: SigmaeeA

    integer :: cascade, ele, g, gpr, i, j, mat, shell
    integer :: Nshells
    integer :: Ncascades
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: Wvals ! MeV
    real, dimension(:), allocatable :: kinetic ! MeV
    integer, dimension(:), allocatable :: cscpershell
    real, dimension(:), allocatable :: Eij ! MeV
    real, dimension(:), allocatable :: eta
    integer, dimension(:,:), allocatable :: deltas
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    integer :: globalcascade
    real :: l_alpha
    real :: l_beta
    real :: integral
    real, dimension(:), allocatable :: integrals

    allocate(SigmaeeA(Ge,Ge,Nmats))
    allocate(integrals(Ge))

    SigmaeeA = 0
    integrals = 0

    call populate_Gauss_Legendre_quadrature(Nint, xq, wq)
    call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
    call populate_Gauss_Legendre_quadrature(8, xq2, wq2)

    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) == 0) cycle
            if (eles_list(mat,ele) .lt. 6) cycle

            ! READ PARAMETERS
            call read_atomic_parameters(eles_list(mat,ele), degen, Wvals, kinetic)
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
            do shell = 1, Nshells
                XSbeta_U2 = beta(kinetic(shell))**2
                XSbeta_W2 = beta(Wvals(shell))**2
                XSw = Wvals(shell)/e_mass_E
                if (adaptive) then
                    integrals = 0
                    do gpr = 1, Ge
                        if (electron_exc_shells) then
                            if (Ee(gpr+1) + 1.0E-6 .lt. Wvals(shell)) cycle
                            l_alpha = Ee(gpr+1)/e_mass_E
                            l_beta = Ee(gpr)/e_mass_E
                        else
                            l_alpha = min(max(Ee(gpr+1),Wvals(shell)),Ee(gpr))/e_mass_E
                            l_beta = Ee(gpr)/e_mass_E
                        end if
                        call Gauss_Legendre_adaptive_quadrature_1D &
                            (unitless_RBEB, l_alpha, l_beta, xq1, wq1, xq2, wq2, tol, integral)
                        integrals(gpr) = twopi*Bohr**2*fsa**4*degen(shell)*e_mass_E*&
                            integral
                    end do
                else
                    print *, "Auger e-> e, non-adaptive WIP, program ending"
                    stop
                end if

                do cascade = 1, cscpershell(shell)
                    globalcascade = globalcascade + 1
                    do g = 1, Ge
                        do gpr = 1, Ge
                            SigmaeeA(gpr,g,mat) = SigmaeeA(gpr,g,mat) + &
                                mat_rho_a(mat,ele)*deltas(g,globalcascade)*eta(globalcascade)*&
                                integrals(gpr)/(fourpi*dEe(gpr))
                        end do
                    end do
                end do
            end do


            deallocate(deltas)
    1       deallocate(degen)
            deallocate(Wvals)
            deallocate(kinetic)
            deallocate(cscpershell)
            deallocate(Eij)
            deallocate(eta)
        end do
    end do
end subroutine FEXS_electron_electron_Auger

!--------------------- Electron-Photon ------------------------------------------------------------

! TOTAL ATTENUATION
!--------------------- Photon ---------------------------------------------------------------------
subroutine FEXS_photon_attn(Ep, dEp, adaptive, tol, Nint, Sigmapt, i_Sigmapt)
    implicit none
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: dEp
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigmapt
    real, dimension(:,:), allocatable, intent(inout) :: i_Sigmapt

    integer :: ele, g, gpr, i, mat, n, np
    integer :: Gpp
    integer :: NNIST
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:,:), allocatable :: Sigmapct
    real, dimension(:,:,:,:), allocatable :: Sigmaptpp
    real, dimension(:,:,:,:), allocatable :: Sigmaptpe
    real :: integral
    real, dimension(:), allocatable :: ENist
    real, dimension(:), allocatable :: ppxs
    real, dimension(:), allocatable :: peWvals
    integer :: NEspec
    real, dimension(:), allocatable :: Espec
    real, dimension(:,:), allocatable :: pexs
    real, dimension(:), allocatable :: pexst
    real, dimension(:), allocatable :: i_pexs
    real, dimension(:,:,:,:), allocatable :: B_Sigmat

    allocate(Sigmapt(2,2,Gp,Nmats))
    allocate(Sigmapct(2,2,Gp))
    allocate(Sigmaptpp(2,2,Gp,Nmats))
    allocate(Sigmaptpe(2,2,Gp,Nmats))
    allocate(i_Sigmapt(Gp+1,Nmats))

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
    i_Sigmapt = 0
    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
        do g = 1, Gp
            XSa = Ep(g+1)/e_mass_E
            XSb = Ep(g)/e_mass_E
            do np = 1, 2
                XSEgp = Ep(g-np+2)/e_mass_E
                do n = 1, 2
                    XSEg = Ep(g-n+2)/e_mass_E
                    call Gauss_Legendre_adaptive_quadrature_1D &
                        (FEXS_photon_total_Compton_integrand, -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                    Sigmapct(np,n,g) = &
                        ((-1)**(n+np))*twopi*(e_radius**2)*(e_mass_E**3)*integral/(dEp(g)**2)
                end do
            end do
        end do
    else
        call populate_Gauss_Legendre_quadrature(Nint, xq1, wq1)
        do g = 1, Gp
            XSa = Ep(g+1)/e_mass_E
            XSb = Ep(g)/e_mass_E
            do np = 1, 2
                XSEgp = Ep(g-np+2)/e_mass_E
                do n = 1, 2
                    XSEg = Ep(g-n+2)/e_mass_E
                    call Gauss_Legendre_quadrature_1D &
                        (FEXS_photon_total_Compton_integrand, -1.0, 1.0, xq1, wq1, integral)
                    Sigmapct(np,n,g) = &
                        ((-1)**(n+np))*twopi*(e_radius**2)*(e_mass_E**3)*integral/(dEp(g)**2)
                end do
            end do
        end do
    end if

    do g = 1, Gp+1
        i_Sigmapt(g,1:Nmats) = i_Sigmapt(g,1:Nmats) + &
            mat_rho_e(1:Nmats)*twopi*(e_radius**2)*unitless_sigmapct(Ep(g)/e_mass_E)
    end do

    Sigmaptpp = 0
    Sigmaptpe = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            if (Epmax .ge. 1.022) then
                call read_pair_production(eles_list(mat,ele), Epmax, ENist, ppxs)
                NNIST = size(ENist)

                do g = 1, Gpp
                    do np = 1, 2
                        MATHA = Ep(g-np+2)/e_mass_E
                        do n = 1, 2
                            MATHB = Ep(g-n+2)/e_mass_E
                            call log_log_integration_1D &
                                (1, ppxs, ENist, Ep(g+1)/e_mass_E, Ep(g)/e_mass_E, integral)
                            Sigmaptpp(np,n,g,mat) = Sigmaptpp(np,n,g,mat) + &
                                ((-1)**(n+np))*mat_rho(mat)*mass_frac(mat,ele)*(e_mass_E**3)*&
                                integral/(dEp(g)**2)
                        end do
                    end do
                end do

                do g = 1, Gpp+1
                    i_Sigmapt(g,mat) = i_Sigmapt(g,mat) + &
                        mat_rho(mat)*mass_frac(mat,ele)*interp1d(ENist,ppxs,Ep(g)/e_mass_E)
                end do
            end if

            call read_photoelectric_effect(eles_list(mat,ele), Epmin, Epmax, peWvals, Espec, pexs)
            NEspec = size(Espec)

            allocate(pexst(NEspec))

            pexst = pexs(:,size(peWvals)+1)

            do g = 1, Gp
                do np = 1, 2
                    MATHA = Ep(g-np+2)/e_mass_E
                    do n = 1, 2
                        MATHB = Ep(g-n+2)/e_mass_E
                        call log_log_integration_1D &
                            (1, pexst, Espec, Ep(g+1)/e_mass_E, Ep(g)/e_mass_E, integral)
                        Sigmaptpe(np,n,g,mat) = Sigmaptpe(np,n,g,mat) + &
                            ((-1)**(n+np))*mat_rho_a(mat,ele)*(e_mass_E**3)*integral/(dEp(g)**2)
                    end do
                end do
            end do

            do g = 1, Gp+1
                i_Sigmapt(g,mat) = i_Sigmapt(g,mat) + &
                    mat_rho_a(mat,ele)*interp1d(Espec,pexst,Ep(g)/e_mass_E)
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
        Sigmapt(1:2,1:2,1:Gp,mat) = &
        mat_rho_e(mat)*Sigmapct + Sigmaptpp(1:2,1:2,1:Gp,mat) + Sigmaptpe(1:2,1:2,1:Gp,mat)
    end do

    allocate(B_Sigmat(2,2,Gp,Nmats))
    do g = 1, Gp
        B_Sigmat(1,1,g,:) = Sigmapt(1,1,g,:)*4/dEp(g) - Sigmapt(1,2,g,:)*2/dEp(g)
        B_Sigmat(2,1,g,:) = Sigmapt(2,1,g,:)*4/dEp(g) - Sigmapt(2,2,g,:)*2/dEp(g)
        B_Sigmat(1,2,g,:) = -Sigmapt(1,1,g,:)*2/dEp(g) + Sigmapt(1,2,g,:)*4/dEp(g)
        B_Sigmat(2,2,g,:) = -Sigmapt(2,1,g,:)*2/dEp(g) + Sigmapt(2,2,g,:)*4/dEp(g)
    end do
    Sigmapt = B_Sigmat
end subroutine FEXS_photon_attn

!--------------------- Electron -------------------------------------------------------------------
subroutine FEXS_electron_attn(Ee, dEe, p_Sigmaet, pi_Sigmaet, adaptive, tol, Nint, Sigmaet, i_Sigmaet)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe
    real, dimension(:,:,:,:), intent(in) :: p_Sigmaet
    real, dimension(:,:), intent(in) :: pi_Sigmaet
    logical, intent(in) :: adaptive
    real, intent(in) :: tol
    integer, intent(in) :: Nint

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigmaet
    real, dimension(:,:), allocatable, intent(inout) :: i_Sigmaet

    integer :: ele, g, mat, np, n, s
    integer :: Nshells
    real, dimension(:), allocatable :: xq1
    real, dimension(:), allocatable :: wq1
    real, dimension(:), allocatable :: xq2
    real, dimension(:), allocatable :: wq2
    real, dimension(:,:,:,:), allocatable :: Sigmaemt
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    real :: integral
    real :: l_alpha
    real :: l_beta
    real :: x
    real, dimension(:,:,:,:), allocatable :: B_Sigmat

    XScut = ccut/e_mass_E

    allocate(Sigmaet(2,2,Ge,Nmats))
    allocate(Sigmaemt(2,2,Ge,Nmats))
    allocate(i_Sigmaet(Ge+1,Nmats))

    Sigmaet = 0
    Sigmaemt = 0
    i_Sigmaet = pi_Sigmaet

    if (adaptive) then
        call populate_Gauss_Legendre_quadrature(5, xq1, wq1)
        call populate_Gauss_Legendre_quadrature(8, xq2, wq2)
    else
        call populate_Gauss_Legendre_quadrature(Nint, xq1, wq1)
    end if

    if (electron_inelastic .eq. "Moller") then
        XSw = 0

        do g = 1, Ge
            XSa = Ee(g+1)/e_mass_E
            XSb = Ee(g)/e_mass_E
            do np = 1, 2
                XSEgp = Ee(g-np+2)/e_mass_E
                do n = 1, 2
                    XSEg = Ee(g-n+2)/e_mass_E

                    if (adaptive) then
                        call Gauss_Legendre_adaptive_quadrature_1D &
                        (FEXS_electron_Moller_R_total_integrand, &
                        -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                    else
                        call Gauss_Legendre_quadrature_1D &
                        (FEXS_electron_Moller_R_total_integrand, &
                        -1.0, 1.0, xq1, wq1, integral)
                    end if

                    do mat = 1, Nmats
                        Sigmaemt(np,n,g,mat) = mat_rho_e(mat)*&
                            ((-1)**(np+n))*twopi*(e_radius**2)*(e_mass_E**3)*integral/(dEe(g)**2)
                    end do
                end do
            end do
        end do

        do mat = 1, Nmats
            do g = 1, Ge+1
                x = Ee(g)/e_mass_E
                i_Sigmaet(g,mat) = i_Sigmaet(g,mat) + &
                    mat_rho_e(mat)*twopi*(e_radius**2)*unitless_Moller_R_Sigmat(x)
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
                    XSbeta_U2 = beta(U(s))**2
                    XSbeta_W2 = beta(W(s))**2
                    do g = 1, Ge
                        XSa = Ee(g+1)/e_mass_E
                        XSb = Ee(g)/e_mass_E
                        do np = 1, 2
                            XSEgp = Ee(g-np+2)/e_mass_E
                            do n = 1, 2
                                XSEg = Ee(g-n+2)/e_mass_E

                                if (adaptive) then
                                    call Gauss_Legendre_adaptive_quadrature_1D &
                                    (FEXS_electron_RBED_total_integrand, &
                                    -1.0, 1.0, xq1, wq1, xq2, wq2, tol, integral)
                                else
                                    call Gauss_Legendre_quadrature_1D &
                                    (FEXS_electron_RBED_total_integrand, &
                                    -1.0, 1.0, xq1, wq1, integral)
                                end if

                                Sigmaemt(np,n,g,mat) = Sigmaemt(np,n,g,mat) + &
                                ((-1)**(np+n))*mat_rho_a(mat,ele)*twopi*degen(s)*(e_radius**2)*&
                                (e_mass_E**3)*integral/(dEe(g)**2)
                            end do
                        end do
                    end do

                    do g = 1, Ge+1
                        x = Ee(g)/e_mass_E
                        if (x .lt. W(s)/e_mass_E) cycle
                        i_Sigmaet(g,mat) = i_Sigmaet(g,mat) + &
                            mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius**2)*&
                            unitless_RBED_Sigmat(x)
                    end do
                end do

                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do
        end do
    end if

    !do g = 1, Ge+1
    !    print *, Ee(g)
    !end do
    !print *, "BREAK"
    !do g = 1, Ge+1
    !    print *, 1.0/i_Sigmaet(g,1)
    !end do
    !stop

    do mat = 1, Nmats
        Sigmaet(1:2,1:2,1:Ge,mat) = Sigmaemt(1:2,1:2,1:Ge,mat) + p_Sigmaet(1:2,1:2,1:Ge,mat)
    end do

    allocate(B_Sigmat(2,2,Ge,Nmats))
    do g = 1, Ge
        B_Sigmat(1,1,g,:) = Sigmaet(1,1,g,:)*4/dEe(g) - Sigmaet(1,2,g,:)*2/dEe(g)
        B_Sigmat(2,1,g,:) = Sigmaet(2,1,g,:)*4/dEe(g) - Sigmaet(2,2,g,:)*2/dEe(g)
        B_Sigmat(1,2,g,:) = -Sigmaet(1,1,g,:)*2/dEe(g) + Sigmaet(1,2,g,:)*4/dEe(g)
        B_Sigmat(2,2,g,:) = -Sigmaet(2,1,g,:)*2/dEe(g) + Sigmaet(2,2,g,:)*4/dEe(g)
    end do
    Sigmaet = B_Sigmat
end subroutine FEXS_electron_attn

! STOPPING POWER
subroutine FEXS_stopping_power(Ee, dEe, Scol)
    implicit none
    real, dimension(:), intent(in) :: Ee
    real, dimension(:), intent(in) :: dEe

    real, dimension(:,:,:), allocatable, intent(inout) :: Scol

    integer :: ele, g, mat, n, s
    real :: Ibar
    real :: x
    integer :: Nshells
    real, dimension(:), allocatable :: degen
    real, dimension(:), allocatable :: W
    real, dimension(:), allocatable :: U
    real, dimension(:,:), allocatable :: g_Scol

    allocate(Scol(2,Ge,Nmats))
    allocate(g_Scol(Ge+1,Nmats))
    g_Scol = 0

    if (electron_inelastic .eq. "Moller") then
        XSw = 0 ! Not necessary here

        do mat = 1, Nmats
            call mean_excitation_energy(eles_list(mat,:), Natoms(mat,:), Ibar)
            call set_density_effect_parameters &
                (mat_rho(mat), total_GAM(mat), eles_list(mat,:), Natoms(mat,:), Ibar, .false.)

            XSI = Ibar*1.0E-6/e_mass_E
            do g = 1, Ge+1
                x = Ee(g)/e_mass_E
                g_Scol(g,mat) = &
                    mat_rho_e(mat)*twopi*(e_radius**2)*e_mass_E*unitless_Bethe(x)/mat_rho(mat)
            end do
        end do
    else if (electron_inelastic .eq. "RBED") then
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
                    do g = 1, Ge+1
                        x = Ee(g)/e_mass_E

                        g_Scol(g,mat) = g_Scol(g,mat) + &
                            mat_rho_a(mat,ele)*degen(s)*twopi*(e_radius**2)*e_mass_E*&
                            unitless_RBED_LET(x)/mat_rho(mat)
                    end do
                end do

                deallocate(degen)
                deallocate(W)
                deallocate(U)
            end do
        end do
    end if

    do g = 1, Ge
        do n = 1, 2
            Scol(n,g,:) = 0.5*g_Scol(g+n-1,:)*dEe(g)
        end do
    end do

    !print *, "BREAK"
    !do g = 1, Ge+1
    !    print *, Ee(g)
    !end do
    !print *, "BREAK"
    !do g = 1, Ge+1
    !    print *, Scol(g,1)
    !end do
    !stop
end subroutine FEXS_stopping_power

end module FEXS