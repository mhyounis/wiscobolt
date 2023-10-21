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

module t2_ele_iteration
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

contains

! SN - general
subroutine t2_ele_SN_read_field &
    (label, l_G, phiu, phi, psi)
    implicit none
    character(*), intent(in) :: label
    integer, intent(in) :: l_G

    real, dimension(:,:,:), allocatable, intent(inout) :: phiu
    real, dimension(:,:,:), allocatable, intent(inout) :: phi
    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: psi

    character(80) :: file_name

    if (FCS .and. &
        ((label .eq. "electron" .and. &
        transport_mode .eq. "external electron beam") .or. &
        (label .eq. "photon" .and. &
        (transport_mode .eq. "external photon beam" .or. &
        transport_mode .eq. "external photon beam coupled")))) allocate(phi(NE,NKe,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(psi(NE,NKe,Nmu,Nphi,l_G))

    if (FCS .and. &
        ((label .eq. "electron" .and. &
        transport_mode .eq. "external electron beam") .or. &
        (label .eq. "photon" .and. &
        (transport_mode .eq. "external photon beam" .or. &
        transport_mode .eq. "external photon beam coupled")))) then
        file_name = trim(adjustl(output_fname))//"/"//label//"_uncollided_fluence.dat"

        open(1, file = file_name, form = "unformatted", action = "read")
        read(1) phiu
        close(1)
    end if

    file_name = trim(adjustl(output_fname))//"/"//label//"_fluence.dat"

    open(1, file = file_name, form = "unformatted", action = "read")
    read(1) phi
    close(1)

    file_name = trim(adjustl(output_fname))//"/"//label//"_angular_fluence.dat"

    open(1, file = file_name, form = "unformatted", action = "read")
    read(1) psi
    close(1)
end subroutine t2_ele_SN_read_field

subroutine t2_ele_SN_write_field &
    (label, psi, phi)
    implicit none
    character(*), intent(in) :: label
    real, dimension(:,:,:,:,:), intent(in) :: psi
    real, dimension(:,:,:), intent(in) :: phi

    character(80) :: file_name

    file_name = trim(adjustl(output_fname))//"/"//label//"_angular_fluence.dat"

    open(1, file = file_name, form = "unformatted")
    write(1) psi
    close(1)

    file_name = trim(adjustl(output_fname))//"/"//label//"_fluence.dat"

    open(1, file = file_name, form = "unformatted")
    write(1) phi
    close(1)
end subroutine t2_ele_SN_write_field

! SN - MGXS
subroutine ele_SN_MGXS_ETC_uncollided_fluence &
    (rglobal, Cekk, bdyel, wL, khat, Pa, factorialmat, Pmlk, cmjk, esweeplist, esweepbounds, &
    enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, Sigma, &
    Sigmat, phiu)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:), intent(in) :: Pa
    real, dimension(:,:), intent(in) :: factorialmat
    real, dimension(:,:,:), intent(in) :: Pmlk
    real, dimension(:,:,:), intent(in) :: cmjk
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:,:,:), intent(in) :: Sigma
    real, dimension(:,:), intent(in) :: Sigmat

    real, dimension(:,:,:), intent(inout) :: phiu

    integer :: e, g, i, j, k, l
    integer :: l_G
    integer :: l_NK
    real, dimension(:,:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: singularity
    real, dimension(:,:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:,:), allocatable :: Tprimeinv
    real, dimension(:,:,:,:), allocatable :: nobdysrc

    l_G = size(phiu,3)
    l_NK = size(nodesinbeam)

    allocate(tempA(NE,NKe))
    allocate(singularity(NK,Nmu,Nphi))
    allocate(t_source(NE,NKe,Nmu,Nphi))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

    singularity = 0
    if (beam_angular_dist .eq. "spherical") then
        !$OMP parallel do collapse(3) shared(singularity) private(i, j, k, l)
        do j = 1, Nphi
            do i = 1, Nmu
                do k = 1, l_NK
                    do l = 0, NL
                        singularity(nodesinbeam(k),i,j) = singularity(nodesinbeam(k),i,j) + &
                            (2*l+1)*sum(factorialmat(l+1,1:l+1)*Pa(i,l+1,1:l+1)*Pmlk(l+1,1:l+1,k)*cmjk(j,1:l+1,k))&
                            /fourpi
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do collapse(2) shared(singularity) private(i, j, k, l)
        do j = 1, Nphi
            do i = 1, Nmu
                do l = 0, NL
                    singularity(nodesinbeam,i,j) = singularity(nodesinbeam,i,j) + &
                        (2*l+1)*sum(factorialmat(l+1,1:l+1)*Pa(i,l+1,1:l+1)*Pmlk(l+1,1:l+1,1)*cmjk(j,1:l+1,1))&
                        /fourpi
                end do
            end do
        end do
        !$OMP end parallel do
    end if

    do g = 2, l_G
        if (Nmats .eq. 1) then
            do k = 1, NKe
                tempA(1:NE,k) = matmul(phiu(1:NE,k,1:g-1),Sigma(1:g-1,g,1))
            end do
            do j = 1, Nphi
                do i = 1, Nmu
                    do k = 1, NKe
                        do e = 1, NE
                            t_source(e,k,i,j) = tempA(e,k)*singularity(Cekk(e,k),i,j)
                        end do
                    end do
                end do
            end do
        else
            ! EXTREMELY UNOPTIMIZED
            tempA = 0
            do k = 1, NKe
                do e = 1, NE
                    tempA(e,k) = tempA(e,k) + &
                        dot_product(phiu(e,k,1:g-1),Sigma(1:g-1,g,eltomat(e)))
                end do
            end do
            do k = 1, NKe
                do e = 1, NE
                    t_source(e,k,1:Nmu,1:Nphi) = tempA(e,k)*singularity(Cekk(e,k),1:Nmu,1:Nphi)
                end do
            end do
        end if

        call t2_ele_SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat(g,:), Tprimeinv)
        call t2_ele_SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
            varsigmadown, I1invI2f, nobdysrc, Tprimeinv, t_source)

        tempA = 0
        do k = 1, NKe
            do j = 1, Nphi
                tempA(1:NE,k) = tempA(1:NE,k) + twopi*matmul(t_source(1:NE,k,1:Nmu,j),wL)/Nphi
            end do
        end do
        print *, norm2(tempA)

        phiu(1:NE,1:NKe,g) = phiu(1:NE,1:NKe,g) + tempA
    end do
end subroutine ele_SN_MGXS_ETC_uncollided_fluence

subroutine temp_TESTING &
    (rglobal, Cekk, bdyel, wL, khat, Pa, factorialmat, Pmlk, cmjk, esweeplist, esweepbounds, &
    enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, Sigma, &
    bdysrc, fEg, Sigmat, phiu, phiu2)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:), intent(in) :: Pa
    real, dimension(:,:), intent(in) :: factorialmat
    real, dimension(:,:,:), intent(in) :: Pmlk
    real, dimension(:,:,:), intent(in) :: cmjk
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:), intent(in) :: bdysrc
    real, dimension(:), intent(in) :: fEg
    real, dimension(:,:), intent(in) :: Sigmat

    real, dimension(:,:,:), allocatable, intent(inout) :: phiu
    real, dimension(:,:,:), allocatable, intent(inout) :: phiu2

    integer :: e, g, i, j, k, l
    integer :: l_G
    integer :: l_NK
    real, dimension(:,:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: singularity
    real, dimension(:,:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:,:), allocatable :: Tprimeinv
    real, dimension(:,:,:,:), allocatable :: t_bdysrc

    l_G = size(fEg)
    l_NK = size(nodesinbeam)

    allocate(phiu(NE,NKe,l_G))
    allocate(phiu2(NE,NKe,l_G))
    allocate(t_bdysrc(NKe,size(bdysrc,2),Nmu,Nphi))
    allocate(tempA(NE,NKe))
    allocate(singularity(NK,Nmu,Nphi))
    allocate(t_source(NE,NKe,Nmu,Nphi))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

    singularity = 0
    if (beam_angular_dist .eq. "spherical") then
        !$OMP parallel do collapse(3) shared(singularity) private(i, j, k, l)
        do j = 1, Nphi
            do i = 1, Nmu
                do k = 1, l_NK
                    do l = 0, NL
                        singularity(nodesinbeam(k),i,j) = singularity(nodesinbeam(k),i,j) + &
                            (2*l+1)*sum(factorialmat(l+1,1:l+1)*Pa(i,l+1,1:l+1)*Pmlk(l+1,1:l+1,k)*cmjk(j,1:l+1,k))&
                            /fourpi
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do collapse(2) shared(singularity) private(i, j, k, l)
        do j = 1, Nphi
            do i = 1, Nmu
                do l = 0, NL
                    singularity(nodesinbeam,i,j) = singularity(nodesinbeam,i,j) + &
                        (2*l+1)*sum(factorialmat(l+1,1:l+1)*Pa(i,l+1,1:l+1)*Pmlk(l+1,1:l+1,1)*cmjk(j,1:l+1,1))&
                        /fourpi
                end do
            end do
        end do
        !$OMP end parallel do
    end if

    t_source = 0
    t_bdysrc = fEg(1)*bdysrc
    call t2_ele_SN_MGXS_construct_inverted_transport_matrix &
        (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat(g,:), Tprimeinv)
    call t2_ele_SN_MGXS_sweep &
        (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
        varsigmadown, I1invI2f, t_bdysrc, Tprimeinv, t_source)
    tempA = 0
    do k = 1, NKe
        do j = 1, Nphi
            tempA(1:NE,k) = tempA(1:NE,k) + twopi*matmul(t_source(1:NE,k,1:Nmu,j),wL)/Nphi
        end do
    end do
    phiu(1:NE,1:NKe,1) = tempA
    phiu2(1:NE,1:NKe,1) = tempA

    do g = 2, l_G
        if (Nmats .eq. 1) then
            do k = 1, NKe
                tempA(1:NE,k) = matmul(phiu(1:NE,k,1:g-1),Sigma(1:g-1,g,1))
            end do
            do k = 1, NKe
                do e = 1, NE
                    t_source(e,k,1:Nmu,1:Nphi) = tempA(e,k)*singularity(Cekk(e,k),1:Nmu,1:Nphi)
                end do
            end do
        else
            ! EXTREMELY UNOPTIMIZED
            do k = 1, NKe
                do e = 1, NE
                    t_source(e,k,1:Nmu,1:Nphi) = &
                        dot_product(phiu(e,k,1:g-1),Sigma(1:g-1,g,eltomat(e)))*&
                        singularity(Cekk(e,k),1:Nmu,1:Nphi)
                end do
            end do
        end if

        t_bdysrc = fEg(g)*bdysrc

        ! MAKING ETC PHIU
        call t2_ele_SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat(g,:), Tprimeinv)
        call t2_ele_SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
            varsigmadown, I1invI2f, t_bdysrc, Tprimeinv, t_source)

        tempA = 0
        do k = 1, NKe
            do j = 1, Nphi
                tempA(1:NE,k) = tempA(1:NE,k) + twopi*matmul(t_source(1:NE,k,1:Nmu,j),wL)/Nphi
            end do
        end do

        phiu(1:NE,1:NKe,g) = tempA

        ! MAKING PHIU ALONE
        t_source = 0
        call t2_ele_SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
            varsigmadown, I1invI2f, t_bdysrc, Tprimeinv, t_source)

        tempA = 0
        do k = 1, NKe
            do j = 1, Nphi
                tempA(1:NE,k) = tempA(1:NE,k) + twopi*matmul(t_source(1:NE,k,1:Nmu,j),wL)/Nphi
            end do
        end do

        phiu2(1:NE,1:NKe,g) = tempA
    end do
end subroutine temp_TESTING

subroutine t2_ele_SN_MGXS_construct_inverted_transport_matrix &
    (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat, Tprimeinv)
    implicit none
    ! Keep? Work into storing the spatial-angular pieces ("alpha" in my notes)?
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    integer, dimension(:,:,:), intent(in) :: enew
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:), intent(in) :: Sigmat

    real, dimension(:,:,:,:,:), intent(inout) :: Tprimeinv

    integer :: dir, e, f, i, j
    real, dimension(:,:), allocatable :: tempC
    real, dimension(:,:), allocatable :: tempD
    real, dimension(:,:), allocatable :: ident

    allocate(tempC(NKe,NKe))
    allocate(tempD(NKe,NKe))
    allocate(ident(NKe,NKe))

    ident = identity(NKe)

    !$OMP parallel do collapse(3) shared(Tprimeinv) private(j, i, e, f, dir) &
    !$OMP private(tempC, tempD)
    do j = 1, Nphi
        do i = 1, Nmu
            do e = 1, NE
                tempD = 0
                do f = 1, NFe
                    tempD = tempD + &
                    varsigmaup(f,e,i,j)*I1invI2(1:NKe,1:NKe,f,e)
                end do
                do dir = 1, 3
                    tempD = tempD - &
                    I1invI3vector(1:NKe,1:NKe,dir,e)*khat(dir,i,j)
                end do
                tempC = tempD + Sigmat(eltomat(e))*ident
                tempD = explicit_4x4_matinv(tempC)
                Tprimeinv(1:NKe,1:NKe,enew(e,i,j),i,j) = transpose(tempD)

                !if (any(abs(matmul(Tprimeinv(1:NKe,1:NKe,enew(e,i,j),i,j),transpose(tempC)) - &
                !    ident) .gt. 1.0E-8)) then
                !    print *, [e,i,j]
                !end if
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine t2_ele_SN_MGXS_construct_inverted_transport_matrix

subroutine t2_ele_SN_MGXS_EI_SI &
    (particle, pmax, l_G, Cekk, bdyel, wL, khat, cosmat, qPa, qfactorialmat, esweeplist, &
    esweepbounds, enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, I1invI2f, &
    I1invI3vector, Sigma, Sigmat, phiu, source, fEg, bdysrc, psi, phi)
    implicit none
    integer, intent(in) :: particle
    integer, intent(in) :: pmax
    integer, intent(in) :: l_G
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:), intent(in) :: Sigmat
    real, dimension(:,:,:), intent(in) :: phiu
    real, dimension(:,:,:,:,:), allocatable, intent(in) :: source
    real, dimension(:), intent(in) :: fEg
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:), allocatable, intent(inout) :: phi

    integer :: e, k, g
    integer :: total
    real, dimension(:,:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:), allocatable :: t_bdysrc
    real, dimension(:,:,:), allocatable :: t_Sigma
    real, dimension(:,:), allocatable :: t_Sigmagg
    real, dimension(:), allocatable :: t_Sigmat
    real, dimension(:,:), allocatable :: t_phiu
    real, dimension(:,:,:,:), allocatable :: t_psi
    real, dimension(:,:,:,:,:), allocatable :: Tprimeinv
    real :: ompstart, ompfinish

    if (particle .eq. 1) then
        print *, "SN MGXS SI STARTED, PARTICLE: PHOTON"
    else if (particle .eq. 2) then
        print *, "SN MGXS SI STARTED, PARTICLE: ELECTRON"
    end if

    allocate(psi(NE,NKe,Nmu,Nphi,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(t_source(NE,NKe,Nmu,Nphi))
    if (allocated(bdysrc)) allocate(t_bdysrc(NKe,size(bdysrc,2),Nmu,Nphi))
    allocate(t_Sigmagg(size(Sigma,3),size(Sigma,4)))
    allocate(t_Sigmat(Nmats))
    allocate(t_phiu(NE,NKe))
    allocate(t_psi(NE,NKe,Nmu,Nphi))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

    do g = 1, l_G
        total = 0
        print *, "ENERGY GROUP:", g

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                t_source = source(:,:,:,:,1)
            else
                t_source = 0
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(1)*bdysrc
            end if
        else
            allocate(t_Sigma, source = Sigma(1:g-1,g,:,:))
            call t2_ele_SN_MGXS_K &
                (wL, cosmat, qPa, qfactorialmat, t_Sigma, psi, t_source)
            if (allocated(source)) then
                t_source = t_source + source(:,:,:,:,g)
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(g)*bdysrc
            end if
            deallocate(t_Sigma)
        end if
        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        t_Sigmagg = Sigma(g,g,:,:)
        t_Sigmat = Sigmat(g,:)
        t_phiu = phiu(:,:,g)

        ompstart = omp_get_wtime()
        call t2_ele_SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()

        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        call t2_ele_SN_MGXS_SI &
            (convergence, pmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, &
            esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, t_Sigmagg, t_source, t_bdysrc, Tprimeinv, total, t_phiu, t_psi)

        psi(:,:,:,:,g) = t_psi

        phi(:,:,g) = t_phiu

        ompfinish = omp_get_wtime()

        print *, "CONVERGED WITH:", total, "ITERATES"
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    phi = phi - phiu

    if (particle .eq. 1) then
        call t2_ele_SN_write_field("photon", psi, phi)
    else if (particle .eq. 2) then
        call t2_ele_SN_write_field("electron", psi, phi)
    end if
end subroutine t2_ele_SN_MGXS_EI_SI

subroutine t2_ele_SN_MGXS_SI &
    (l_convergence, pmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, &
    enew, fold, Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, source, bdysrc, Tprimeinv, total, &
    phi, psi)
    implicit none
    real, intent(in) :: l_convergence
    integer, intent(in) :: pmax
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2f ! (kp,k,f,e)
    real, dimension(:,:), intent(in) :: Sigma ! (l+1,mat)
    real, dimension(:,:,:,:), intent(inout) :: source ! (e,k,i,j)
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc
    real, dimension(:,:,:,:,:), intent(in) :: Tprimeinv ! (kp,k,enew,i,j)

    integer, intent(inout) :: total
    real, dimension(:,:), intent(inout) :: phi ! (e,k) ! phiu initially, leaves as phic + phiu
    real, dimension(:,:,:,:), intent(inout) :: psi ! (e,k,i,j)

    integer :: e, j, k, p
    real, dimension(:,:,:), allocatable :: jterm
    real, dimension(:,:), allocatable :: tempE
    real :: prevnorm
    real :: nextnorm
    real :: ompstart, ompfinish

    allocate(jterm(NE,NKe,Nphi))
    allocate(tempE(NE,NKe))

    ompstart = omp_get_wtime()

    prevnorm = norm2(phi)

    call t2_ele_SN_MGXS_sweep &
        (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
        I1invI2f, bdysrc, Tprimeinv, source)

    psi = source
    !$OMP parallel do shared(jterm) private(j, k)
    do j = 1, Nphi
        do k = 1, NKe
            jterm(1:NE,k,j) = twopi*matmul(source(1:NE,k,1:Nmu,j),wL)/Nphi
        end do
    end do
    !$OMP end parallel do
    tempE = sum(jterm,3)
    phi = phi + tempE

    nextnorm = norm2(tempE)

    total = 1

    ompfinish = omp_get_wtime()
    if (prevnorm .eq. 0) then
        print *, "ITERATE:", 1, "CONVERGENCE:", "Arbitrary", "TIME (s):", ompfinish - ompstart
    else
        print *, "ITERATE:", 1, "CONVERGENCE:", nextnorm/prevnorm, "TIME (s):", ompfinish - ompstart
        if (nextnorm/prevnorm .le. l_convergence) return
    end if

    do p = 2, pmax
        ompstart = omp_get_wtime()

        prevnorm = norm2(phi)

        call t2_ele_SN_MGXS_Kgg &
            (wL, cosmat, qPa, qfactorialmat, Sigma, source)

        call t2_ele_SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, bdysrc, Tprimeinv, source)

        psi = psi + source
        !$OMP parallel do shared(jterm) private(j, k)
        do j = 1, Nphi
            do k = 1, NKe
                jterm(1:NE,k,j) = twopi*matmul(source(1:NE,k,1:Nmu,j),wL)/Nphi
            end do
        end do
        !$OMP end parallel do
        tempE = sum(jterm,3)
        phi = phi + tempE

        nextnorm = norm2(tempE)

        total = p

        ompfinish = omp_get_wtime()
        print *, "ITERATE:", p, "CONVERGENCE:", nextnorm/prevnorm, "TIME (s):", ompfinish - ompstart
        if (nextnorm/prevnorm .le. l_convergence) return
        if (prevnorm .eq. 0.0) return ! Is this necessary given that this conditional was already invoked??
        if (isnan(nextnorm/prevnorm)) then
            print *, "NAN"
            stop
        end if
    end do

    print *, "CONVERGENCE NOT ACHIEVED WITH", pmax, "ITERATES."
end subroutine t2_ele_SN_MGXS_SI

subroutine t2_ele_SN_MGXS_K(wL, cosmat, qPa, qfactorialmat, Sigma, scatterers, vector)
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat ! (j,jp,m+1)
    real, dimension(:,:), intent(in) :: qPa ! (i,q)
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:), intent(in) :: Sigma ! (gp,l+1)
    real, dimension(:,:,:,:,:), intent(in) :: scatterers ! (e,k,i,j,g)

    real, dimension(:,:,:,:), intent(inout) :: vector ! (e,k,i,j)

    integer :: e, i, ip, j, jp, k, l, mat, q, g
    integer :: NQ
    integer :: group
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:,:), allocatable :: Amat ! (e,k,ip,jp,l)
    real, dimension(:,:,:,:), allocatable :: Bmat ! (e,k,jp,q)
    real, dimension(:,:,:,:), allocatable :: Cmat ! (e,k,jp,q)
    !real, dimension(Nphi) :: tempA

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    allocate(Amat(NE,NKe,Nmu,Nphi,NL+1))
    allocate(Bmat(NE,NKe,Nphi,NQ))
    allocate(Cmat(NE,NKe,Nphi,NQ))

    group = size(Sigma,1)
    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(4) shared(Amat) private(l, jp, ip, k)
        do l = 0, NL
            do jp = 1, Nphi
                do ip = 1, Nmu
                    do k = 1, NKe
                        Amat(1:NE,k,ip,jp,l+1) = twopi*matmul(scatterers(1:NE,k,ip,jp,1:group),Sigma(1:group,l+1,1))/Nphi
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do collapse(4) shared(Amat) private(e, l, ip, k)
        do l = 0, NL
            do ip = 1, Nmu
                do k = 1, NKe
                    do e = 1, NE
                        Amat(e,k,ip,1:Nphi,l+1) = twopi*matmul(scatterers(e,k,ip,1:Nphi,1:group),&
                            Sigma(1:group,l+1,eltomat(e)))/Nphi
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    end if

    !$OMP parallel do shared(Bmat) private(k, jp, q)
    do q = 1, NQ
        do jp = 1, Nphi
            do k = 1, NKe
                Bmat(1:NE,k,jp,q) = matmul(Amat(1:NE,k,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
            end do
        end do
        Bmat(1:NE,1:NKe,1:Nphi,q) = qfactorialmat(q)*Bmat(1:NE,1:NKe,1:Nphi,q)
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(2) shared(Cmat) private(k, q)
    do q = 1, NQ
        do k = 1, NKe
            Cmat(1:NE,k,1:Nphi,q) = matmul(Bmat(1:NE,k,1:Nphi,q),cosmat(1:Nphi,1:Nphi,emm(q)+1))
        end do
    end do
    !$OMP end parallel do

    if (size(Sigma,2) .eq. NL+2) then
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(3) shared(vector) private(i, j, k)
            do j = 1, Nphi
                do i = 1, Nmu
                    do k = 1, NKe
                        vector(1:NE,k,i,j) = matmul(scatterers(1:NE,k,i,j,1:group),&
                            Sigma(1:group,NL+2,1))
                    end do
                end do
            end do
            !$OMP end parallel do
        else
            !$OMP parallel do collapse(3) shared(vector) private(e, i, j, k)
            do i = 1, Nmu
                do k = 1, NKe
                    do e = 1, NE
                        vector(e,k,i,1:Nphi) = matmul(scatterers(e,k,i,1:Nphi,1:group),&
                            Sigma(1:group,NL+2,eltomat(e)))
                    end do
                end do
            end do
            !$OMP end parallel do
        end if
    else
        vector = 0
    end if

    !$OMP parallel do shared(vector) private(i, q)
    do i = 1, Nmu
        do q = 1, NQ
            vector(1:NE,1:NKe,i,1:Nphi) = vector(1:NE,1:NKe,i,1:Nphi) + &
                Cmat(1:NE,1:NKe,1:Nphi,q)*qPa(i,q)
        end do
    end do
    !$OMP end parallel do
end subroutine t2_ele_SN_MGXS_K

subroutine t2_ele_SN_MGXS_Kgg(wL, cosmat, qPa, qfactorialmat, Sigma, vector)
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat ! (j,jp,m+1)
    real, dimension(:,:), intent(in) :: qPa ! (i,q)
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:), intent(in) :: Sigma ! (l)

    real, dimension(:,:,:,:), intent(inout) :: vector ! (e,k,i,j)

    integer :: e, i, ip, j, jp, k, l, mat, q
    integer :: NQ
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:,:), allocatable :: Amat ! (e,k,i,j,l)
    real, dimension(:,:,:,:), allocatable :: Bmat ! (e,k,j,q)
    real, dimension(:,:,:,:), allocatable :: Cmat ! (e,k,j,q)
    real, dimension(Nphi) :: tempA

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))

    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    allocate(Amat(NE,NKe,Nmu,Nphi,NL+1))
    if (Nmats .eq. 1) then
        !$OMP parallel do shared(Amat) private(l)
        do l = 0, NL
            Amat(1:NE,1:NKe,1:Nmu,1:Nphi,l+1) = twopi*Sigma(l+1,1)*vector/Nphi
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do collapse(2) shared(Amat) private(e, l)
        do l = 0, NL
            do e = 1, NE
                Amat(e,1:NKe,1:Nmu,1:Nphi,l+1) = twopi*Sigma(l+1,eltomat(e))*&
                    vector(e,1:NKe,1:Nmu,1:Nphi)/Nphi
            end do
        end do
        !$OMP end parallel do
    end if

    allocate(Bmat(NE,NKe,Nphi,NQ))
    !$OMP parallel do shared(Bmat) private(jp, k, q)
    do q = 1, NQ
        do jp = 1, Nphi
            do k = 1, NKe
                Bmat(1:NE,k,jp,q) = matmul(Amat(1:NE,k,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
            end do
        end do
        Bmat(1:NE,1:NKe,1:Nphi,q) = qfactorialmat(q)*Bmat(1:NE,1:NKe,1:Nphi,q)
    end do
    !$OMP end parallel do
    deallocate(Amat)

    allocate(Cmat(NE,NKe,Nphi,NQ))
    !$OMP parallel do collapse(2) shared(Cmat) private(k, q)
    do q = 1, NQ
        do k = 1, NKe
            Cmat(1:NE,k,1:Nphi,q) = matmul(Bmat(1:NE,k,1:Nphi,q),cosmat(1:Nphi,1:Nphi,emm(q)+1))
        end do
    end do
    !$OMP end parallel do
    deallocate(Bmat)

    vector = 0
    !$OMP parallel do shared(vector) private(i, q) ! collapse(1) this?
    do i = 1, Nmu
        do q = 1, NQ
            vector(1:NE,1:NKe,i,1:Nphi) = vector(1:NE,1:NKe,i,1:Nphi) + &
                Cmat(1:NE,1:NKe,1:Nphi,q)*qPa(i,q)
        end do
    end do
    !$OMP end parallel do
end subroutine t2_ele_SN_MGXS_Kgg

subroutine t2_ele_SN_MGXS_sweep &
    (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
    I1invI2f, bdysrc, Tprimeinv, vector)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc
    real, dimension(:,:,:,:,:), intent(in) :: Tprimeinv

    real, dimension(:,:,:,:), intent(inout) :: vector

    integer :: e, f, i, ip, j, k
    integer, dimension(:), allocatable :: l_esweeplist
    integer, dimension(:), allocatable :: l_enew
    integer, dimension(:,:), allocatable :: l_fold
    integer, dimension(:), allocatable :: l_Ia
    integer, dimension(:), allocatable :: l_Ja
    real, dimension(:), allocatable :: l_varsigmadown
    real, dimension(NKe) :: tempA
    real, dimension(:,:), allocatable :: l_vector
    integer :: SS1
    logical :: use_bdy
    integer :: NBE

    allocate(l_esweeplist(NE))
    allocate(l_enew(NE))
    allocate(l_fold(NFe,NE))
    allocate(l_Ia(NE+1))
    allocate(l_Ja(noJa))
    allocate(l_varsigmadown(noJa))
    allocate(l_vector(NKe,NE))

    use_bdy = allocated(bdysrc)
    if (use_bdy) then
        NBE = size(bdysrc,2)
    end if

    !$OMP parallel do collapse(2) shared(vector) private(j, i, k, ip, f, e) &
    !$OMP private(l_esweeplist, l_enew, l_vector, l_fold, l_Ia, l_Ja, l_varsigmadown, SS1, tempA)
    do j = 1, Nphi
        do i = 1, Nmu
            l_esweeplist = esweeplist(1:NE,i,j) ! Don't use local arrays with parallelization?
            l_enew = enew(1:NE,i,j)

            do k = 1, NKe
                do e = 1, NE
                    l_vector(k,l_enew(e)) = vector(e,k,i,j)
                end do
            end do

            if (use_bdy) then
                do e = 1, NBE
                    l_vector(1:NKe,l_enew(bdyel(e))) = l_vector(1:NKe,l_enew(bdyel(e))) + &
                        bdysrc(1:NKe,e,i,j)
                end do
            end if

            l_fold = fold(1:NFe,1:NE,i,j)
            l_Ia = Ia(1:NE+1,i,j) ! Why isn't the NE+1 term ever used? I forgot
            l_Ja = Ja(1:noJa,i,j)
            l_varsigmadown = varsigmadown(1:noJa,i,j)

            SS1 = esweepbounds(2,i,j)
            do ip = 1, SS1
                tempA = l_vector(1:NKe,ip)
                do k = 1, NKe
                    l_vector(k,ip) = sum(Tprimeinv(1:NKe,k,ip,i,j)*tempA)
                end do
            end do

            do ip = SS1+1, NE
                tempA = l_vector(1:NKe,ip)
                do f = 1, l_Ia(ip+1)-l_Ia(ip)
                    do k = 1, NKe
                        tempA(k) = tempA(k) - l_varsigmadown(l_Ia(ip)-1+f)*&
                        sum(I1invI2f(1:NKe,k,l_fold(f,ip),l_esweeplist(ip))*l_vector(1:NKe,l_Ja(l_Ia(ip)-1+f))) ! Make this faster by using a temporary for l_vector(eprime)? I can know its size before f iteration
                        ! And/or use l_I1invI2f to avoid repeated indirect addressing??
                    end do
                end do
                do k = 1, NKe
                    l_vector(k,ip) = sum(Tprimeinv(1:NKe,k,ip,i,j)*tempA)
                end do
            end do

            do e = 1, NE
                vector(e,1:NKe,i,j) = l_vector(1:NKe,l_enew(e))
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine t2_ele_SN_MGXS_sweep

subroutine t2_ele_SN_MGXS_EI_GMRESm &
    (particle, pmax, mmax, l_G, Cekk, bdyel, wL, khat, cosmat, qPa, qfactorialmat, esweeplist, &
    esweepbounds, enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, I1invI2f, &
    I1invI3vector, Sigma, Sigmat, source, fEg, bdysrc, psi, phi)
    implicit none
    integer, intent(in) :: particle
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    integer, intent(in) :: l_G
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:), intent(in) :: Sigmat
    real, dimension(:,:,:,:,:), allocatable, intent(in) :: source
    real, dimension(:), intent(in) :: fEg
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:), allocatable, intent(inout) :: phi

    integer :: e, g, i, j, k
    integer :: total
    real, dimension(:,:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:), allocatable :: t_bdysrc
    real, dimension(:,:,:), allocatable :: t_Sigma
    real, dimension(:,:), allocatable :: t_Sigmagg
    real, dimension(:), allocatable :: t_Sigmat
    real, dimension(:,:,:,:), allocatable :: t_psi
    real, dimension(:,:,:,:,:), allocatable :: Tprimeinv
    real :: ompstart, ompfinish

    if (particle .eq. 1) then
        print *, "SN MGXS GMRESm STARTED, PARTICLE: PHOTON"
    else if (particle .eq. 2) then
        print *, "SN MGXS GMRESm STARTED, PARTICLE: ELECTRON"
    end if

    allocate(psi(NE,NKe,Nmu,Nphi,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(t_source(NE,NKe,Nmu,Nphi))
    if (allocated(bdysrc)) allocate(t_bdysrc(NKe,size(bdysrc,2),Nmu,Nphi))
    allocate(t_Sigmagg(size(Sigma,3),size(Sigma,4)))
    allocate(t_Sigmat(Nmats))
    allocate(t_psi(NE,NKe,Nmu,Nphi))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

    do g = 1, l_G
        total = 0
        print *, "ENERGY GROUP:", g

        ompstart = omp_get_wtime()
        t_Sigmagg = Sigma(g,g,:,:)
        t_Sigmat = Sigmat(g,:)

        call t2_ele_SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()
        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                t_source = source(:,:,:,:,1)
            else
                t_source = 0
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(1)*bdysrc
            end if
        else
            allocate(t_Sigma, source = Sigma(1:g-1,g,:,:))
            call t2_ele_SN_MGXS_K &
                (wL, cosmat, qPa, qfactorialmat, t_Sigma, psi, t_source)
            if (allocated(source)) then
                t_source = t_source + source(:,:,:,:,g)
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(g)*bdysrc
            end if
            deallocate(t_Sigma)
        end if

        call t2_ele_SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
            varsigmadown, I1invI2f, t_bdysrc, Tprimeinv, t_source)

        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        t_psi = t_source ! GUESS
        call t2_ele_SN_MGXS_GMRESm &
            (convergence, pmax, mmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, &
            esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, t_Sigmagg, t_source, t_bdysrc, Tprimeinv, total, t_psi)

        psi(:,:,:,:,g) = t_psi

        print *, "CONVERGED WITH:", total, "ITERATES"

        ompfinish = omp_get_wtime()
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    phi = 0
    do g = 1, l_G
        do k = 1, NKe
            do j = 1, Nphi
                phi(1:NE,k,g) = phi(1:NE,k,g) + twopi*matmul(psi(1:NE,k,1:Nmu,j,g),wL)/Nphi
            end do
        end do
    end do

    if (particle .eq. 1) then
        call t2_ele_SN_write_field("photon", psi, phi)
    else if (particle .eq. 2) then
        call t2_ele_SN_write_field("electron", psi, phi)
    end if
end subroutine t2_ele_SN_MGXS_EI_GMRESm

subroutine t2_ele_SN_MGXS_GMRESm &
    (l_convergence, pmax, mmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, &
    esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, source, bdysrc, &
    Tprimeinv, total, psi)
    implicit none
    real, intent(in) :: l_convergence
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2f ! (kp,k,f,e)
    real, dimension(:,:), intent(in) :: Sigma ! (l+1,mat)
    real, dimension(:,:,:,:), intent(in) :: source ! (e,k,i,j)
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc
    real, dimension(:,:,:,:,:), intent(in) :: Tprimeinv ! (kp,k,enew,i,j)

    integer, intent(inout) :: total
    real, dimension(:,:,:,:), intent(inout) :: psi ! (e,k,i,j)

    integer :: i, j
    integer :: m
    integer :: max_iter
    real :: resid
    real, dimension(:,:,:,:), allocatable :: tempA
    real, dimension(:,:,:,:), allocatable :: tempB
    real, dimension(:,:,:,:), allocatable :: tempD
    real, dimension(:,:,:,:), allocatable :: resid0
    real, dimension(:,:,:,:,:), allocatable :: vec ! Name it V?
    real, dimension(:,:), allocatable :: Qk
    real, dimension(:,:), allocatable :: local_I
    real, dimension(:,:), allocatable :: Hess
    real, dimension(:,:), allocatable :: Rk
    real, dimension(:,:), allocatable :: Rot
    real, dimension(:,:), allocatable :: tempC
    real, dimension(:), allocatable :: gk
    real, dimension(:), allocatable :: tvec
    real, dimension(:), allocatable :: e1
    real, dimension(:), allocatable :: y
    real :: denom
    real :: ompstart, ompfinish

    m = pmax
    max_iter = mmax

    allocate(tempA(NE,NKe,Nmu,Nphi))
    allocate(tempB(NE,NKe,Nmu,Nphi))
    allocate(tempD(NE,NKe,Nmu,Nphi))
    allocate(resid0(NE,NKe,Nmu,Nphi))
    allocate(vec(NE,NKe,Nmu,Nphi,m))
    allocate(Qk(m+1,m+1))
    allocate(local_I(m+1,m+1))
    allocate(Hess(m+1,m))
    allocate(Rk(m+1,m))
    allocate(Rot(2,2))
    allocate(tempC(2,m+1))
    allocate(gk(m+1))
    allocate(tvec(m+1))
    allocate(e1(m+1))
    allocate(y(m))

    local_I = identity(m+1)
    e1 = 0
    e1(1) = 1

    tempA = psi

    call t2_ele_SN_MGXS_swept_Lgg &
        (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
        Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempD, tempA)

    resid0 = source - tempA

    GMRESm: do
        resid = norm2(resid0)

        if (resid .lt. l_convergence) then
            print *, total, resid
            exit GMRESm
        end if

        vec(1:NE,1:NKe,1:Nmu,1:Nphi,1) = resid0/resid

        Qk = local_I

        do j = 1, m-1
            ompstart = omp_get_wtime()
            total = total + 1
            ! Gram-Schmidt orthogonalization to construct vec(:,j) vectors and Hess array
            tempA = vec(1:NE,1:NKe,1:Nmu,1:Nphi,j)

            call t2_ele_SN_MGXS_swept_Lgg &
                (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempD, tempA)

            tempB = tempA
            !$OMP parallel do shared(Hess, tempB) private(i)
            do i = 1, j
                Hess(i,j) = sum(tempA*vec(1:NE,1:NKe,1:Nmu,1:Nphi,i))
            end do
            !$OMP end parallel do
            do i = 1, j
                tempB = tempB - Hess(i,j)*vec(1:NE,1:NKe,1:Nmu,1:Nphi,i)
            end do

            Hess(j+1,j) = norm2(tempB)
            vec(1:NE,1:NKe,1:Nmu,1:Nphi,i) = tempB/Hess(j+1,j)

            ! Append latest updated column of Hess to R, with cumulative rotations applied
            Rk(j+1,j) = Hess(j+1,j)
            Rk(1:j,j) = matmul(Qk(1:j,1:j),Hess(1:j,j)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

            ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
            denom = hypot(Rk(j,j),Rk(j+1,j))
            Rot(:,1) = [Rk(j,j)/denom, -Rk(j+1,j)/denom]
            Rot(:,2) = [-Rot(2,1), Rot(1,1)]

            Rk(j:j+1,j) = [denom, 0.0]

            ! Update cumulative rotation matrix
            tempC(1:2,1:j+1) = Qk(j:j+1,1:j+1)
            Qk(j:j+1,1:j+1) = matmul(Rot,tempC(1:2,1:j+1))

            ! Check for convergence
            ompfinish = omp_get_wtime()
            print *, "ITERATE:", j, "CONVERGENCE:", abs(resid*Qk(j+1,1)), &
                "TIME (s):", ompfinish - ompstart
            if (abs(resid*Qk(j+1,1)) <= l_convergence) then
                gk(1:j+1) = resid*Qk(1:j+1,1)
                do i = j, 1, -1
                    y(i) = (gk(i) - sum(Rk(i,i+1:j)*y(i+1:j)))/Rk(i,i)
                end do

                do i = 1, j ! THIS STEP MAY BE EASIER IF RESHAPES AND PACKS ARE USED
                    psi = psi + vec(1:NE,1:NKe,1:Nmu,1:Nphi,i)*y(i)
                end do

                exit GMRESm
            end if
            if (isnan(abs(resid*Qk(j+1,1)))) then
                print *, "NAN"
                stop
            end if
        end do

        total = total + 1

        ! Do mth step differently
        tempA = vec(1:NE,1:NKe,1:Nmu,1:Nphi,m)

        call t2_ele_SN_MGXS_swept_Lgg &
            (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
            Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempD, tempA)

        !$OMP parallel do shared(Hess) private(i)
        do i = 1, m
            Hess(i,m) = sum(tempA*vec(1:NE,1:NKe,1:Nmu,1:Nphi,i))
        end do
        !$OMP end parallel do
        !Hess(m+1,m) = sqrt(norm2(tempA)**2 - sum(Hess(1:m,m)**2)) ! Maybe relies on elemental??
        tempB = 0
        do i = 1, m
            tempB = tempB - Hess(i,m)*vec(1:NE,1:NKe,1:Nmu,1:Nphi,i)
        end do
        Hess(m+1,m) = norm2(tempA - tempB)

        ! Append latest updated column of Hess to R, with cumulative rotations applied
        Rk(m+1,m) = Hess(m+1,m)
        Rk(1:m,m) = matmul(Qk(1:m,1:m),Hess(1:m,m)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

        ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
        denom = hypot(Rk(m,m),Rk(m+1,m))
        Rot(:,1) = [Rk(m,m)/denom, -Rk(m+1,m)/denom] ! Should I transpose these? were they supposed to be transposed? (~3/31/23)
        Rot(:,2) = [-Rot(2,1), Rot(1,1)]

        Rk(m:m+1,m) = [denom, 0.0]

        ! Update cumulative rotation matrix
        tempC(1:2,1:m+1) = Qk(m:m+1,1:m+1)
        Qk(m:m+1,1:m+1) = matmul(Rot,tempC(1:2,1:m+1))

        ! Construct current solution
        gk(1:m+1) = resid*Qk(1:m+1,1)
        do i = m, 1, -1
            y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
        end do

        do i = 1, m
            psi = psi + vec(1:NE,1:NKe,1:Nmu,1:Nphi,i)*y(i)
        end do

        tempA = psi

        call t2_ele_SN_MGXS_swept_Lgg &
            (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
            Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempD, tempA)

        resid0 = source - tempA

        !! Slower?? Doesn't work??
        !! Calculate residual for new vector without computing f - Ax0
        !tvec = resid*e1 - matmul(Hess,y)
        !
        !resid0 = tvec(m+1)*tempA/Hess(m+1,m)
        !!!$OMP parallel do shared(resid0) private(i)
        !do i = 1, m
        !    resid0 = resid0 + &
        !    vec(1:NE,1:NKe,1:NQ,i)*(tvec(i)-tvec(m+1)*Hess(i,m)/Hess(m+1,m))
        !end do
        !!!$OMP end parallel do
        print *, "ITERATE:", m, "CONVERGENCE:", abs(gk(m+1)), &
            "TIME (s):", ompfinish - ompstart
        if (abs(gk(m+1)) .le. convergence) then
            exit GMRESm
        end if
        if (isnan(abs(gk(m+1)))) then
            print *, "NAN"
            stop
        end if
        if (total .ge. max_iter) then
            gk(1:m+1) = resid*Qk(1:m+1,1)
            do i = m, 1, -1
                y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
            end do

            do i = 1, m
                psi = psi + vec(1:NE,1:NKe,1:Nmu,1:Nphi,i)*y(i)
            end do
            print *, "GMRESm: Did not converge within", max_iter, "iterations."
            print *, "Solution being given with:", abs(gk(j+1)), "residual."
            exit GMRESm
        end if

        print *, "GMRESm: RESTARTING"
    end do GMRESm
end subroutine t2_ele_SN_MGXS_GMRESm

subroutine t2_ele_SN_MGXS_swept_Lgg &
    (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, Ia, Ja, &
    noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempA, vector)
    implicit none
    ! NON-ELEMENTALNESS GONNA BE A PROBLEM?? SURE HOPE NOT.
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:), intent(in) :: esweepbounds
    integer, dimension(:,:,:), intent(in) :: enew
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, dimension(:,:,:), intent(in) :: Ja
    integer, intent(in) :: noJa
    real, dimension(:,:,:), intent(in) :: varsigmadown
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:,:), intent(in) :: Tprimeinv

    real, dimension(:,:,:,:), intent(inout) :: tempA ! Pre-allocated temporary array, saves time. Make sure I don't forget it if I do MMS stuff.
    real, dimension(:,:,:,:), intent(inout) :: vector

    real, dimension(:,:,:,:), allocatable :: nobdysrc

    tempA = vector

    call t2_ele_SN_MGXS_Kgg(wL, cosmat, qPa, qfactorialmat, Sigma, tempA)

    call t2_ele_SN_MGXS_sweep &
        (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
        varsigmadown, I1invI2f, nobdysrc, Tprimeinv, tempA)

    vector = vector - tempA
end subroutine t2_ele_SN_MGXS_swept_Lgg

! MISC
!!! wiscoslab stuff
subroutine wiscoslab_SN_MGXS_construct_inverted_transport_matrix(mu, GpFterm, Sigmat, Tprimeinv)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:,:,:,:), intent(in) :: GpFterm
    real, dimension(:), intent(in) :: Sigmat

    real, dimension(:,:,:,:), intent(inout) :: Tprimeinv ! (k,kp,e,i)

    integer :: e, i, k, kp
    real, dimension(:,:), allocatable :: ident
    real, dimension(:,:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: tempB
    real :: factor
    integer :: estart
    integer :: efin
    integer :: einc

    allocate(ident(NKe,NKe))
    allocate(tempA(NKe,NKe))
    allocate(tempB(NKe,NKe,NE))

    ident = identity(NKe)

    do i = 1, Nmu
        do e = 1, NE
            Tprimeinv(1:NKe,1:NKe,e,i) = GpFterm(1:NKe,1:NKe,e,i) + Sigmat(eltomat(e))*ident
        end do
    end do

    do i = 1, Nmu
        do e = 1, NE
            tempA = Tprimeinv(1:NKe,1:NKe,e,i)
            factor = tempA(1,1)*tempA(2,2) - tempA(2,1)*tempA(1,2)
            ! (1,1)
            Tprimeinv(1,1,e,i) = tempA(2,2)/factor
            ! (1,2)
            Tprimeinv(1,2,e,i) = -tempA(1,2)/factor
            ! (2,1)
            Tprimeinv(2,1,e,i) = -tempA(2,1)/factor
            ! (2,2)
            Tprimeinv(2,2,e,i) = tempA(1,1)/factor
        end do
    end do

    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            estart = NE
            efin = 1 ! Make it NE-estart+1
            einc = -1 ! Make it sign(1,int(efin-estart)) or something
        else
            estart = 1
            efin = NE
            einc = 1
        end if

        ! Sweep order stuff
        tempB = Tprimeinv(1:NKe,1:NKe,1:NE,i)
        Tprimeinv(1:NKe,1:NKe,1:NE,i) = tempB(1:NKe,1:NKe,estart:efin:einc)
    end do
end subroutine wiscoslab_SN_MGXS_construct_inverted_transport_matrix

subroutine wiscoslab_SN_MGXS_sweep(mu, Fsweep, bdy, Tprimeinv, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:,:), allocatable, intent(in) :: bdy
    real, dimension(:,:,:,:), intent(in) :: Tprimeinv

    real, dimension(:,:,:), intent(inout) :: vector ! (e,k,i)

    integer :: i, ip, j, k
    integer :: estart
    integer :: efin
    integer :: einc
    integer :: SS1
    real, dimension(:,:,:), allocatable :: l_Fsweep
    real, dimension(:,:), allocatable :: l_vector
    real, dimension(:), allocatable :: tempA
    logical :: do_bdy

    allocate(l_Fsweep(NKe,NKe,NE))
    allocate(l_vector(NKe,NE))
    allocate(tempA(NKe))

    do_bdy = allocated(bdy)

    !$OMP parallel do shared(vector) private(i, ip) private(estart, efin, einc, l_Fsweep, l_vector) &
    !$OMP private(SS1, tempA)
    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            estart = NE
            efin = 1 ! Make it NE-estart+1
            einc = -1 ! Make it sign(1,int(efin-estart)) or something
            l_Fsweep = Fsweep(1:NKe,1:NKe,1:NE,i,2)
        else
            estart = 1
            efin = NE
            einc = 1
            l_Fsweep = Fsweep(1:NKe,1:NKe,1:NE,i,1)
        end if

        l_vector(1:NKe,1:NE) = transpose(vector(estart:efin:einc,1:NKe,i))

        if (do_bdy) then
            !!! USE THIS STRUCTURE IN 3D. NOT USING HERE BECAUSE SITUATION IS MUCH SIMPLER
            !!! MUST USE THIS BECAUSE SOME BOUNDARY TERMS AREN'T FIRST TERM IN SWEEP.
            !do ip = 1, size(bdy,2)
            !    l_vector(1:NKe,bdyupstr(ip)) = l_vector(1,NKe,bdyupstr(ip)) - bdysrc(1:NKe,ip,i)
            !end do

            l_vector(1:NKe,1) = l_vector(1:NKe,1) - &
                merge(bdy(1:NKe,2,i), bdy(1:NKe,1,i), mu(i) .gt. 0)
        end if

        SS1 = 1
        do ip = 1, SS1
            tempA = l_vector(1:NKe,ip)

            l_vector(1:NKe,ip) = matmul(Tprimeinv(1:NKe,1:NKe,ip,i),tempA)
        end do

        do ip = SS1+1, NE
            tempA = l_vector(1:NKe,ip) - matmul(l_Fsweep(1:NKe,1:NKe,ip),l_vector(1:NKe,ip-1))

            l_vector(1:NKe,ip) = matmul(Tprimeinv(1:NKe,1:NKe,ip,i),tempA)
        end do

        vector(estart:efin:einc,1:NKe,i) = transpose(l_vector(1:NKe,1:NE))
    end do
    !$OMP end parallel do
end subroutine wiscoslab_SN_MGXS_sweep

subroutine ele_SN_MGXS_T &
    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, Sigmat, &
    bdysrc, vector)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:,:), intent(in) :: eprime
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:,:,:), intent(in) :: normal
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:), intent(in) :: Sigmat
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc

    real, dimension(:,:,:,:), intent(inout) :: vector

    integer :: dir, e, f, i, ip, j, k, kp
    real, dimension(:,:,:,:), allocatable :: Mterm
    real, dimension(:,:,:,:,:), allocatable :: Tprime
    real, dimension(:,:), allocatable :: tempD
    real, dimension(:,:), allocatable :: ident
    integer :: l_e
    real :: l_varsigmadown
    logical :: use_bdy

    allocate(Mterm(NE,NKe,Nmu,Nphi))
    allocate(Tprime(NKe,NKe,NE,Nmu,Nphi))
    allocate(tempD(NKe,NKe))
    allocate(ident(NKe,NKe))

    use_bdy = allocated(bdysrc)

    ident = identity(NKe)

    !$OMP parallel do collapse(3) shared(Tprime) private(j, i, e, f, dir) &
    !$OMP private(tempD)
    do j = 1, Nphi
        do i = 1, Nmu
            do e = 1, NE
                tempD = 0
                do f = 1, NFe
                    tempD = tempD + &
                    varsigmaup(f,e,i,j)*I1invI2(1:NKe,1:NKe,f,e)
                end do
                do dir = 1, 3
                    tempD = tempD - &
                    I1invI3vector(1:NKe,1:NKe,dir,e)*khat(dir,i,j)
                end do
                Tprime(1:NKe,1:NKe,e,i,j) = tempD + Sigmat(eltomat(e))*ident
            end do
        end do
    end do
    !$OMP end parallel do

    Mterm = 0
    !$OMP parallel do collapse(3) shared(Mterm) private(j, i, e)
    do j = 1, Nphi
        do i = 1, Nmu
            do e = 1, NE
                Mterm(e,1:NKe,i,j) = matmul(Tprime(1:NKe,1:NKe,e,i,j),vector(e,1:NKe,i,j))
            end do
        end do
    end do
    !$OMP end parallel do

    if (use_bdy) then
        do ip = 1, size(bdyel)
            Mterm(bdyel(ip),1:NKe,1:Nmu,1:Nphi) = Mterm(bdyel(ip),1:NKe,1:Nmu,1:Nphi) - & ! Uses minus because bdysrc is on RHS generally
                bdysrc(1:NKe,ip,1:Nmu,1:Nphi)
        end do
    end if

    !!$OMP parallel do collapse(2) shared(Mterm) private(f, e, j, i, k, kp) private(l_varsigmadown)
    do f = 1, NFe
        do e = 1, NE
            l_e = eprime(f,e)
            if (l_e .eq. 0) cycle
            do j = 1, Nphi
                do i = 1, Nmu
                    l_varsigmadown = dot_product(khat(1:3,i,j),normal(1:3,f,e))
                    if (l_varsigmadown .ge. 0) cycle
                    do k = 1, NKe
                        do kp = 1, NKe
                            Mterm(e,k,i,j) = Mterm(e,k,i,j) + &
                                l_varsigmadown*I1invI2f(kp,k,f,e)*vector(l_e,kp,i,j)
                        end do
                    end do
                end do
            end do
        end do
    end do
    !!$OMP end parallel do

    vector = Mterm
end subroutine ele_SN_MGXS_T

subroutine wiscoslab_SN_MGXS_EI_GMRESm &
    (particle, pmax, mmax, l_G, mu, wL, GpFterm, Fsweep, Sigma, Sigmat, source, bdy, fEg, psi, phi)
    implicit none
    integer, intent(in) :: particle
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    integer, intent(in) :: l_G
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:), intent(in) :: GpFterm
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:), intent(in) :: Sigmat
    real, dimension(:,:,:,:), allocatable, intent(in) :: source
    real, dimension(:,:,:), allocatable, intent(in) :: bdy
    real, dimension(:), intent(in) :: fEg

    real, dimension(:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:), allocatable, intent(inout) :: phi

    integer :: e, g, i, j, k
    integer :: total
    real, dimension(:,:,:), allocatable :: t_source
    real, dimension(:,:,:), allocatable :: t_Sigma
    real, dimension(:,:), allocatable :: t_Sigmagg
    real, dimension(:), allocatable :: t_Sigmat
    real, dimension(:,:,:), allocatable :: t_psi
    real, dimension(:,:,:,:), allocatable :: Tprimeinv
    real, dimension(:,:,:), allocatable :: l_bdy
    real :: ompstart, ompfinish

    if (particle .eq. 1) then
        print *, "wiscoslab: SN MGXS GMRESm STARTED, PARTICLE: PHOTON"
    else if (particle .eq. 2) then
        print *, "wiscoslab: SN MGXS GMRESm STARTED, PARTICLE: ELECTRON"
    end if

    allocate(psi(NE,NKe,Nmu,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(t_source(NE,NKe,Nmu))
    allocate(t_Sigmagg(size(Sigma,3),size(Sigma,4)))
    allocate(t_Sigmat(Nmats))
    allocate(t_psi(NE,NKe,Nmu))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu))
    if (allocated(bdy)) allocate(l_bdy(size(bdy,1),size(bdy,2),size(bdy,3)))

    do g = 1, l_G
        total = 0
        print *, "ENERGY GROUP:", g

        ompstart = omp_get_wtime()
        t_Sigmagg = Sigma(g,g,:,:)
        t_Sigmat = Sigmat(g,:)

        call wiscoslab_SN_MGXS_construct_inverted_transport_matrix &
            (mu, GpFterm, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()
        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                t_source = source(:,:,:,1)
            else
                t_source = 0
            end if
            if (allocated(bdy)) then
                l_bdy = fEg(1)*bdy
            end if
        else
            allocate(t_Sigma, source = Sigma(1:g-1,g,:,:))
            call wiscoslab_SN_MGXS_K &
                (mu, wL, t_Sigma, psi, t_source)
            if (allocated(source)) then
                t_source = t_source + source(:,:,:,g)
            end if
            if (allocated(bdy)) then
                l_bdy = fEg(g)*bdy
            end if
            deallocate(t_Sigma)
        end if

        call wiscoslab_SN_MGXS_sweep(mu, Fsweep, l_bdy, Tprimeinv, t_source)

        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        t_psi = t_source ! GUESS
        call wiscoslab_SN_MGXS_GMRESm &
            (convergence, pmax, mmax, mu, wL, &
            Fsweep, t_Sigmagg, Tprimeinv, t_source, l_bdy, total, t_psi)

        psi(:,:,:,g) = t_psi

        print *, "CONVERGED WITH:", total, "ITERATES"

        ompfinish = omp_get_wtime()
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    do g = 1, l_G
        do k = 1, NKe
            phi(1:NE,k,g) = twopi*matmul(psi(1:NE,k,1:Nmu,g),wL)
        end do
    end do
end subroutine wiscoslab_SN_MGXS_EI_GMRESm

subroutine wiscoslab_SN_MGXS_GMRESm &
    (l_convergence, pmax, mmax, mu, wL, Fsweep, Sigma, Tprimeinv, source, bdy, total, psi)
    implicit none
    real, intent(in) :: l_convergence
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:), intent(in) :: Sigma ! (l+1,mat)
    real, dimension(:,:,:,:), intent(in) :: Tprimeinv ! (kp,k,ip,i)
    real, dimension(:,:,:), intent(in) :: source ! (e,k,i)
    real, dimension(:,:,:), allocatable, intent(in) :: bdy

    integer, intent(inout) :: total
    real, dimension(:,:,:), intent(inout) :: psi ! (e,k,i)

    integer :: i, j
    integer :: m
    integer :: max_iter
    real :: resid
    real, dimension(:,:,:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: tempB
    real, dimension(:,:,:), allocatable :: tempD
    real, dimension(:,:,:), allocatable :: resid0
    real, dimension(:,:,:,:), allocatable :: vec ! Name it V?
    real, dimension(:,:), allocatable :: Qk
    real, dimension(:,:), allocatable :: local_I
    real, dimension(:,:), allocatable :: Hess
    real, dimension(:,:), allocatable :: Rk
    real, dimension(:,:), allocatable :: Rot
    real, dimension(:,:), allocatable :: tempC
    real, dimension(:), allocatable :: gk
    real, dimension(:), allocatable :: tvec
    real, dimension(:), allocatable :: e1
    real, dimension(:), allocatable :: y
    real :: denom
    real :: ompstart, ompfinish

    m = pmax
    max_iter = mmax

    allocate(tempA(NE,NKe,Nmu))
    allocate(tempB(NE,NKe,Nmu))
    allocate(tempD(NE,NKe,Nmu))
    allocate(resid0(NE,NKe,Nmu))
    allocate(vec(NE,NKe,Nmu,m))
    allocate(Qk(m+1,m+1))
    allocate(local_I(m+1,m+1))
    allocate(Hess(m+1,m))
    allocate(Rk(m+1,m))
    allocate(Rot(2,2))
    allocate(tempC(2,m+1))
    allocate(gk(m+1))
    allocate(tvec(m+1))
    allocate(e1(m+1))
    allocate(y(m))

    local_I = identity(m+1)
    e1 = 0
    e1(1) = 1

    tempA = psi

    call wiscoslab_SN_MGXS_swept_Lgg &
        (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

    resid0 = source - tempA

    GMRESm: do
        resid = norm2(resid0)

        if (resid .lt. l_convergence) then
            print *, total, resid
            exit GMRESm
        end if

        vec(1:NE,1:NKe,1:Nmu,1) = resid0/resid

        Qk = local_I

        do j = 1, m-1
            ompstart = omp_get_wtime()
            total = total + 1
            ! Gram-Schmidt orthogonalization to construct vec(:,j) vectors and Hess array
            tempA = vec(1:NE,1:NKe,1:Nmu,j)

            call wiscoslab_SN_MGXS_swept_Lgg &
                (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

            tempB = tempA
            !$OMP parallel do shared(Hess, tempB) private(i)
            do i = 1, j
                Hess(i,j) = sum(tempA*vec(1:NE,1:NKe,1:Nmu,i))
            end do
            !$OMP end parallel do
            do i = 1, j
                tempB = tempB - Hess(i,j)*vec(1:NE,1:NKe,1:Nmu,i)
            end do

            Hess(j+1,j) = norm2(tempB)
            vec(1:NE,1:NKe,1:Nmu,i) = tempB/Hess(j+1,j)

            ! Append latest updated column of Hess to R, with cumulative rotations applied
            Rk(j+1,j) = Hess(j+1,j)
            Rk(1:j,j) = matmul(Qk(1:j,1:j),Hess(1:j,j)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

            ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
            denom = hypot(Rk(j,j),Rk(j+1,j))
            Rot(:,1) = [Rk(j,j)/denom, -Rk(j+1,j)/denom]
            Rot(:,2) = [-Rot(2,1), Rot(1,1)]

            Rk(j:j+1,j) = [denom, 0.0]

            ! Update cumulative rotation matrix
            tempC(1:2,1:j+1) = Qk(j:j+1,1:j+1)
            Qk(j:j+1,1:j+1) = matmul(Rot,tempC(1:2,1:j+1))

            ! Check for convergence
            ompfinish = omp_get_wtime()
            print *, "ITERATE:", j, "CONVERGENCE:", abs(resid*Qk(j+1,1)), &
                "TIME (s):", ompfinish - ompstart
            if (abs(resid*Qk(j+1,1)) <= l_convergence) then
                gk(1:j+1) = resid*Qk(1:j+1,1)
                do i = j, 1, -1
                    y(i) = (gk(i) - sum(Rk(i,i+1:j)*y(i+1:j)))/Rk(i,i)
                end do

                do i = 1, j ! THIS STEP MAY BE EASIER IF RESHAPES AND PACKS ARE USED
                    psi = psi + vec(1:NE,1:NKe,1:Nmu,i)*y(i)
                end do

                exit GMRESm
            end if
            if (isnan(abs(resid*Qk(j+1,1)))) then
                print *, "NAN"
                stop
            end if
        end do

        total = total + 1

        ! Do mth step differently
        tempA = vec(1:NE,1:NKe,1:Nmu,m)

        call wiscoslab_SN_MGXS_swept_Lgg &
            (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

        !$OMP parallel do shared(Hess) private(i)
        do i = 1, m
            Hess(i,m) = sum(tempA*vec(1:NE,1:NKe,1:Nmu,i))
        end do
        !$OMP end parallel do
        !Hess(m+1,m) = sqrt(norm2(tempA)**2 - sum(Hess(1:m,m)**2)) ! Maybe relies on elemental??
        tempB = 0
        do i = 1, m
            tempB = tempB - Hess(i,m)*vec(1:NE,1:NKe,1:Nmu,i)
        end do
        Hess(m+1,m) = norm2(tempA - tempB)

        ! Append latest updated column of Hess to R, with cumulative rotations applied
        Rk(m+1,m) = Hess(m+1,m)
        Rk(1:m,m) = matmul(Qk(1:m,1:m),Hess(1:m,m)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

        ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
        denom = hypot(Rk(m,m),Rk(m+1,m))
        Rot(:,1) = [Rk(m,m)/denom, -Rk(m+1,m)/denom] ! Should I transpose these? were they supposed to be transposed? (~3/31/23)
        Rot(:,2) = [-Rot(2,1), Rot(1,1)]

        Rk(m:m+1,m) = [denom, 0.0]

        ! Update cumulative rotation matrix
        tempC(1:2,1:m+1) = Qk(m:m+1,1:m+1)
        Qk(m:m+1,1:m+1) = matmul(Rot,tempC(1:2,1:m+1))

        ! Construct current solution
        gk(1:m+1) = resid*Qk(1:m+1,1)
        do i = m, 1, -1
            y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
        end do

        do i = 1, m
            psi = psi + vec(1:NE,1:NKe,1:Nmu,i)*y(i)
        end do

        tempA = psi

        call wiscoslab_SN_MGXS_swept_Lgg &
            (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

        resid0 = source - tempA

        !! Slower?? Doesn't work??
        !! Calculate residual for new vector without computing f - Ax0
        !tvec = resid*e1 - matmul(Hess,y)
        !
        !resid0 = tvec(m+1)*tempA/Hess(m+1,m)
        !!!$OMP parallel do shared(resid0) private(i)
        !do i = 1, m
        !    resid0 = resid0 + &
        !    vec(1:NE,1:NKe,1:NQ,i)*(tvec(i)-tvec(m+1)*Hess(i,m)/Hess(m+1,m))
        !end do
        !!!$OMP end parallel do
        print *, "ITERATE:", m, "CONVERGENCE:", abs(gk(m+1)), &
            "TIME (s):", ompfinish - ompstart
        if (abs(gk(m+1)) .le. convergence) then
            exit GMRESm
        end if
        if (isnan(abs(gk(m+1)))) then
            print *, "NAN"
            stop
        end if
        if (total .ge. max_iter) then
            gk(1:m+1) = resid*Qk(1:m+1,1)
            do i = m, 1, -1
                y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
            end do

            do i = 1, m
                psi = psi + vec(1:NE,1:NKe,1:Nmu,i)*y(i)
            end do
            print *, "GMRESm: Did not converge within", max_iter, "iterations."
            print *, "Solution being given with:", abs(gk(j+1)), "residual."
            exit GMRESm
        end if

        print *, "GMRESm: RESTARTING"
    end do GMRESm
end subroutine wiscoslab_SN_MGXS_GMRESm

subroutine wiscoslab_SN_MGXS_swept_Lgg &
    (mu, wL, Fsweep, Sigma, Tprimeinv, tempA, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:), intent(in) :: Tprimeinv

    real, dimension(:,:,:), intent(inout) :: tempA ! Pre-allocated temporary array, saves time. Make sure I don't forget it if I do MMS stuff.
    real, dimension(:,:,:), intent(inout) :: vector

    real, dimension(:,:,:), allocatable :: bdy

    tempA = vector

    call wiscoslab_SN_MGXS_Kgg(mu, wL, Sigma, tempA)

    !!! BIG NOTE: Assumes boundary condition of K*psi term is zero. Not generally true for MMS, but generally true for external beam problems.
    call wiscoslab_SN_MGXS_sweep(mu, Fsweep, bdy, Tprimeinv, tempA)

    vector = vector - tempA
end subroutine wiscoslab_SN_MGXS_swept_Lgg

subroutine wiscoslab_SN_MGXS_Kgg(mu, wL, Sigma, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:), intent(in) :: Sigma ! (l,mat)

    real, dimension(:,:,:), intent(inout) :: vector ! (e,k,i)

    integer :: e, i, ip, k, l
    real, dimension(:,:), allocatable :: P
    real, dimension(:,:,:), allocatable :: tempA ! (e,k,i)

    allocate(P(Nmu,0:NL))
    allocate(tempA(NE,NKe,Nmu))

    P = Legendre_polynomials(NL,mu)
    tempA = vector

    vector = 0
    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(2) shared(vector) private(i, ip, k, l)
        do i = 1, Nmu
            do k = 1, NKe
                do l = 0, NL
                    do ip = 1, Nmu
                        vector(1:NE,k,i) = vector(1:NE,k,i) + &
                            twopi*wL(ip)*P(i,l)*P(ip,l)*tempA(1:NE,k,ip)*Sigma(l+1,1)
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do collapse(2) shared(vector) private(e, i, ip, l)
        do i = 1, Nmu
            do e = 1, NE
                do l = 0, NL
                    do ip = 1, Nmu
                        vector(e,1:NKe,i) = vector(e,1:NKe,i) + &
                            twopi*wL(ip)*P(i,l)*P(ip,l)*tempA(e,1:NKe,ip)*Sigma(l+1,eltomat(e))
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    end if
end subroutine wiscoslab_SN_MGXS_Kgg

subroutine wiscoslab_SN_MGXS_K(mu, wL, Sigma, scatterers, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: Sigma ! (g,l,mat)
    real, dimension(:,:,:,:), intent(in) :: scatterers ! (e,k,i,g)

    real, dimension(:,:,:), intent(inout) :: vector ! (e,k,i)

    integer :: e, i, ip, k, l
    integer :: g
    real, dimension(:,:), allocatable :: P

    g = size(Sigma,1)
    allocate(P(Nmu,0:NL))

    P = Legendre_polynomials(NL,mu)

    vector = 0
    if (Nmats .eq. 1) then ! Would I be able to use this structure if I merged [i,j]? This is much better than how I've been doing it.
        do i = 1, Nmu
            do k = 1, NKe
                do l = 0, NL
                    do ip = 1, Nmu
                        vector(1:NE,k,i) = vector(1:NE,k,i) + &
                            twopi*wL(ip)*P(i,l)*P(ip,l)*matmul(scatterers(1:NE,k,ip,1:g),Sigma(1:g,l+1,1))
                    end do
                end do
            end do
        end do

        if (size(Sigma,2) .eq. NL+2) then
            do i = 1, Nmu
                do k = 1, NKe
                    vector(1:NE,k,i) = vector(1:NE,k,i) + &
                        matmul(scatterers(1:NE,k,i,1:g),Sigma(1:g,NL+2,1))
                end do
            end do
        end if
    else
        do i = 1, Nmu
            do e = 1, NE
                do l = 0, NL
                    do ip = 1, Nmu
                        vector(e,1:NKe,i) = vector(e,1:NKe,i) + &
                            twopi*wL(ip)*P(i,l)*P(ip,l)*&
                            matmul(scatterers(e,1:NKe,ip,1:g),Sigma(1:g,l+1,eltomat(e)))
                    end do
                end do
            end do
        end do

        if (size(Sigma,2) .eq. NL+2) then
            do k = 1, NKe
                do e = 1, NE
                    vector(e,k,1:Nmu) = vector(e,k,1:Nmu) + &
                        matmul(scatterers(e,k,1:Nmu,1:g),Sigma(1:g,NL+2,eltomat(e)))
                end do
            end do
        end if
    end if
end subroutine wiscoslab_SN_MGXS_K

subroutine wiscoslab_SN_MGXS_ETC_uncollided_fluence &
    (rglobal, mu, wL, GpFterm, Fsweep, Sigma, Sigmat, phiu)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:), intent(in) :: GpFterm
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:,:), intent(in) :: Sigma
    real, dimension(:,:), intent(in) :: Sigmat

    real, dimension(:,:,:), intent(inout) :: phiu

    integer :: e, g, i, j, k, l
    integer :: l_G
    real, dimension(:,:), allocatable :: P
    real, dimension(:), allocatable :: singularity
    real, dimension(:,:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:), allocatable :: Tprimeinv
    real, dimension(:,:,:), allocatable :: bdy

    l_G = size(phiu,3)

    allocate(P(Nmu,0:NL))
    allocate(singularity(Nmu))
    allocate(tempA(NE,NKe))
    allocate(t_source(NE,NKe,Nmu))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu))

    P = Legendre_polynomials(NL,mu)
    singularity = 0
    do i = 1, Nmu
        do l = 0, NL
            singularity(i) = singularity(i) + &
                (2*l+1)*P(i,l)*((-1)**l)/fourpi
        end do
    end do

    do g = 2, l_G
        if (Nmats .eq. 1) then
            do i = 1, Nmu
                do k = 1, NKe
                    t_source(1:NE,k,i) = singularity(i)*matmul(phiu(1:NE,k,1:g-1),Sigma(1:g-1,g,1))
                end do
            end do
        else
            do i = 1, Nmu
                do e = 1, NE
                    t_source(e,1:NKe,i) = singularity(i)*matmul(phiu(e,1:NKe,1:g-1),Sigma(1:g-1,g,eltomat(e)))
                end do
            end do
        end if

        call wiscoslab_SN_MGXS_construct_inverted_transport_matrix &
            (mu, GpFterm, Sigmat(g,:), Tprimeinv)
        call wiscoslab_SN_MGXS_sweep(mu, Fsweep, bdy, Tprimeinv, t_source)

        do k = 1, NKe
            tempA(1:NE,k) = twopi*matmul(t_source(1:NE,k,1:Nmu),wL)
        end do

        phiu(1:NE,1:NKe,g) = phiu(1:NE,1:NKe,g) + tempA
    end do
end subroutine wiscoslab_SN_MGXS_ETC_uncollided_fluence

subroutine wiscoslab_MMS(r, dz, mu, wL, Fsweep, Sigma, psi, phi, source, bdy)
    implicit none
    real, dimension(:,:,:), intent(in) :: r
    real, dimension(:), intent(in) :: dz
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigma
    real, dimension(:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:), allocatable, intent(inout) :: phi
    real, dimension(:,:,:,:), allocatable, intent(inout) :: source
    real, dimension(:,:,:), allocatable, intent(inout) :: bdy

    integer :: e, i, ip, k, l
    real, dimension(:,:), allocatable :: P
    real, dimension(:), allocatable :: Sigmai
    real, dimension(:,:,:), allocatable :: ang_fl

    allocate(Sigma(1,1,NL+1,1))
    allocate(psi(NE,NKe,Nmu,1))
    allocate(phi(NE,NKe,1))
    allocate(source(NE,NKe,Nmu,1))
    allocate(P(Nmu,0:NL))
    allocate(Sigmai(Nmu))

    P = Legendre_polynomials(NL,mu)

    do i = 1, Nmu
        Sigmai(i) = full_MMS_scattering_XS(mu(i))
    end do

    Sigma = 0
    do l = 0, NL
        do i = 1, Nmu
            Sigma(:,:,l+1,:) = Sigma(:,:,l+1,:) + &
                (2*l+1)*wL(i)*P(i,l)*Sigmai(i)/2
        end do
    end do

    do i = 1, Nmu
        do k = 1, NKe
            do e = 1, NE
                !psi(e,k,i,1) = exp(mu(i))*r(1,k,e)*(1-r(1,k,e))
                !source(e,k,i,1) = mu(i)*exp(mu(i))*(1-2*r(1,k,e)) + &
                !    MMS_attn*psi(e,k,i,1)
                psi(e,k,i,1) = exp(mu(i))*(r(1,k,e)-0.5)**2
                source(e,k,i,1) = mu(i)*exp(mu(i))*(2*r(1,k,e)-1.0) + &
                    MMS_attn*psi(e,k,i,1)
            end do
        end do
    end do

    do i = 1, Nmu
        do k = 1, NKe
            do e = 1, NE
                do l = 0, NL
                    do ip = 1, Nmu
                        source(e,k,i,1) = source(e,k,i,1) - &
                            twopi*wL(ip)*P(i,l)*P(ip,l)*Sigma(1,1,l+1,1)*psi(e,k,ip,1)
                    end do
                end do
            end do
        end do
    end do

    phi = 0
    do k = 1, NKe
        do e = 1, NE
            do i = 1, Nmu
                phi(e,k,1) = phi(e,k,1) + &
                    twopi*wL(i)*psi(e,k,i,1)
            end do
        end do
    end do

    allocate(bdy(NKe,2,Nmu)) ! (k,bdyels,i)
    allocate(ang_fl(Nmu,2,NKe)) ! (i,bdyels,k)

    ang_fl = 0
    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            ang_fl(i,2,2) = psi(NE,2,i,1)
        else
            ang_fl(i,1,1) = psi(1,1,i,1)
        end if
    end do

    bdy = 0
    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            bdy(1:NKe,2,i) = matmul(Fsweep(1:NKe,1:NKe,1,i,2),ang_fl(i,2,1:NKe))
        else
            bdy(1:NKe,1,i) = matmul(Fsweep(1:NKe,1:NKe,1,i,1),ang_fl(i,1,1:NKe))
        end if
    end do
end subroutine wiscoslab_MMS

subroutine wiscoslab_SN_FEXS_construct_inverted_transport_matrix(mu, GpFterm, Sigmat, Tprimeinv)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:,:,:,:), intent(in) :: GpFterm
    real, dimension(:,:,:), intent(in) :: Sigmat ! (np,n,mat)

    real, dimension(:,:,:,:,:,:), intent(inout) :: Tprimeinv ! (k,kp,np,n,enew,i)

    integer :: dir, e, f, i, j, k, kp, n, np, mat
    real, dimension(:,:,:,:), allocatable :: temp_Tprimeinv
    real, dimension(:,:), allocatable :: ident
    real, dimension(:,:), allocatable :: tempZ
    real, dimension(:,:), allocatable :: tempY
    real, dimension(:,:), allocatable :: inv1
    real, dimension(:,:), allocatable :: inv2
    real, dimension(:,:), allocatable :: inv2inv1
    real, dimension(:,:), allocatable :: inv1inv2
    integer :: estart
    integer :: efin
    integer :: einc
    real, dimension(:,:,:,:,:), allocatable :: tempB

    allocate(temp_Tprimeinv(NKe,NKe,2,2))
    allocate(ident(NKe,NKe))
    allocate(tempZ(NKe,NKe))
    allocate(tempY(NKe,NKe))
    allocate(inv1(NKe,NKe))
    allocate(inv2(NKe,NKe))
    allocate(inv2inv1(NKe,NKe))
    allocate(inv1inv2(NKe,NKe))
    allocate(tempB(NKe,NKe,2,2,NE))

    ident = identity(NKe)

    do i = 1, Nmu
        do e = 1, NE
            mat = eltomat(e)
            inv1 = 0
            inv2 = 0
            tempZ = GpFterm(1:NKe,1:NKe,e,i) + Sigmat(1,1,mat)*ident ! F + G + M11*I, or A
            tempY = GpFterm(1:NKe,1:NKe,e,i) + Sigmat(2,2,mat)*ident ! F + G + M22*I, or D

            inv1 = explicit_2x2_matinv(tempZ) ! (F + G + M11)^(-1)
            tempY = tempY - Sigmat(2,1,mat)*Sigmat(1,2,mat)*inv1
            inv2 = explicit_2x2_matinv(tempY) ! (D-CAinvB)^(-1)
            inv2inv1 = matmul(inv2,inv1)
            inv1inv2 = matmul(inv1,inv2)

            tempZ = tempZ - Sigmat(2,1,mat)*Sigmat(1,2,mat)*inv2
            temp_Tprimeinv(1:NKe,1:NKe,1,1) = inv1 + Sigmat(2,1,mat)*Sigmat(1,2,mat)*matmul(inv1,inv2inv1)
            temp_Tprimeinv(1:NKe,1:NKe,1,2) = -Sigmat(2,1,mat)*inv1inv2
            temp_Tprimeinv(1:NKe,1:NKe,2,1) = -Sigmat(1,2,mat)*inv2inv1
            temp_Tprimeinv(1:NKe,1:NKe,2,2) = inv2

            do k = 1, NKe
                do kp = 1, NKe
                    Tprimeinv(k,kp,:,:,e,i) = transpose(temp_Tprimeinv(k,kp,:,:))
                end do
            end do
        end do
    end do

    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            estart = NE
            efin = 1 ! Make it NE-estart+1
            einc = -1 ! Make it sign(1,int(efin-estart)) or something
        else
            estart = 1
            efin = NE
            einc = 1
        end if

        ! Sweep order stuff
        tempB = Tprimeinv(1:NKe,1:NKe,1:2,1:2,1:NE,i)
        Tprimeinv(1:NKe,1:NKe,1:2,1:2,1:NE,i) = tempB(1:NKe,1:NKe,1:2,1:2,estart:efin:einc)
    end do
end subroutine wiscoslab_SN_FEXS_construct_inverted_transport_matrix

subroutine wiscoslab_SN_FEXS_EI_GMRESm &
    (particle, pmax, mmax, l_G, mu, wL, GpFterm, Fsweep, Sigma, Sigmat, source, bdy, fEg, psi, phi)
    implicit none
    integer, intent(in) :: particle
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    integer, intent(in) :: l_G
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:), intent(in) :: GpFterm
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:,:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:), intent(in) :: Sigmat
    real, dimension(:,:,:,:,:), allocatable, intent(in) :: source
    real, dimension(:,:,:), allocatable, intent(in) :: bdy
    real, dimension(:,:), intent(in) :: fEg

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:,:), allocatable, intent(inout) :: phi

    integer :: e, g, i, j, k, n
    integer :: total
    real, dimension(:,:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:,:), allocatable :: t_Sigma
    real, dimension(:,:,:,:), allocatable :: t_Sigmagg
    real, dimension(:,:,:), allocatable :: t_Sigmat
    real, dimension(:,:,:,:), allocatable :: t_psi
    real, dimension(:,:,:,:,:,:), allocatable :: Tprimeinv
    real, dimension(:,:,:,:), allocatable :: l_bdy
    real :: ompstart, ompfinish

    if (particle .eq. 1) then
        print *, "wiscoslab: SN MGXS GMRESm STARTED, PARTICLE: PHOTON"
    else if (particle .eq. 2) then
        print *, "wiscoslab: SN MGXS GMRESm STARTED, PARTICLE: ELECTRON"
    end if

    allocate(psi(NE,NKe,Nmu,2,l_G))
    allocate(phi(NE,NKe,2,l_G))
    allocate(t_source(NE,NKe,Nmu,2))
    allocate(t_Sigmagg(2,2,size(Sigma,5),size(Sigma,6)))
    allocate(t_Sigmat(2,2,Nmats))
    allocate(t_psi(NE,NKe,Nmu,2))
    allocate(Tprimeinv(NKe,NKe,2,2,NE,Nmu))
    if (allocated(bdy)) allocate(l_bdy(size(bdy,1),size(bdy,2),size(bdy,3),2))

    do g = 1, l_G
        total = 0
        print *, "ENERGY GROUP:", g

        ompstart = omp_get_wtime()
        t_Sigmagg = Sigma(:,:,g,g,:,:)
        t_Sigmat = Sigmat(:,:,g,:)

        call wiscoslab_SN_FEXS_construct_inverted_transport_matrix &
            (mu, GpFterm, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()
        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                t_source = source(:,:,:,:,1)
            else
                t_source = 0
            end if
            if (allocated(bdy)) then
                do n = 1, 2
                    l_bdy(:,:,:,n) = fEg(n,1)*bdy
                end do
            end if
        else
            allocate(t_Sigma, source = Sigma(:,:,1:g-1,g,:,:))
            call wiscoslab_SN_FEXS_K &
                (mu, wL, t_Sigma, psi, t_source)
            if (allocated(source)) then
                t_source = t_source + source(:,:,:,:,g)
            end if
            if (allocated(bdy)) then
                do n = 1, 2
                    l_bdy(:,:,:,n) = fEg(n,g)*bdy
                end do
            end if
            deallocate(t_Sigma)
        end if

        call wiscoslab_SN_FEXS_sweep(mu, Fsweep, l_bdy, Tprimeinv, t_source)

        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        t_psi = t_source ! GUESS
        call wiscoslab_SN_FEXS_GMRESm &
            (convergence, pmax, mmax, mu, wL, &
            Fsweep, t_Sigmagg, Tprimeinv, t_source, l_bdy, total, t_psi)

        psi(:,:,:,:,g) = t_psi

        print *, "CONVERGED WITH:", total, "ITERATES"

        ompfinish = omp_get_wtime()
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    do g = 1, l_G
        do n = 1, 2
            do k = 1, NKe
                phi(1:NE,k,n,g) = twopi*matmul(psi(1:NE,k,1:Nmu,n,g),wL)
            end do
        end do
    end do
end subroutine wiscoslab_SN_FEXS_EI_GMRESm

subroutine wiscoslab_SN_FEXS_GMRESm &
    (l_convergence, pmax, mmax, mu, wL, Fsweep, Sigma, Tprimeinv, source, bdy, total, psi)
    implicit none
    real, intent(in) :: l_convergence
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (np,n,l+1,mat)
    real, dimension(:,:,:,:,:,:), intent(in) :: Tprimeinv ! (k,kp,np,n,ip,i)
    real, dimension(:,:,:,:), intent(in) :: source ! (e,k,i,n)
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdy

    integer, intent(inout) :: total
    real, dimension(:,:,:,:), intent(inout) :: psi ! (e,k,i,n)

    integer :: i, j
    integer :: m
    integer :: max_iter
    real :: resid
    real, dimension(:,:,:,:), allocatable :: tempA
    real, dimension(:,:,:,:), allocatable :: tempB
    real, dimension(:,:,:,:), allocatable :: tempD
    real, dimension(:,:,:,:), allocatable :: resid0
    real, dimension(:,:,:,:,:), allocatable :: vec ! Name it V?
    real, dimension(:,:), allocatable :: Qk
    real, dimension(:,:), allocatable :: local_I
    real, dimension(:,:), allocatable :: Hess
    real, dimension(:,:), allocatable :: Rk
    real, dimension(:,:), allocatable :: Rot
    real, dimension(:,:), allocatable :: tempC
    real, dimension(:), allocatable :: gk
    real, dimension(:), allocatable :: tvec
    real, dimension(:), allocatable :: e1
    real, dimension(:), allocatable :: y
    real :: denom
    real :: ompstart, ompfinish

    m = pmax
    max_iter = mmax

    allocate(tempA(NE,NKe,Nmu,2))
    allocate(tempB(NE,NKe,Nmu,2))
    allocate(tempD(NE,NKe,Nmu,2))
    allocate(resid0(NE,NKe,Nmu,2))
    allocate(vec(NE,NKe,Nmu,2,m))
    allocate(Qk(m+1,m+1))
    allocate(local_I(m+1,m+1))
    allocate(Hess(m+1,m))
    allocate(Rk(m+1,m))
    allocate(Rot(2,2))
    allocate(tempC(2,m+1))
    allocate(gk(m+1))
    allocate(tvec(m+1))
    allocate(e1(m+1))
    allocate(y(m))

    local_I = identity(m+1)
    e1 = 0
    e1(1) = 1

    tempA = psi

    call wiscoslab_SN_FEXS_swept_Lgg &
        (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

    resid0 = source - tempA

    GMRESm: do
        resid = norm2(resid0)

        if (resid .lt. l_convergence) then
            print *, total, resid
            exit GMRESm
        end if

        vec(1:NE,1:NKe,1:Nmu,1:2,1) = resid0/resid

        Qk = local_I

        do j = 1, m-1
            ompstart = omp_get_wtime()
            total = total + 1
            ! Gram-Schmidt orthogonalization to construct vec(:,j) vectors and Hess array
            tempA = vec(1:NE,1:NKe,1:Nmu,1:2,j)

            call wiscoslab_SN_FEXS_swept_Lgg &
                (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

            tempB = tempA
            !$OMP parallel do shared(Hess, tempB) private(i)
            do i = 1, j
                Hess(i,j) = sum(tempA*vec(1:NE,1:NKe,1:Nmu,1:2,i))
            end do
            !$OMP end parallel do
            do i = 1, j
                tempB = tempB - Hess(i,j)*vec(1:NE,1:NKe,1:Nmu,1:2,i)
            end do

            Hess(j+1,j) = norm2(tempB)
            vec(1:NE,1:NKe,1:Nmu,1:2,i) = tempB/Hess(j+1,j)

            ! Append latest updated column of Hess to R, with cumulative rotations applied
            Rk(j+1,j) = Hess(j+1,j)
            Rk(1:j,j) = matmul(Qk(1:j,1:j),Hess(1:j,j)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

            ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
            denom = hypot(Rk(j,j),Rk(j+1,j))
            Rot(:,1) = [Rk(j,j)/denom, -Rk(j+1,j)/denom]
            Rot(:,2) = [-Rot(2,1), Rot(1,1)]

            Rk(j:j+1,j) = [denom, 0.0]

            ! Update cumulative rotation matrix
            tempC(1:2,1:j+1) = Qk(j:j+1,1:j+1)
            Qk(j:j+1,1:j+1) = matmul(Rot,tempC(1:2,1:j+1))

            ! Check for convergence
            ompfinish = omp_get_wtime()
            print *, "ITERATE:", j, "CONVERGENCE:", abs(resid*Qk(j+1,1)), &
                "TIME (s):", ompfinish - ompstart
            if (abs(resid*Qk(j+1,1)) <= l_convergence) then
                gk(1:j+1) = resid*Qk(1:j+1,1)
                do i = j, 1, -1
                    y(i) = (gk(i) - sum(Rk(i,i+1:j)*y(i+1:j)))/Rk(i,i)
                end do

                do i = 1, j ! THIS STEP MAY BE EASIER IF RESHAPES AND PACKS ARE USED
                    psi = psi + vec(1:NE,1:NKe,1:Nmu,1:2,i)*y(i)
                end do

                exit GMRESm
            end if
            if (isnan(abs(resid*Qk(j+1,1)))) then
                print *, "NAN"
                stop
            end if
        end do

        total = total + 1

        ! Do mth step differently
        tempA = vec(1:NE,1:NKe,1:Nmu,1:2,m)

        call wiscoslab_SN_FEXS_swept_Lgg &
            (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

        !$OMP parallel do shared(Hess) private(i)
        do i = 1, m
            Hess(i,m) = sum(tempA*vec(1:NE,1:NKe,1:Nmu,1:2,i))
        end do
        !$OMP end parallel do
        !Hess(m+1,m) = sqrt(norm2(tempA)**2 - sum(Hess(1:m,m)**2)) ! Maybe relies on elemental??
        tempB = 0
        do i = 1, m
            tempB = tempB - Hess(i,m)*vec(1:NE,1:NKe,1:Nmu,1:2,i)
        end do
        Hess(m+1,m) = norm2(tempA - tempB)

        ! Append latest updated column of Hess to R, with cumulative rotations applied
        Rk(m+1,m) = Hess(m+1,m)
        Rk(1:m,m) = matmul(Qk(1:m,1:m),Hess(1:m,m)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

        ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
        denom = hypot(Rk(m,m),Rk(m+1,m))
        Rot(:,1) = [Rk(m,m)/denom, -Rk(m+1,m)/denom] ! Should I transpose these? were they supposed to be transposed? (~3/31/23)
        Rot(:,2) = [-Rot(2,1), Rot(1,1)]

        Rk(m:m+1,m) = [denom, 0.0]

        ! Update cumulative rotation matrix
        tempC(1:2,1:m+1) = Qk(m:m+1,1:m+1)
        Qk(m:m+1,1:m+1) = matmul(Rot,tempC(1:2,1:m+1))

        ! Construct current solution
        gk(1:m+1) = resid*Qk(1:m+1,1)
        do i = m, 1, -1
            y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
        end do

        do i = 1, m
            psi = psi + vec(1:NE,1:NKe,1:Nmu,1:2,i)*y(i)
        end do

        tempA = psi

        call wiscoslab_SN_FEXS_swept_Lgg &
            (mu, wL, Fsweep, Sigma, Tprimeinv, tempD, tempA)

        resid0 = source - tempA

        !! Slower?? Doesn't work??
        !! Calculate residual for new vector without computing f - Ax0
        !tvec = resid*e1 - matmul(Hess,y)
        !
        !resid0 = tvec(m+1)*tempA/Hess(m+1,m)
        !!!$OMP parallel do shared(resid0) private(i)
        !do i = 1, m
        !    resid0 = resid0 + &
        !    vec(1:NE,1:NKe,1:NQ,i)*(tvec(i)-tvec(m+1)*Hess(i,m)/Hess(m+1,m))
        !end do
        !!!$OMP end parallel do
        print *, "ITERATE:", m, "CONVERGENCE:", abs(gk(m+1)), &
            "TIME (s):", ompfinish - ompstart
        if (abs(gk(m+1)) .le. convergence) then
            exit GMRESm
        end if
        if (isnan(abs(gk(m+1)))) then
            print *, "NAN"
            stop
        end if
        if (total .ge. max_iter) then
            gk(1:m+1) = resid*Qk(1:m+1,1)
            do i = m, 1, -1
                y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
            end do

            do i = 1, m
                psi = psi + vec(1:NE,1:NKe,1:Nmu,1:2,i)*y(i)
            end do
            print *, "GMRESm: Did not converge within", max_iter, "iterations."
            print *, "Solution being given with:", abs(gk(j+1)), "residual."
            exit GMRESm
        end if

        print *, "GMRESm: RESTARTING"
    end do GMRESm
end subroutine wiscoslab_SN_FEXS_GMRESm

subroutine wiscoslab_SN_FEXS_swept_Lgg &
    (mu, wL, Fsweep, Sigma, Tprimeinv, tempA, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:,:,:), intent(in) :: Tprimeinv

    real, dimension(:,:,:,:), intent(inout) :: tempA ! Pre-allocated temporary array, saves time. Make sure I don't forget it if I do MMS stuff.
    real, dimension(:,:,:,:), intent(inout) :: vector

    real, dimension(:,:,:,:), allocatable :: bdy

    tempA = vector

    call wiscoslab_SN_FEXS_Kgg(mu, wL, Sigma, tempA)

    !!! BIG NOTE: Assumes boundary condition of K*psi term is zero. Not generally true for MMS, but generally true for external beam problems.
    call wiscoslab_SN_FEXS_sweep(mu, Fsweep, bdy, Tprimeinv, tempA)

    vector = vector - tempA
end subroutine wiscoslab_SN_FEXS_swept_Lgg

subroutine wiscoslab_SN_FEXS_sweep(mu, Fsweep, bdy, Tprimeinv, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdy
    real, dimension(:,:,:,:,:,:), intent(in) :: Tprimeinv

    real, dimension(:,:,:,:), intent(inout) :: vector ! (e,k,i,n)

    integer :: i, ip, j, k, n, np
    integer :: estart
    integer :: efin
    integer :: einc
    integer :: SS1
    real, dimension(:,:,:), allocatable :: l_Fsweep
    real, dimension(:,:,:), allocatable :: l_vector
    real, dimension(:,:), allocatable :: tempA
    logical :: do_bdy

    allocate(l_Fsweep(NKe,NKe,NE))
    allocate(l_vector(NKe,NE,2))
    allocate(tempA(NKe,2))

    do_bdy = allocated(bdy)

    !$OMP parallel do shared(vector) private(i, ip, n, np) &
    !$OMP private(estart, efin, einc, l_Fsweep, l_vector, SS1, tempA)
    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            estart = NE
            efin = 1 ! Make it NE-estart+1
            einc = -1 ! Make it sign(1,int(efin-estart)) or something
            l_Fsweep = Fsweep(1:NKe,1:NKe,1:NE,i,2)
        else
            estart = 1
            efin = NE
            einc = 1
            l_Fsweep = Fsweep(1:NKe,1:NKe,1:NE,i,1)
        end if

        do n = 1, 2
            l_vector(1:NKe,1:NE,n) = transpose(vector(estart:efin:einc,1:NKe,i,n))
        end do

        if (do_bdy) then
            !!! USE THIS STRUCTURE IN 3D. NOT USING HERE BECAUSE SITUATION IS MUCH SIMPLER
            !!! MUST USE THIS BECAUSE SOME BOUNDARY TERMS AREN'T FIRST TERM IN SWEEP.
            !do ip = 1, size(bdy,2)
            !    l_vector(1:NKe,bdyupstr(ip)) = l_vector(1,NKe,bdyupstr(ip)) - bdysrc(1:NKe,ip,i)
            !end do

            l_vector(1:NKe,1,1:2) = l_vector(1:NKe,1,1:2) - &
                merge(bdy(1:NKe,2,i,1:2), bdy(1:NKe,1,i,1:2), mu(i) .gt. 0)
        end if

        SS1 = 1
        do ip = 1, SS1
            tempA = l_vector(1:NKe,ip,1:2)

            l_vector(1:NKe,ip,1:2) = 0
            do n = 1, 2
                do np = 1, 2
                    l_vector(1:NKe,ip,n) = l_vector(1:NKe,ip,n) + &
                        matmul(Tprimeinv(1:NKe,1:NKe,np,n,ip,i),tempA(1:NKe,np))
                end do
            end do
        end do

        do ip = SS1+1, NE
            do n = 1, 2
                tempA(1:NKe,n) = l_vector(1:NKe,ip,n) - matmul(l_Fsweep(1:NKe,1:NKe,ip),l_vector(1:NKe,ip-1,n))
            end do

            l_vector(1:NKe,ip,1:2) = 0
            do n = 1, 2
                do np = 1, 2
                    l_vector(1:NKe,ip,n) = l_vector(1:NKe,ip,n) + &
                        matmul(Tprimeinv(1:NKe,1:NKe,np,n,ip,i),tempA(1:NKe,np))
                end do
            end do
        end do

        do n = 1, 2
            vector(estart:efin:einc,1:NKe,i,n) = transpose(l_vector(1:NKe,1:NE,n))
        end do
    end do
    !$OMP end parallel do
end subroutine wiscoslab_SN_FEXS_sweep

subroutine wiscoslab_SN_FEXS_Kgg(mu, wL, Sigma, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (np,n,l,mat)

    real, dimension(:,:,:,:), intent(inout) :: vector ! (e,k,i)

    integer :: e, i, ip, k, l, n, np
    real, dimension(:,:), allocatable :: P
    real, dimension(:,:,:,:), allocatable :: tempA ! (e,k,i)

    allocate(P(Nmu,0:NL))
    allocate(tempA(NE,NKe,Nmu,2))

    P = Legendre_polynomials(NL,mu)
    tempA = vector

    vector = 0
    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(2) shared(vector) private(i, ip, k, l, n, np)
        do n = 1, 2
            do i = 1, Nmu
                do k = 1, NKe
                    do l = 0, NL
                        do np = 1, 2
                            do ip = 1, Nmu
                                vector(1:NE,k,i,n) = vector(1:NE,k,i,n) + &
                                    twopi*wL(ip)*P(i,l)*P(ip,l)*tempA(1:NE,k,ip,np)*Sigma(np,n,l+1,1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do collapse(2) shared(vector) private(e, i, ip, l, n, np)
        do n = 1, 2
            do i = 1, Nmu
                do e = 1, NE
                    do l = 0, NL
                        do ip = 1, Nmu
                            do np = 1, 2
                                vector(e,1:NKe,i,n) = vector(e,1:NKe,i,n) + &
                                    twopi*wL(ip)*P(i,l)*P(ip,l)*tempA(e,1:NKe,ip,np)*&
                                    Sigma(np,n,l+1,eltomat(e))
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    end if
end subroutine wiscoslab_SN_FEXS_Kgg

subroutine wiscoslab_SN_FEXS_K(mu, wL, Sigma, scatterers, vector)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:,:), intent(in) :: Sigma ! (np,n,g,l,mat)
    real, dimension(:,:,:,:,:), intent(in) :: scatterers ! (e,k,i,n,g)

    real, dimension(:,:,:,:), intent(inout) :: vector ! (e,k,i,n)

    integer :: e, i, ip, k, l, n, np
    integer :: g
    real, dimension(:,:), allocatable :: P

    g = size(Sigma,3)
    allocate(P(Nmu,0:NL))

    P = Legendre_polynomials(NL,mu)

    vector = 0
    if (Nmats .eq. 1) then ! Would I be able to use this structure if I merged [i,j]? This is much better than how I've been doing it.
        do n = 1, 2
            do i = 1, Nmu
                do k = 1, NKe
                    do l = 0, NL
                        do ip = 1, Nmu
                            do np = 1, 2
                                vector(1:NE,k,i,n) = vector(1:NE,k,i,n) + &
                                    twopi*wL(ip)*P(i,l)*P(ip,l)*&
                                    matmul(scatterers(1:NE,k,ip,np,1:g),Sigma(np,n,1:g,l+1,1))
                            end do
                        end do
                    end do
                end do
            end do
        end do

        if (size(Sigma,4) .eq. NL+2) then
            do n = 1, 2
                do i = 1, Nmu
                    do k = 1, NKe
                        do np = 1, 2
                            vector(1:NE,k,i,n) = vector(1:NE,k,i,n) + &
                                matmul(scatterers(1:NE,k,i,np,1:g),Sigma(np,n,1:g,NL+2,1))
                        end do
                    end do
                end do
            end do
        end if
    else
        do n = 1, 2
            do i = 1, Nmu
                do e = 1, NE
                    do l = 0, NL
                        do ip = 1, Nmu
                            do np = 1, 2
                                vector(e,1:NKe,i,n) = vector(e,1:NKe,i,n) + &
                                    twopi*wL(ip)*P(i,l)*P(ip,l)*&
                                    matmul(scatterers(e,1:NKe,ip,np,1:g),&
                                    Sigma(np,n,1:g,l+1,eltomat(e)))
                            end do
                        end do
                    end do
                end do
            end do
        end do

        if (size(Sigma,4) .eq. NL+2) then
            do n = 1, 2
                do k = 1, NKe
                    do e = 1, NE
                        do np = 1, 2
                            vector(e,k,1:Nmu,n) = vector(e,k,1:Nmu,n) + &
                                matmul(scatterers(e,k,1:Nmu,np,1:g),Sigma(np,n,1:g,NL+2,eltomat(e)))
                        end do
                    end do
                end do
            end do
        end if
    end if
end subroutine wiscoslab_SN_FEXS_K

end module t2_ele_iteration