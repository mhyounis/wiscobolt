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

module iteration
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

contains

! SN - general
subroutine SN_read_field &
    (label, l_G, phiu, phi, psi)
    implicit none
    character(*), intent(in) :: label
    integer, intent(in) :: l_G

    real, dimension(:,:,:), allocatable, intent(inout) :: phi
    real, dimension(:,:,:), allocatable, intent(inout) :: phiu
    real, dimension(:,:,:,:), allocatable, intent(inout) :: psi

    character(80) :: file_name

    if (FCS .and. &
        ((label .eq. "electron" .and. &
        transport_mode .eq. "external electron beam") .or. &
        (label .eq. "photon" .and. &
        (transport_mode .eq. "external photon beam" .or. &
        transport_mode .eq. "external photon beam coupled")))) allocate(phiu(NE,NKe,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(psi(NK,Nmu,Nphi,l_G))

    if (FCS .and. &
        ((label .eq. "electron" .and. &
        transport_mode .eq. "external electron beam") .or. &
        (label .eq. "photon" .and. &
        (transport_mode .eq. "external photon beam" .or. &
        transport_mode .eq. "external photon beam coupled")))) then
        file_name = "Stored solutions/"//label//"_uncollided_fluence.dat"

        open(1, file = file_name, form = "unformatted", action = "read")
        read(1) phiu
        close(1)
    end if

    file_name = "Stored solutions/"//label//"_fluence.dat"

    open(1, file = file_name, form = "unformatted", action = "read")
    read(1) phi
    close(1)

    file_name = "Stored solutions/"//label//"_angular_fluence.dat"

    open(1, file = file_name, form = "unformatted", action = "read")
    read(1) psi
    close(1)
end subroutine SN_read_field

subroutine SN_write_field &
    (label, phi, psi)
    implicit none
    character(*), intent(in) :: label
    real, dimension(:,:,:), intent(in) :: phi
    real, dimension(:,:,:,:), intent(in) :: psi

    character(80) :: file_name

    file_name = trim(adjustl(output_fname))//"/"//label//"_fluence.dat"

    open(1, file = file_name, form = "unformatted")
    write(1) phi
    close(1)

    file_name = trim(adjustl(output_fname))//"/"//label//"_angular_fluence.dat"

    open(1, file = file_name, form = "unformatted")
    write(1) psi
    close(1)
end subroutine SN_write_field

! SN - MGXS
!! MOSTLY OPTIMIZED
subroutine SN_MGXS_ETC_uncollided_fluence &
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

    real, dimension(:,:), intent(inout) :: phiu

    integer :: e, g, i, j, k, l
    integer :: l_G
    integer :: l_NK
    real, dimension(:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: singularity
    real, dimension(:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:,:), allocatable :: Tprimeinv
    real, dimension(:,:,:,:), allocatable :: nobdysrc

    l_G = size(phiu,2)
    l_NK = size(nodesinbeam)

    allocate(tempA(NK))
    allocate(singularity(NK,Nmu,Nphi))
    allocate(t_source(NK,Nmu,Nphi))
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
            tempA = matmul(phiu(1:NK,1:g-1),Sigma(1:g-1,g,1))
            do j = 1, Nphi
                do i = 1, Nmu
                    t_source(1:NK,i,j) = tempA*singularity(1:NK,i,j)
                end do
            end do
        else
            ! EXTREMELY UNOPTIMIZED
            tempA = 0
            do k = 1, NKe
                do e = 1, NE
                    tempA(Cekk(e,k)) = tempA(Cekk(e,k)) + &
                        dot_product(phiu(Cekk(e,k),1:g-1),Sigma(1:g-1,g,eltomat(e)))
                end do
            end do
            do j = 1, Nphi
                do i = 1, Nmu
                    t_source(1:NK,i,j) = tempA*singularity(1:NK,i,j)
                end do
            end do
        end if

        call SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat(g,:), Tprimeinv)
        call SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
            varsigmadown, I1invI2f, nobdysrc, Tprimeinv, t_source)

        tempA = 0
        do j = 1, Nphi
            tempA = tempA + twopi*matmul(t_source(1:NK,1:Nmu,j),wL)/Nphi
        end do

        phiu(1:NK,g) = phiu(1:NK,g) + tempA
    end do
end subroutine SN_MGXS_ETC_uncollided_fluence

subroutine SN_MGXS_construct_inverted_transport_matrix &
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
end subroutine SN_MGXS_construct_inverted_transport_matrix

subroutine SN_MGXS_EI_SI &
    (particle, pmax, l_G, Cekk, bdyel, eprime, normal, wL, khat, cosmat, qPa, qfactorialmat, &
    esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
    I1invI2f, I1invI3vector, Sigma, Sigmat, phiu, source, fEg, bdysrc, psi, phi)
    implicit none
    integer, intent(in) :: particle
    integer, intent(in) :: pmax
    integer, intent(in) :: l_G
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: bdyel
    integer, dimension(:,:), intent(in) :: eprime
    real, dimension(:,:,:), intent(in) :: normal
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
    real, dimension(:,:,:), allocatable, intent(in) :: phiu
    real, dimension(:,:,:,:), allocatable, intent(in) :: source
    real, dimension(:), intent(in) :: fEg
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc

    real, dimension(:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:), allocatable, intent(inout) :: phi

    integer :: e, g, k
    integer :: total
    real, dimension(:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:), allocatable :: t_bdysrc
    real, dimension(:,:,:), allocatable :: t_Sigma
    real, dimension(:,:), allocatable :: t_Sigmagg
    real, dimension(:), allocatable :: t_Sigmat
    real, dimension(:), allocatable :: t_phi
    real, dimension(:,:,:), allocatable :: t_psi
    real, dimension(:,:,:,:,:), allocatable :: Tprimeinv
    real, dimension(:,:,:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: tempB
    real :: ompstart, ompfinish

    if (particle .eq. 1) then
        print *, "SN MGXS SI STARTED, PARTICLE: PHOTON"
    else if (particle .eq. 2) then
        print *, "SN MGXS SI STARTED, PARTICLE: ELECTRON"
    end if

    allocate(psi(NK,Nmu,Nphi,l_G))
    allocate(phi(NE,NKe,l_G))

    allocate(t_source(NK,Nmu,Nphi))
    if (allocated(bdysrc)) allocate(t_bdysrc(NKe,size(bdysrc,2),Nmu,Nphi))
    allocate(t_Sigmagg(size(Sigma,3),size(Sigma,4)))
    allocate(t_Sigmat(Nmats))
    allocate(t_phi(NK))
    allocate(t_psi(NK,Nmu,Nphi))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

    if (check_residuals) then
        allocate(tempA(NK,Nmu,Nphi))
        allocate(tempB(NK,Nmu,Nphi))
    end if

    do g = 1, l_G
        total = 0
        print *, "ENERGY GROUP:", g

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                t_source = source(:,:,:,1)
            else
                t_source = 0
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(1)*bdysrc
            end if
        else
            allocate(t_Sigma, source = Sigma(1:g-1,g,:,:))
            call SN_MGXS_K &
                (wL, cosmat, qPa, qfactorialmat, t_Sigma, psi, t_source)
            if (allocated(source)) then
                t_source = t_source + source(:,:,:,g)
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(g)*bdysrc
            end if
            deallocate(t_Sigma)
        end if
        if (check_residuals) tempB = t_source

        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        t_Sigmagg = Sigma(g,g,:,:)
        t_Sigmat = Sigmat(g,:)
        if (allocated(phiu)) then
            t_phi = 0
            do k = 1, NKe
                do e = 1, NE
                    t_phi(Cekk(e,k)) = t_phi(Cekk(e,k)) + &
                        phiu(e,k,g)/Nglobal(Cekk(e,k))
                end do
            end do
        else
            t_phi = 0
        end if

        ompstart = omp_get_wtime()
        call SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()

        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        call SN_MGXS_SI &
            (convergence, pmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, &
            esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, t_Sigmagg, t_source, t_bdysrc, Tprimeinv, total, t_phi, t_psi)

        psi(:,:,:,g) = t_psi
        if (allocated(phiu)) then
            do k = 1, NKe
                do e = 1, NE
                    phi(e,k,g) = t_phi(Cekk(e,k)) - phiu(e,k,g)
                end do
            end do
        else
            do k = 1, NKe
                do e = 1, NE
                    phi(e,k,g) = t_phi(Cekk(e,k))
                end do
            end do
        end if

        print *, "CONVERGED WITH:", total, "ITERATES"

        if (check_residuals) then
            print *, "CHECKING RESIDUAL:"
            print *, "CHECKING RESIDUAL IS WIP"
            exit
            tempA = t_psi
            !call SN_MGXS_T &
            !    (Cekk, eprime, normal, khat, varsigmaup, I1invI2, I1invI2f, &
            !    I1invI3vector, t_Sigmat, tempA) ! tempA = T*psi
            call SN_MGXS_Kgg &
                (wL, cosmat, qPa, qfactorialmat, t_Sigmagg, t_psi) ! t_psi = K*psi

            tempA = tempA - t_psi
            tempA = tempB - tempA
            print *, norm2(tempA)
        end if

        ompfinish = omp_get_wtime()
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    if (particle .eq. 1) then
        call SN_write_field("photon", phi, psi)
    else if (particle .eq. 2) then
        call SN_write_field("electron", phi, psi)
    end if
end subroutine SN_MGXS_EI_SI

subroutine SN_MGXS_SI &
    (l_convergence, pmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, &
    enew, fold, Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, source, bdysrc, Tprimeinv, total, phi, psi)
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
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc
    real, dimension(:,:,:,:,:), intent(in) :: Tprimeinv ! (kp,k,enew,i,j)

    integer, intent(inout) :: total
    real, dimension(:,:,:), intent(inout) :: source ! (k,i,j)
    real, dimension(:), intent(inout) :: phi ! (k) ! phiu initially, leaves as phic + phiu
    real, dimension(:,:,:), intent(inout) :: psi ! (k,i,j)

    integer :: j, p
    real, dimension(:,:), allocatable :: jterm
    real, dimension(:), allocatable :: tempE
    real :: prevnorm
    real :: nextnorm
    real :: ompstart, ompfinish

    allocate(jterm(NK,Nphi))
    allocate(tempE(NK))

    ompstart = omp_get_wtime()

    prevnorm = norm2(phi)

    call SN_MGXS_sweep &
        (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
        I1invI2f, bdysrc, Tprimeinv, source)

    psi = source
    !$OMP parallel do shared(jterm) private(j)
    do j = 1, Nphi
        jterm(1:NK,j) = twopi*matmul(source(1:NK,1:Nmu,j),wL)/Nphi
    end do
    !$OMP end parallel do
    tempE = sum(jterm,2)
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

        call SN_MGXS_Kgg &
            (wL, cosmat, qPa, qfactorialmat, Sigma, source)

        call SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, bdysrc, Tprimeinv, source)

        psi = psi + source
        !$OMP parallel do shared(jterm) private(j)
        do j = 1, Nphi
            jterm(1:NK,j) = twopi*matmul(source(1:NK,1:Nmu,j),wL)/Nphi
        end do
        !$OMP end parallel do
        tempE = sum(jterm,2)
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
end subroutine SN_MGXS_SI

subroutine SN_MGXS_EI_GMRESm &
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
    real, dimension(:,:,:,:), allocatable, intent(in) :: source
    real, dimension(:), intent(in) :: fEg
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc

    real, dimension(:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:), allocatable, intent(inout) :: phi

    integer :: g, i, j, k
    integer :: total
    real, dimension(:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:), allocatable :: t_bdysrc
    real, dimension(:,:,:), allocatable :: t_Sigma
    real, dimension(:,:), allocatable :: t_Sigmagg
    real, dimension(:), allocatable :: t_Sigmat
    real, dimension(:,:,:), allocatable :: t_psi
    real, dimension(:,:,:,:,:), allocatable :: Tprimeinv
    real :: ompstart, ompfinish

    if (particle .eq. 1) then
        print *, "SN MGXS GMRESm STARTED, PARTICLE: PHOTON"
    else if (particle .eq. 2) then
        print *, "SN MGXS GMRESm STARTED, PARTICLE: ELECTRON"
    else if (particle .eq. 3) then
        print *, "SN MGXS GMRESm STARTED, PARTICLE: MMS"
    end if

    allocate(psi(NK,Nmu,Nphi,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(t_source(NK,Nmu,Nphi))
    if (allocated(bdysrc)) allocate(t_bdysrc(NKe,size(bdysrc,2),Nmu,Nphi))
    allocate(t_Sigmagg(size(Sigma,3),size(Sigma,4)))
    allocate(t_Sigmat(Nmats))
    allocate(t_psi(NK,Nmu,Nphi))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

    do g = 1, l_G
        total = 0
        print *, "ENERGY GROUP:", g

        ompstart = omp_get_wtime()
        t_Sigmagg = Sigma(g,g,:,:)
        t_Sigmat = Sigmat(g,:)

        call SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()
        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                t_source = source(:,:,:,1)
            else
                t_source = 0
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(1)*bdysrc
            end if
        else
            allocate(t_Sigma, source = Sigma(1:g-1,g,:,:))
            call SN_MGXS_K &
                (wL, cosmat, qPa, qfactorialmat, t_Sigma, psi, t_source)
            if (allocated(source)) then
                t_source = t_source + source(:,:,:,g)
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(g)*bdysrc
            end if
            deallocate(t_Sigma)
        end if

        call SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
            varsigmadown, I1invI2f, t_bdysrc, Tprimeinv, t_source)

        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        !t_psi = t_source/t_Sigmat(:,1) ! GUESS
        t_psi = t_source
        call SN_MGXS_GMRESm &
            (convergence, pmax, mmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, &
            esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, t_Sigmagg, t_source, t_bdysrc, Tprimeinv, total, t_psi)

        psi(:,:,:,g) = t_psi

        print *, "CONVERGED WITH:", total, "ITERATES"

        ompfinish = omp_get_wtime()
        print *, "TOTAL TIME (s):", ompfinish - ompstart

        !allocate(t_Sigma(NK,Nmu,Nphi))
        !call SN_MGXS_swept_Lgg &
        !    (Cekk, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
        !    varsigmadown, I1invI2f, t_Sigmagg, Tprimeinv, t_Sigma, t_psi)
        !print *, norm2(t_source - t_psi)
        !stop
    end do

    phi = 0
    do g = 1, l_G
        do j = 1, Nphi
            do k = 1, NKe
                phi(1:NE,k,g) = phi(1:NE,k,g) + twopi*matmul(psi(Cekk(1:NE,k),1:Nmu,j,g),wL)/Nphi
            end do
        end do
    end do

    if (particle .eq. 1) then
        call SN_write_field("photon", phi, psi)
    else if (particle .eq. 2) then
        call SN_write_field("electron", phi, psi)
    end if
end subroutine SN_MGXS_EI_GMRESm

subroutine SN_MGXS_GMRESm &
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
    real, dimension(:,:,:), intent(in) :: source ! (k,i,j)
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc
    real, dimension(:,:,:,:,:), intent(in) :: Tprimeinv ! (kp,k,enew,i,j)

    integer, intent(inout) :: total
    real, dimension(:,:,:), intent(inout) :: psi ! (k,i,j)

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

    allocate(tempA(NK,Nmu,Nphi))
    allocate(tempB(NK,Nmu,Nphi))
    allocate(tempD(NK,Nmu,Nphi))
    allocate(resid0(NK,Nmu,Nphi))
    allocate(vec(NK,Nmu,Nphi,m))
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

    call SN_MGXS_swept_Lgg &
        (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
        Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempD, tempA)

    resid0 = source - tempA

    GMRESm: do
        resid = norm2(resid0)

        if (resid .lt. l_convergence) then
            print *, "ITERATE:", 1, "CONVERGENCE:", abs(gk(m+1)), &
                "TIME (s):", ompfinish - ompstart
            exit GMRESm
        end if

        vec(1:NK,1:Nmu,1:Nphi,1) = resid0/resid

        Qk = local_I

        do j = 1, m-1
            ompstart = omp_get_wtime()
            total = total + 1
            ! Gram-Schmidt orthogonalization to construct vec(:,j) vectors and Hess array
            tempA = vec(1:NK,1:Nmu,1:Nphi,j)

            call SN_MGXS_swept_Lgg &
                (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempD, tempA)

            tempB = tempA
            !$OMP parallel do shared(Hess, tempB) private(i)
            do i = 1, j
                Hess(i,j) = sum(tempA*vec(1:NK,1:Nmu,1:Nphi,i))
            end do
            !$OMP end parallel do
            do i = 1, j
                tempB = tempB - Hess(i,j)*vec(1:NK,1:Nmu,1:Nphi,i)
            end do

            Hess(j+1,j) = norm2(tempB)
            vec(1:NK,1:Nmu,1:Nphi,i) = tempB/Hess(j+1,j)

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
                    psi = psi + vec(1:NK,1:Nmu,1:Nphi,i)*y(i)
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
        tempA = vec(1:NK,1:Nmu,1:Nphi,m)

        call SN_MGXS_swept_Lgg &
            (Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
            Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, Tprimeinv, tempD, tempA)

        !$OMP parallel do shared(Hess) private(i)
        do i = 1, m
            Hess(i,m) = sum(tempA*vec(1:NK,1:Nmu,1:Nphi,i))
        end do
        !$OMP end parallel do
        !Hess(m+1,m) = sqrt(norm2(tempA)**2 - sum(Hess(1:m,m)**2)) ! Maybe relies on elemental??
        tempB = 0
        do i = 1, m
            tempB = tempB - Hess(i,m)*vec(1:NK,1:Nmu,1:Nphi,i)
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
            psi = psi + vec(1:NK,1:Nmu,1:Nphi,i)*y(i)
        end do

        tempA = psi

        call SN_MGXS_swept_Lgg &
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
                psi = psi + vec(1:NK,1:Nmu,1:Nphi,i)*y(i)
            end do
            print *, "GMRESm: Did not converge within", max_iter, "iterations."
            print *, "Solution being given with:", abs(gk(j+1)), "residual."
            exit GMRESm
        end if

        print *, "GMRESm: RESTARTING"
    end do GMRESm
end subroutine SN_MGXS_GMRESm

subroutine SN_MGXS_K &
    (wL, cosmat, qPa, qfactorialmat, Sigma, scatterers, vector)
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat ! (j,jp,m+1)
    real, dimension(:,:), intent(in) :: qPa ! (i,q)
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:), intent(in) :: Sigma ! (gp,l+1)
    real, dimension(:,:,:,:), intent(in) :: scatterers ! (k,i,j,g)

    real, dimension(:,:,:), intent(inout) :: vector ! (k,i,j)

    integer :: e, i, ip, j, jp, k, l, mat, q, g
    integer :: NQ
    integer :: group
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:), allocatable :: Amat ! (k,ip,jp,l)
    real, dimension(:,:,:), allocatable :: Bmat ! (k,jp,q)
    real, dimension(:,:,:), allocatable :: Cmat ! (k,jp,q)
    !real, dimension(Nphi) :: tempA

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    allocate(Amat(NK,Nmu,Nphi,NL+1))
    allocate(Bmat(NK,Nphi,NQ))
    allocate(Cmat(NK,Nphi,NQ))

    group = size(Sigma,1)
    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(3) shared(Amat) private(l, jp, ip)
        do l = 0, NL
            do jp = 1, Nphi
                do ip = 1, Nmu
                    Amat(1:NK,ip,jp,l+1) = twopi*matmul(scatterers(1:NK,ip,jp,1:group),Sigma(1:group,l+1,1))/Nphi
                end do
            end do
        end do
        !$OMP end parallel do
    else
        Amat = 0
        do mat = 1, Nmats
            !$OMP parallel do collapse(4) shared(Amat) private(ip, jp, k, l)
            do l = 0, NL
                do jp = 1, Nphi
                    do ip = 1, Nmu
                        do k = 1, NKmat(mat)
                            Amat(nodesinmat(k,mat),ip,jp,l+1) = Amat(nodesinmat(k,mat),ip,jp,l+1) + &
                                twopi*Nelnodes(nodesinmat(k,mat),mat)*&
                                sum(scatterers(nodesinmat(k,mat),ip,jp,1:group)*Sigma(1:group,l+1,mat))&
                                /(Nphi*Nglobal(nodesinmat(k,mat)))
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        end do
    end if

    !$OMP parallel do shared(Bmat) private(q, jp)
    do q = 1, NQ
        do jp = 1, Nphi
            Bmat(1:NK,jp,q) = matmul(Amat(1:NK,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
        end do
        Bmat(1:NK,1:Nphi,q) = qfactorialmat(q)*Bmat(1:NK,1:Nphi,q)
    end do
    !$OMP end parallel do

    !$OMP parallel do shared(Cmat) private(q)
    do q = 1, NQ
        Cmat(1:NK,1:Nphi,q) = matmul(Bmat(1:NK,1:Nphi,q),cosmat(1:Nphi,1:Nphi,emm(q)+1))
    end do
    !$OMP end parallel do

    if (size(Sigma,2) .eq. NL+2) then
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(2) shared(vector) private(i, j)
            do j = 1, Nphi
                do i = 1, Nmu
                    vector(1:NK,i,j) = matmul(scatterers(1:NK,i,j,1:group),Sigma(1:group,NL+2,1))
                end do
            end do
            !$OMP end parallel do
        else
            vector = 0
            do mat = 1, Nmats
                !$OMP parallel do collapse(3) shared(vector) private(i, j, k)
                do j = 1, Nphi
                    do i = 1, Nmu
                        do k = 1, NKmat(mat)
                            vector(nodesinmat(k,mat),i,j) = vector(nodesinmat(k,mat),i,j) + &
                                Nelnodes(nodesinmat(k,mat),mat)*&
                                sum(scatterers(nodesinmat(k,mat),i,j,1:group)*Sigma(1:group,NL+2,mat))/&
                                Nglobal(nodesinmat(k,mat))
                        end do
                    end do
                end do
                !$OMP end parallel do
            end do
        end if
    else
       vector = 0
    end if

    !$OMP parallel do shared(vector) private(i, q)
    do i = 1, Nmu
        do q = 1, NQ
            vector(1:NK,i,1:Nphi) = vector(1:NK,i,1:Nphi) + Cmat(1:NK,1:Nphi,q)*qPa(i,q)
        end do
    end do
    !$OMP end parallel do
end subroutine SN_MGXS_K

subroutine SN_MGXS_Kgg &
    (wL, cosmat, qPa, qfactorialmat, Sigma, vector)
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat ! (j,jp,m+1)
    real, dimension(:,:), intent(in) :: qPa ! (i,q)
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:), intent(in) :: Sigma ! (l,mat)

    real, dimension(:,:,:), intent(inout) :: vector ! (k,i,j)

    integer :: e, i, ip, j, jp, k, l, mat, q
    integer :: NQ
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:), allocatable :: Amat ! (j,i,k,l)
    real, dimension(:,:,:), allocatable :: Bmat ! (j,k,q)
    real, dimension(:,:,:), allocatable :: Cmat ! (j,k,q)
    real, dimension(:,:,:), allocatable :: gold
    real, dimension(Nphi) :: tempA

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    allocate(Amat(NK,Nmu,Nphi,NL+1))
    allocate(Bmat(NK,Nphi,NQ))
    allocate(Cmat(NK,Nphi,NQ))

    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    if (Nmats .eq. 1) then
        !$OMP parallel do shared(Amat) private(l)
        do l = 0, NL
            Amat(1:NK,1:Nmu,1:Nphi,l+1) = twopi*Sigma(l+1,1)*vector/Nphi
        end do
        !$OMP end parallel do
    else
        Amat = 0
        do mat = 1, Nmats
            !$OMP parallel do collapse(2) shared(Amat) private(k, l)
            do l = 0, NL
                do k = 1, NKmat(mat)
                    Amat(nodesinmat(k,mat),1:Nmu,1:Nphi,l+1) = Amat(nodesinmat(k,mat),1:Nmu,1:Nphi,l+1) + &
                        twopi*Nelnodes(nodesinmat(k,mat),mat)*&
                        Sigma(l+1,mat)*vector(nodesinmat(k,mat),1:Nmu,1:Nphi)/(Nphi*Nglobal(nodesinmat(k,mat)))
                end do
            end do
            !$OMP end parallel do
        end do
    end if

    !$OMP parallel do shared(Bmat) private(q, jp)
    do q = 1, NQ
        do jp = 1, Nphi
            Bmat(1:NK,jp,q) = matmul(Amat(1:NK,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
        end do
        Bmat(1:NK,1:Nphi,q) = qfactorialmat(q)*Bmat(1:NK,1:Nphi,q)
    end do
    !$OMP end parallel do

    !$OMP parallel do shared(Cmat) private(q)
    do q = 1, NQ
        Cmat(1:NK,1:Nphi,q) = matmul(Bmat(1:NK,1:Nphi,q),cosmat(1:Nphi,1:Nphi,emm(q)+1))
    end do
    !$OMP end parallel do

    vector = 0
    !$OMP parallel do shared(vector) private(i, q) ! collapse(1) this?
    do i = 1, Nmu
        do q = 1, NQ
            vector(1:NK,i,1:Nphi) = vector(1:NK,i,1:Nphi) + Cmat(1:NK,1:Nphi,q)*qPa(i,q)
        end do
    end do
    !$OMP end parallel do
end subroutine SN_MGXS_Kgg

subroutine SN_MGXS_sweep &
    (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
    varsigmadown, I1invI2f, bdysrc, Tprimeinv, vector)
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

    real, dimension(:,:,:), intent(inout) :: vector

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
                    l_vector(k,l_enew(e)) = vector(Cekk(e,k),i,j)
                end do
            end do

            if (use_bdy) then
                do e = 1, NBE
                    l_vector(1:NKe,l_enew(bdyel(e))) = l_vector(1:NKe,l_enew(bdyel(e))) + &
                        bdysrc(1:NKe,e,i,j)
                end do
            end if

            l_fold = fold(1:NFe,1:NE,i,j)
            l_Ia = Ia(1:NE+1,i,j) ! why isn't the NE+1 term ever used? i forgot. there was a good reason though.
            l_Ja = Ja(1:noJa,i,j)
            l_varsigmadown = varsigmadown(1:noJa,i,j)

            SS1 = esweepbounds(2,i,j)
            do ip = 1, SS1
                tempA = l_vector(1:NKe,ip)
                do k = 1, NKe
                    l_vector(k,ip) = sum(Tprimeinv(1:NKe,k,ip,i,j)*tempA)
                end do
            end do

            !!! MUST ENSURE THAT I AM USING I2f'S SELF-SOURCE TERM. I NEGLECTED IT IN SN_MGXS_T... ! I THINK I AM NOT!!! HUGE... DOESN'T AFFECT PREVIOUS MMS THOUGH...
            !!! Check T*T^(-1) to prove?? Never tried this with the sphere, so may be huge change...
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

            vector(1:NK,i,j) = 0
            do k = 1, NKe
                do e = 1, NE
                    vector(Cekk(l_esweeplist(e),k),i,j) = vector(Cekk(l_esweeplist(e),k),i,j) + l_vector(k,e) ! MAKE THIS FASTER!!
                    ! Apparently doing addition like this is slow.
                end do
            end do

            vector(1:NK,i,j) = vector(1:NK,i,j)/Nglobal
        end do
    end do
    !$OMP end parallel do
end subroutine SN_MGXS_sweep

subroutine SN_MGXS_swept_Lgg &
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

    real, dimension(:,:,:), intent(inout) :: tempA ! Pre-allocated temporary array, saves time.
    real, dimension(:,:,:), intent(inout) :: vector

    real, dimension(:,:,:,:), allocatable :: nobdysrc

    tempA = vector

    call SN_MGXS_Kgg(wL, cosmat, qPa, qfactorialmat, Sigma, tempA)

    call SN_MGXS_sweep &
        (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
        varsigmadown, I1invI2f, nobdysrc, Tprimeinv, tempA)

    vector = vector - tempA
end subroutine SN_MGXS_swept_Lgg

subroutine SN_MGXS_T &
    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, Sigmat, &
    bdysrc, vector)
    ! NOT OPTIMIZED, NOT RECOMMENDED, NOT USED FREQUENTLY DUE TO SUPERIORITY OF SWEPT L IN SN_MGXS_GMRES.
    ! This is primarily for MMS purposes, to verify Tf = s for a given f and s.
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

    real, dimension(:,:,:), intent(inout) :: vector

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

    !$OMP parallel do collapse(3) shared(Tprime) private(dir, e, f, i, j) private(tempD)
    do j = 1, Nphi
        do i = 1, Nmu
            do e = 1, NE
                tempD = 0
                ! FTERM
                do f = 1, NFe
                    tempD = tempD + &
                    varsigmaup(f,e,i,j)*I1invI2(1:NKe,1:NKe,f,e)
                end do
                ! GTERM
                do dir = 1, 3
                    tempD = tempD - &
                    I1invI3vector(1:NKe,1:NKe,dir,e)*khat(dir,i,j)
                end do
                ! ADD TO MTERM
                Tprime(1:NKe,1:NKe,e,i,j) = tempD + Sigmat(eltomat(e))*ident
            end do
        end do
    end do
    !$OMP end parallel do

    Mterm = 0
    if (use_bdy) then
        do ip = 1, size(bdyel)
            Mterm(bdyel(ip),1:NKe,1:Nmu,1:Nphi) = Mterm(bdyel(ip),1:NKe,1:Nmu,1:Nphi) - & ! Uses minus because bdysrc is on RHS generally
                bdysrc(1:NKe,ip,1:Nmu,1:Nphi)
        end do
    end if

    !$OMP parallel do collapse(4) shared(Mterm) private(j, i, e, k, kp)
    do j = 1, Nphi
        do i = 1, Nmu
            do e = 1, NE
                do k = 1, NKe
                    do kp = 1, NKe
                        Mterm(e,k,i,j) = Mterm(e,k,i,j) + Tprime(k,kp,e,i,j)*vector(Cekk(e,kp),i,j)
                    end do
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    ! SWEPT FTERM
    !!$OMP parallel do collapse(2) shared(Mterm) private(f, e, j, i, k, kp) &
    !!$OMP private(l_varsigmadown)
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
                                l_varsigmadown*I1invI2f(kp,k,f,e)*vector(Cekk(l_e,kp),i,j)
                        end do
                    end do
                end do
            end do
        end do
    end do
    !!$OMP end parallel do

    vector = 0
    do k = 1, NKe
        do e = 1, NE
            vector(Cekk(e,k),1:Nmu,1:Nphi) = vector(Cekk(e,k),1:Nmu,1:Nphi) + &
                Mterm(e,k,1:Nmu,1:Nphi)/Nglobal(Cekk(e,k))
        end do
    end do
end subroutine SN_MGXS_T

subroutine SN_MGXS_residual &
    (Cekk, eprime, normal, wL, khat, varsigmaup, cosmat, qPa, qfactorialmat, &
    I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, source, psi, residual)
    ! CURRENTLY SLOW AND REDUNDANT
    ! OPTIMIZE, ALSO HAVE IT CREATE THE SOURCE FROM SCRATCH SO THIS DOESN'T RUN DURING SI.
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:,:), intent(in) :: eprime
    real, dimension(:,:,:), intent(in) :: normal
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:), intent(in) :: Sigmat
    real, dimension(:,:,:,:), intent(in) :: source
    real, dimension(:,:,:,:), intent(in) :: psi

    real, dimension(:), allocatable, intent(inout) :: residual
    !real, dimension(:,:,:), intent(inout) :: residual

    integer :: g
    integer :: l_G
    real, dimension(:,:,:), allocatable :: tempA
    real, dimension(:,:,:,:), allocatable :: tempAA
    real, dimension(:,:,:), allocatable :: temp_Sigma

    l_G = size(psi,4)
    allocate(tempA(NK,Nmu,Nphi))
    allocate(tempAA(NK,Nmu,Nphi,l_G))
    allocate(temp_Sigma(l_G,NL+1,Nmats))
    allocate(residual(l_G))

    do g = 1, l_G
        tempA = psi(:,:,:,g)
        !call SN_MGXS_T &
        !(Cekk, eprime, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, Sigmat, tempA)
        tempAA(:,:,:,g) = tempA
    end do

    deallocate(tempA)

    do g = 1, l_G
        temp_Sigma(1:g,1:NL+1,:) = Sigma(1:g,g,1:NL+1,:)
        print *, "I CHANGED THIS SUBROUTINE"
        stop
        !call SN_MGXS_K(g, wL, cosmat, qPa, qfactorialmat, temp_Sigma(1:g,1:NL+1,:), psi(:,:,:,1:g), tempA)
        tempAA(:,:,:,g) = tempAA(:,:,:,g) - tempA
        deallocate(tempA)
    end do

    !residual = (1 - (tempA - tempB)/source)
    do g = 1, l_G
        residual(g) = norm2(source(:,:,:,g) - tempAA(:,:,:,g))
    end do
end subroutine SN_MGXS_residual

! SN - FEXS
!! MOSTLY UNOPTIMIZED
subroutine SN_FEXS_EI_SI &
    (particle, pmax, l_G, Cekk, bdyel, wL, khat, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, &
    enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, Sigma, &
    Sigmat, source, fEg, bdysrc, psi, phi)
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
    real, dimension(:,:,:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:), intent(in) :: Sigmat
    real, dimension(:,:,:,:,:), allocatable, intent(in) :: source ! (k,i,j,n,g)
    real, dimension(:,:), intent(in) :: fEg
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc

    real, dimension(:,:,:,:), allocatable, intent(inout) :: psi
    real, dimension(:,:,:), allocatable, intent(inout) :: phi

    integer :: e, g, k, n
    integer :: total
    real, dimension(:,:,:,:), allocatable :: t_source
    real, dimension(:,:,:,:,:), allocatable :: t_bdysrc
    real, dimension(:,:,:,:,:), allocatable :: t_Sigma
    real, dimension(:,:,:,:), allocatable :: t_Sigmagg
    real, dimension(:,:,:), allocatable :: t_Sigmat
    real, dimension(:,:), allocatable :: t_phi
    real, dimension(:,:,:,:), allocatable :: t_psi
    real, dimension(:,:,:,:,:,:,:), allocatable :: Tprimeinv
    real :: ompstart, ompfinish

    if (particle .eq. 1) then
        print *, "SN FEXS SI STARTED, PARTICLE: PHOTON"
    else if (particle .eq. 2) then
        print *, "SN FEXS SI STARTED, PARTICLE: ELECTRON"
    end if

    allocate(psi(NK,Nmu,Nphi,l_G+1))
    allocate(phi(NE,NKe,l_G+1))
    allocate(t_source(NK,Nmu,Nphi,2))
    if (allocated(bdysrc)) allocate(t_bdysrc(NKe,size(bdysrc,2),Nmu,Nphi,2))
    allocate(t_Sigmagg(2,2,size(Sigma,3),size(Sigma,4)))
    allocate(t_Sigmat(2,2,Nmats))
    allocate(t_phi(NK,2))
    allocate(t_psi(NK,Nmu,Nphi,2))
    allocate(Tprimeinv(NKe,NKe,2,2,NE,Nmu,Nphi))

    do g = 1, l_G
        total = 0
        print *, "ENERGY ELEMENT:", g

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                t_source = source(:,:,:,:,1)
            else
                t_source = 0
            end if
            if (allocated(bdysrc)) then
                do n = 1, 2
                    t_bdysrc(:,:,:,:,n) = fEg(n,1)*bdysrc
                end do
            end if
        else
            allocate(t_Sigma, source = Sigma(:,:,1:g-1,g,:,:))
            call SN_FEXS_K &
                (wL, cosmat, qPa, qfactorialmat, t_Sigma, psi, t_source)
            if (allocated(source)) then
                t_source = t_source + source(:,:,:,:,g)
            end if
            if (allocated(bdysrc)) then
                do n = 1, 2
                    t_bdysrc(:,:,:,:,n) = fEg(n,g)*bdysrc
                end do
            end if
            deallocate(t_Sigma)
        end if
        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        t_Sigmagg = Sigma(:,:,g,g,:,:)
        t_Sigmat = Sigmat(:,:,g,:)
        t_phi = 0

        ompstart = omp_get_wtime()
        call SN_FEXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()

        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        call SN_FEXS_SI &
            (convergence, pmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, &
            esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, t_Sigmagg, t_source, t_bdysrc, Tprimeinv, total, t_phi, t_psi)

        if (g .eq. 1) then
            psi(:,:,:,1:2) = t_psi
            do k = 1, NKe
                do e = 1, NE
                    phi(e,k,1:2) = t_phi(Cekk(e,k),1:2)
                end do
            end do
        else
            psi(:,:,:,g) = 0.5*psi(:,:,:,g) + 0.5*t_psi(:,:,:,1)
            psi(:,:,:,g+1) = t_psi(:,:,:,2)
            do k = 1, NKe
                do e = 1, NE
                    phi(e,k,g) = 0.5*phi(e,k,g) + 0.5*t_phi(Cekk(e,k),1)
                    phi(e,k,g+1) = t_phi(Cekk(e,k),2)
                end do
            end do
        end if

        ompfinish = omp_get_wtime()

        print *, "CONVERGED WITH:", total, "ITERATES"
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    if (particle .eq. 1) then
        call SN_write_field("photon", phi, psi)
    else if (particle .eq. 2) then
        call SN_write_field("electron", phi, psi)
    end if
end subroutine SN_FEXS_EI_SI

subroutine SN_FEXS_construct_inverted_transport_matrix &
    (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat, Tprimeinv)
    implicit none
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    integer, dimension(:,:,:), intent(in) :: enew
    real, dimension(:,:,:,:), intent(in) :: varsigmaup
    real, dimension(:,:,:), intent(in) :: Sigmat

    real, dimension(:,:,:,:,:,:,:), intent(inout) :: Tprimeinv ! (kp,k,np,n,enew,i,j)

    integer :: dir, e, f, i, j, k, kp, n, np, mat
    real, dimension(:,:,:,:), allocatable :: temp_Tprimeinv
    real, dimension(:,:), allocatable :: ident
    real, dimension(:,:), allocatable :: tempZ
    real, dimension(:,:), allocatable :: tempY
    real, dimension(:,:), allocatable :: inv1
    real, dimension(:,:), allocatable :: inv2
    real, dimension(:,:), allocatable :: inv2inv1
    real, dimension(:,:), allocatable :: inv1inv2

    allocate(temp_Tprimeinv(NKe,NKe,2,2))
    allocate(ident(NKe,NKe))
    allocate(tempZ(NKe,NKe))
    allocate(tempY(NKe,NKe))
    allocate(inv1(NKe,NKe))
    allocate(inv2(NKe,NKe))
    allocate(inv2inv1(NKe,NKe))
    allocate(inv1inv2(NKe,NKe))

    ident = identity(NKe)

    !$OMP parallel do collapse(3) shared(Tprimeinv) private(dir, e, f, i,j, k, kp, n, np, mat) &
    !$OMP private(inv1, inv2, tempZ, tempY, inv2inv1, inv1inv2, temp_Tprimeinv)
    do j = 1, Nphi
        do i = 1, Nmu
            do e = 1, NE
                mat = eltomat(e)
                ! NOTE: CAN STORE F+G. THEN WOULD DO:
                !tempZ = alpha + Sigmat(1,1,eltomat(e))*ident
                !tempY = alpha + Sigmat(2,2,eltomat(e))*ident
                inv1 = 0
                inv2 = 0
                do f = 1, NFe ! Fterm
                    inv1 = inv1 + varsigmaup(f,e,i,j)*I1invI2(1:NKe,1:NKe,f,e)
                end do
                do dir = 1, 3 ! Gterm
                    inv2 = inv2 - I1invI3vector(1:NKe,1:NKe,dir,e)*khat(dir,i,j)
                end do
                tempZ = inv1 + inv2 + Sigmat(1,1,mat)*ident ! F + G + M11*I, or A
                tempY = inv1 + inv2 + Sigmat(2,2,mat)*ident ! F + G + M22*I, or D

                    !verification(1:NKe,1:NKe,1,1) = tempZ
                    !verification(1:NKe,1:NKe,2,2) = tempY

                inv1 = explicit_4x4_matinv(tempZ) ! (F + G + M11)^(-1)
                tempY = tempY - Sigmat(2,1,mat)*Sigmat(1,2,mat)*inv1
                inv2 = explicit_4x4_matinv(tempY) ! (D-CAinvB)^(-1)
                inv2inv1 = matmul(inv2,inv1)
                inv1inv2 = matmul(inv1,inv2)

                tempZ = tempZ - Sigmat(2,1,mat)*Sigmat(1,2,mat)*inv2
                temp_Tprimeinv(1:NKe,1:NKe,1,1) = inv1 + Sigmat(2,1,mat)*Sigmat(1,2,mat)*matmul(inv1,inv2inv1)
                temp_Tprimeinv(1:NKe,1:NKe,1,2) = -Sigmat(2,1,mat)*inv1inv2
                temp_Tprimeinv(1:NKe,1:NKe,2,1) = -Sigmat(1,2,mat)*inv2inv1
                temp_Tprimeinv(1:NKe,1:NKe,2,2) = inv2

                    !verification(1:NKe,1:NKe,1,2) = test(1,2,g,1)*ident ! M12*I, or B
                    !verification(1:NKe,1:NKe,2,1) = test(2,1,g,1)*ident ! M21*I, or C
                    !verify_matmul = 0
                    !verify_matmul(:,:,1,1) = matmul(verification(:,:,1,1),temp_Tprimeinv(:,:,1,1)) + &
                    !    matmul(verification(:,:,1,2),temp_Tprimeinv(:,:,2,1))
                    !verify_matmul(:,:,1,2) = matmul(verification(:,:,1,1),temp_Tprimeinv(:,:,1,2)) + &
                    !    matmul(verification(:,:,1,2),temp_Tprimeinv(:,:,2,2))
                    !verify_matmul(:,:,2,1) = matmul(verification(:,:,2,1),temp_Tprimeinv(:,:,1,1)) + &
                    !    matmul(verification(:,:,2,2),temp_Tprimeinv(:,:,2,1))
                    !verify_matmul(:,:,2,2) = matmul(verification(:,:,2,1),temp_Tprimeinv(:,:,1,2)) + &
                    !    matmul(verification(:,:,2,2),temp_Tprimeinv(:,:,2,2))
                    !
                    !if (any(abs(verify_matmul - verify_ident) > 1.0E-7)) then
                    !    print *, [g,j,i,e]
                    !    print *, abs(verify_matmul - verify_ident)
                    !    print *, "BREAK"
                    !    print *, verify_matmul
                    !    stop
                    !end if

                do k = 1, NKe
                    do kp = 1, NKe
                        Tprimeinv(kp,k,:,:,enew(e,i,j),i,j) = transpose(temp_Tprimeinv(k,kp,:,:))
                    end do
                end do
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine SN_FEXS_construct_inverted_transport_matrix

subroutine SN_FEXS_SI &
    (l_convergence, pmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, &
    fold, Ia, Ja, noJa, varsigmadown, I1invI2f, Sigma, source, bdysrc, Tprimeinv, total, phi, psi)
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
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (np,n,l+1,mat)
    real, dimension(:,:,:,:,:), allocatable, intent(in) :: bdysrc
    real, dimension(:,:,:,:), intent(inout) :: source ! (k,i,j,n)
    real, dimension(:,:,:,:,:,:,:), intent(in) :: Tprimeinv ! (kp,k,np,n,enew,i,j)

    integer, intent(inout) :: total
    real, dimension(:,:), intent(inout) :: phi ! (k,n) ! phiu initially, leaves as phic + phiu
    real, dimension(:,:,:,:), intent(inout) :: psi ! (k,i,j,n)

    integer :: j, n, p
    real, dimension(:,:,:), allocatable :: jterm
    real, dimension(:,:), allocatable :: tempE
    real :: prevnorm
    real :: nextnorm
    real :: ompstart, ompfinish

    allocate(jterm(NK,Nphi,2))
    allocate(tempE(NK,2))

    ompstart = omp_get_wtime()

    prevnorm = norm2(phi)

    call SN_FEXS_sweep &
        (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
        I1invI2f, bdysrc, Tprimeinv, source)

    psi = source
    !$OMP parallel do collapse(2) shared(jterm) private(n, j)
    do n = 1, 2
        do j = 1, Nphi
            jterm(1:NK,j,n) = twopi*matmul(source(1:NK,1:Nmu,j,n),wL)/Nphi
        end do
    end do
    !$OMP end parallel do
    tempE = sum(jterm,2)
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

        call SN_FEXS_Kgg &
            (wL, cosmat, qPa, qfactorialmat, Sigma, source)

        call SN_FEXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, bdysrc, Tprimeinv, source)

        psi = psi + source
        !$OMP parallel do collapse(2) shared(jterm) private(n, j)
        do n = 1, 2
            do j = 1, Nphi
                jterm(1:NK,j,n) = twopi*matmul(source(1:NK,1:Nmu,j,n),wL)/Nphi
            end do
        end do
        !$OMP end parallel do
        tempE = sum(jterm,2)
        phi = phi + tempE

        nextnorm = norm2(tempE)

        total = p

        ompfinish = omp_get_wtime()
        if (nextnorm .eq. 0.0) then
            print *, "ITERATE:", p, "CONVERGENCE:", 0.0, "TIME (s):", ompfinish - ompstart
            return
        end if
        print *, "ITERATE:", p, "CONVERGENCE:", nextnorm/prevnorm, "TIME (s):", ompfinish - ompstart
        if (nextnorm/prevnorm .le. l_convergence) return
        if (prevnorm .eq. 0.0) return ! Is this necessary given that this conditional was already invoked??
        if (isnan(nextnorm/prevnorm)) then
            print *, "NAN"
            stop
        end if
    end do

    print *, "CONVERGENCE NOT ACHIEVED WITH", pmax, "ITERATES."
end subroutine SN_FEXS_SI

subroutine SN_FEXS_K &
    (wL, cosmat, qPa, qfactorialmat, Sigma, scatterers, vector)
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat ! (j,jp,m+1)
    real, dimension(:,:), intent(in) :: qPa ! (i,q)
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:,:,:), intent(in) :: Sigma ! (np,n,gp,l+1,mat)
    real, dimension(:,:,:,:), intent(in) :: scatterers ! (k,i,j,g)

    real, dimension(:,:,:,:), intent(inout) :: vector ! (k,i,j,n)

    integer :: i, ip, j, jp, k, l, mat, n, np, q, g, gpr
    integer :: NQ
    integer :: group
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:,:), allocatable :: Amat ! (k,n,ip,jp,l)
    real, dimension(:,:,:,:), allocatable :: Bmat ! (k,jp,q,n)
    real, dimension(:,:,:,:), allocatable :: Cmat ! (k,jp,q,n)

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    allocate(Amat(NK,2,Nmu,Nphi,NL+1))
    allocate(Bmat(NK,Nphi,NQ,2))
    allocate(Cmat(NK,Nphi,NQ,2))

    group = size(Sigma,3)
    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    Amat = 0
    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(4) shared(Amat) private(l, jp, ip, n, np, gpr)
        do l = 0, NL
            do jp = 1, Nphi
                do ip = 1, Nmu
                    do n = 1, 2
                        do gpr = 1, group
                            do np = 1, 2
                                Amat(1:NK,n,ip,jp,l+1) = Amat(1:NK,n,ip,jp,l+1) + &
                                twopi*scatterers(1:NK,ip,jp,gpr+np-1)*Sigma(np,n,gpr,l+1,1)/Nphi
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        do l = 0, NL
            do jp = 1, Nphi
                do ip = 1, Nmu
                    do n = 1, 2
                        do gpr = 1, group
                            do np = 1, 2
                                do mat = 1, Nmats
                                    do k = 1, NKmat(mat)
                                        Amat(nodesinmat(k,mat),n,ip,jp,l+1) = Amat(nodesinmat(k,mat),n,ip,jp,l+1) + &
                                            twopi*Nelnodes(nodesinmat(k,mat),mat)*&
                                            scatterers(nodesinmat(k,mat),ip,jp,gpr+np-1)*Sigma(np,n,gpr,l+1,mat)/&
                                            (Nphi*Nglobal(nodesinmat(k,mat)))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if

    !$OMP parallel do collapse(2) shared(Bmat) private(n, q, jp)
    do n = 1, 2
        do q = 1, NQ
            do jp = 1, Nphi
                Bmat(1:NK,jp,q,n) = matmul(Amat(1:NK,n,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
            end do
            Bmat(1:NK,1:Nphi,q,n) = qfactorialmat(q)*Bmat(1:NK,1:Nphi,q,n)
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(2) shared(Cmat) private(n, q)
    do n = 1, 2
        do q = 1, NQ
            Cmat(1:NK,1:Nphi,q,n) = matmul(Bmat(1:NK,1:Nphi,q,n),cosmat(1:Nphi,1:Nphi,emm(q)+1))
        end do
    end do
    !$OMP end parallel do

    if (size(Sigma,4) .eq. NL+2) then
        if (Nmats .eq. 1) then
            vector = 0
            !$OMP parallel do collapse(3) shared(vector) private(n, i, j, np, gpr)
            do n = 1, 2
                do j = 1, Nphi
                    do i = 1, Nmu
                        do gpr = 1, group
                            do np = 1, 2
                                vector(1:NK,i,j,n) = vector(1:NK,i,j,n) + &
                                    scatterers(1:NK,i,j,gpr+np-1)*Sigma(np,n,gpr,NL+2,1)
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        else
            vector = 0
            do n = 1, 2
                do j = 1, Nphi
                    do i = 1, Nmu
                        do gpr = 1, group
                            do np = 1, 2
                                do mat = 1, Nmats
                                    do k = 1, NKmat(mat)
                                        vector(nodesinmat(k,mat),i,j,n) = vector(nodesinmat(k,mat),i,j,n) + &
                                            Nelnodes(nodesinmat(k,mat),mat)*&
                                            scatterers(nodesinmat(k,mat),i,j,gpr+np-1)*Sigma(np,n,gpr,NL+2,mat)/&
                                            Nglobal(nodesinmat(k,mat))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end if
    else
        vector = 0
    end if

    !$OMP parallel do shared(vector) private(i, q)
    do i = 1, Nmu
        do q = 1, NQ
            vector(1:NK,i,1:Nphi,1:2) = vector(1:NK,i,1:Nphi,1:2) + Cmat(1:NK,1:Nphi,q,1:2)*qPa(i,q)
        end do
    end do
    !$OMP end parallel do
end subroutine SN_FEXS_K

subroutine SN_FEXS_Kgg &
    (wL, cosmat, qPa, qfactorialmat, Sigma, vector)
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat ! (j,jp,m+1)
    real, dimension(:,:), intent(in) :: qPa ! (i,q)
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (np,n,l,mat)

    real, dimension(:,:,:,:), intent(inout) :: vector ! (k,i,j,n)

    integer :: i, ip, j, jp, k, l, mat, n, np, q
    integer :: NQ
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:,:), allocatable :: Amat ! (k,n,i,j,l)
    real, dimension(:,:,:,:), allocatable :: Bmat ! (j,k,q,n)
    real, dimension(:,:,:,:), allocatable :: Cmat ! (j,k,q,n)
    real, dimension(:,:,:,:), allocatable :: gold
    real, dimension(Nphi) :: tempA

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    allocate(Amat(NK,2,Nmu,Nphi,NL+1))
    allocate(Bmat(NK,Nphi,NQ,2))
    allocate(Cmat(NK,Nphi,NQ,2))

    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    Amat = 0
    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(2) shared(Amat) private(l, n, np)
        do l = 0, NL
            do n = 1, 2
                do np = 1, 2
                    Amat(1:NK,n,1:Nmu,1:Nphi,l+1) = Amat(1:NK,n,1:Nmu,1:Nphi,l+1) + &
                    twopi*Sigma(np,n,l+1,1)*vector(1:NK,1:Nmu,1:Nphi,np)/Nphi
                end do
            end do
        end do
        !$OMP end parallel do
    else
        do l = 0, NL
            do n = 1, 2
                do np = 1, 2
                    do mat = 1, Nmats
                        do k = 1, NKmat(mat)
                            Amat(nodesinmat(k,mat),n,1:Nmu,1:Nphi,l+1) = Amat(nodesinmat(k,mat),n,1:Nmu,1:Nphi,l+1) + &
                            twopi*Nelnodes(nodesinmat(k,mat),mat)*&
                            Sigma(np,n,l+1,mat)*vector(nodesinmat(k,mat),1:Nmu,1:Nphi,np)/&
                            (Nphi*Nglobal(nodesinmat(k,mat)))
                        end do
                    end do
                end do
            end do
        end do
    end if

    !$OMP parallel do collapse(2) shared(Bmat) private(n, q, jp)
    do n = 1, 2
        do q = 1, NQ
            do jp = 1, Nphi
                Bmat(1:NK,jp,q,n) = matmul(Amat(1:NK,n,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
            end do
            Bmat(1:NK,1:Nphi,q,n) = qfactorialmat(q)*Bmat(1:NK,1:Nphi,q,n)
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(2) shared(Cmat) private(n, q)
    do n = 1, 2
        do q = 1, NQ
            Cmat(1:NK,1:Nphi,q,n) = matmul(Bmat(1:NK,1:Nphi,q,n),cosmat(1:Nphi,1:Nphi,emm(q)+1))
        end do
    end do
    !$OMP end parallel do

    if (size(Sigma,3) .eq. NL+2) then
        if (Nmats .eq. 1) then
            !!$OMP parallel do collapse(2) shared(vector) private(i, j)
            do j = 1, Nphi
                do i = 1, Nmu
                    vector(1:NK,i,j,1:2) = matmul(vector(1:NK,i,j,1:2),Sigma(1:2,1:2,NL+2,1))
                    ! NOTE: This is not a mistake. But, is this slow since it's rewriting vector?
                end do
            end do
            !!$OMP end parallel do
        else
            print *, "SN_FEXS_Kgg DELTA DOWN MATMUL IS WIP. FIX THIS!"
            stop
        end if
    else
        vector = 0
    end if

    !$OMP parallel do collapse(2) shared(vector) private(n, i, q)
    do n = 1, 2
        do i = 1, Nmu
            do q = 1, NQ
                vector(1:NK,i,1:Nphi,n) = vector(1:NK,i,1:Nphi,n) + Cmat(1:NK,1:Nphi,q,n)*qPa(i,q)
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine SN_FEXS_Kgg

subroutine SN_FEXS_sweep &
    (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, I1invI2f, &
    bdysrc, Tprimeinv, vector)
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
    real, dimension(:,:,:,:,:), allocatable, intent(in) :: bdysrc
    real, dimension(:,:,:,:,:,:,:), intent(in) :: Tprimeinv

    real, dimension(:,:,:,:), intent(inout) :: vector

    integer :: e, f, i, ip, j, k, n, np
    integer, dimension(:), allocatable :: l_esweeplist
    integer, dimension(:), allocatable :: l_enew
    integer, dimension(:,:), allocatable :: l_fold
    integer, dimension(:), allocatable :: l_Ia
    integer, dimension(:), allocatable :: l_Ja
    real, dimension(:), allocatable :: l_varsigmadown
    real, dimension(NKe,2) :: tempA
    real, dimension(:,:,:), allocatable :: l_vector
    integer :: SS1
    logical :: use_bdy
    integer :: NBE

    allocate(l_esweeplist(NE))
    allocate(l_enew(NE))
    allocate(l_fold(NFe,NE))
    allocate(l_Ia(NE+1))
    allocate(l_Ja(noJa))
    allocate(l_varsigmadown(noJa))
    allocate(l_vector(NKe,2,NE))

    use_bdy = allocated(bdysrc)

    if (use_bdy) then
        NBE = size(bdysrc,2)
    end if

    !$OMP parallel do collapse(2) shared(vector) private(j, i, n, k, ip, f, e) &
    !$OMP private(l_esweeplist, l_enew, l_vector, l_fold, l_Ia, l_Ja, l_varsigmadown, SS1, tempA)
    do j = 1, Nphi
        do i = 1, Nmu
            l_esweeplist = esweeplist(1:NE,i,j) ! Don't use local arrays with parallelization?
            l_enew = enew(1:NE,i,j)

            do k = 1, NKe
                do e = 1, NE
                    l_vector(k,1:2,l_enew(e)) = vector(Cekk(e,k),i,j,1:2)
                end do
            end do

            if (use_bdy) then
                do e = 1, NBE
                    l_vector(1:NKe,1:2,l_enew(bdyel(e))) = l_vector(1:NKe,1:2,l_enew(bdyel(e))) + &
                        bdysrc(1:NKe,e,i,j,1:2)
                end do
            end if

            l_fold = fold(1:NFe,1:NE,i,j)
            l_Ia = Ia(1:NE+1,i,j)
            l_Ja = Ja(1:noJa,i,j)
            l_varsigmadown = varsigmadown(1:noJa,i,j)

            SS1 = esweepbounds(2,i,j)
            do ip = 1, SS1
                tempA(1:NKe,1:2) = l_vector(1:NKe,1:2,ip)

                do n = 1, 2
                    do k = 1, NKe
                        l_vector(k,n,ip) = sum(Tprimeinv(1:NKe,k,1:2,n,ip,i,j)*tempA)
                    end do
                end do
            end do

            do ip = SS1+1, NE
                tempA(1:NKe,1:2) = l_vector(1:NKe,1:2,ip)
                do n = 1, 2
                    do f = 1, l_Ia(ip+1)-l_Ia(ip)
                        do k = 1, NKe
                            tempA(k,n) = tempA(k,n) - l_varsigmadown(l_Ia(ip)-1+f)*&
                            sum(I1invI2f(1:NKe,k,l_fold(f,ip),l_esweeplist(ip))*l_vector(1:NKe,n,l_Ja(l_Ia(ip)-1+f))) ! Make this faster by using a temporary for l_vector(eprime)? I can know its size before f iteration
                            ! And/or use l_I1invI2f to avoid repeated indirect addressing??
                        end do
                    end do
                end do

                do n = 1, 2
                    do k = 1, NKe
                        l_vector(k,n,ip) = sum(Tprimeinv(1:NKe,k,1:2,n,ip,i,j)*tempA)
                    end do
                end do
            end do

            vector(1:NK,i,j,1:2) = 0
            do k = 1, NKe
                do e = 1, NE
                    vector(Cekk(l_esweeplist(e),k),i,j,1:2) = vector(Cekk(l_esweeplist(e),k),i,j,1:2) + l_vector(k,1:2,e) ! MAKE THIS FASTER!!
                    ! Apparently doing addition like this is slow.
                end do
            end do
            do n = 1, 2
                vector(1:NK,i,j,n) = vector(1:NK,i,j,n)/Nglobal
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine SN_FEXS_sweep

! PN - general
subroutine PN_read_field &
    (label, l_G, l_NL, phiu, phi, psi)
    implicit none
    character(*), intent(in) :: label
    integer, intent(in) :: l_G
    integer, intent(in) :: l_NL

    real, dimension(:,:,:), allocatable, intent(inout) :: phiu
    real, dimension(:,:,:), allocatable, intent(inout) :: phi
    real, dimension(:,:,:), allocatable, intent(inout) :: psi

    character(80) :: file_name

    if (FCS .and. &
        ((label .eq. "electron" .and. &
        transport_mode .eq. "external electron beam") .or. &
        (label .eq. "photon" .and. &
        (transport_mode .eq. "external photon beam" .or. &
        transport_mode .eq. "external photon beam coupled")))) allocate(phiu(NE,NKe,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(psi(NK,(l_NL+1)**2,l_G))

    if (FCS .and. &
        ((label .eq. "electron" .and. &
        transport_mode .eq. "external electron beam") .or. &
        (label .eq. "photon" .and. &
        (transport_mode .eq. "external photon beam" .or. &
        transport_mode .eq. "external photon beam coupled")))) then
        file_name = "Stored solutions/"//label//"_uncollided_fluence.dat"

        open(1, file = file_name, form = "unformatted", action = "read")
        read(1) phiu
        close(1)
    end if

    file_name = "Stored solutions/"//label//"_fluence.dat"

    open(1, file = file_name, form = "unformatted", action = "read")
    read(1) phi
    close(1)

    file_name = "Stored solutions/"//label//"_angular_fluence.dat"

    open(1, file = file_name, form = "unformatted", action = "read")
    read(1) psi
    close(1)
end subroutine PN_read_field

    ! Had the fnames overwrite eachother. This may have screwed things up earlier... Phi file was used for both
subroutine PN_write_field &
    (label, phi, psi)
    implicit none
    character(*), intent(in) :: label
    real, dimension(:,:,:), intent(in) :: phi
    real, dimension(:,:,:), intent(in) :: psi

    character(80) :: file_name

    file_name = trim(adjustl(output_fname))//"/"//label//"_fluence.dat"

    open(1, file = file_name, form = "unformatted")
    write(1) phi
    close(1)

    file_name = trim(adjustl(output_fname))//"/"//label//"_angular_fluence.dat"

    open(1, file = file_name, form = "unformatted")
    write(1) psi
    close(1)
end subroutine PN_write_field

! PN - MGXS
!! FULLY UNOPTIMIZED
subroutine PN_MGXS_prepare_source &
    (l_NL, Cekk, Sigma, source, scatterers, e_source)
    implicit none
    integer, intent(in) :: l_NL
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:), intent(in) :: Sigma
    real, dimension(:,:), intent(inout) :: source
    real, dimension(:,:,:), intent(in) :: scatterers

    real, dimension(:,:,:), intent(inout) :: e_source

    integer :: e, k

    call PN_MGXS_K(l_NL, Sigma, scatterers, source)

    do k = 1, NKe
        do e = 1, NE
            e_source(e,k,:) = source(Cekk(e,k),:)
        end do
    end do
end subroutine PN_MGXS_prepare_source

subroutine PN_MGXS_GMRESm &
    (l_NL, l_convergence, pmax, mmax, eprime, Ceff, sgn, Awl, Aup, Adn, &
    I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, source, total, vector)
    implicit none
    integer, intent(in) :: l_NL
    real :: l_convergence
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    integer, dimension(:,:), intent(in) :: eprime
    integer, dimension(:,:), intent(in) :: Ceff
    integer, dimension(:,:), intent(in) :: sgn
    real, dimension(:,:,:), intent(in) :: Awl
    real, dimension(:,:,:), intent(in) :: Aup
    real, dimension(:,:,:), intent(in) :: Adn
    !real, dimension(:,:,:,:), intent(in) :: Aup
    !real, dimension(:,:,:,:), intent(in) :: Adn
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:,:), intent(in) :: Sigma
    real, dimension(:), intent(in) :: Sigmat
    real, dimension(:,:,:), intent(inout) :: source

    integer, intent(inout) :: total
    real, dimension(:,:,:), intent(inout) :: vector

    integer :: e, f, i, j
    integer :: m
    integer :: max_iter
    integer :: NQ
    real, dimension(:,:,:), allocatable :: tempA
    real, dimension(:,:,:), allocatable :: tempB
    real, dimension(:,:,:), allocatable :: resid0
    real :: resid
    real, dimension(:,:,:,:), allocatable :: vec
    real, dimension(:,:), allocatable :: Qk
    integer, dimension(:,:), allocatable :: local_ident
    real, dimension(:,:), allocatable :: Hess
    real, dimension(:,:), allocatable :: Rk
    real :: denom
    real, dimension(:,:), allocatable :: Rot
    real, dimension(:,:), allocatable :: tempC
    real, dimension(:), allocatable :: gk
    real, dimension(:), allocatable :: tvec
    real, dimension(:), allocatable :: e1
    real, dimension(:), allocatable :: y
    real :: ompstart, ompfinish
    real, dimension(:,:,:,:,:,:), allocatable :: A
    real, dimension(:,:,:), allocatable :: storpsi
    real, dimension(:,:,:), allocatable :: storpsi1
    real, dimension(:,:,:,:), allocatable :: tA

    NQ = (l_NL+1)**2
    m = pmax
    max_iter = mmax

    allocate(A(NQ,NKe,NQ,NKe,NFe,NE))
    allocate(tA(NQ,NQ,NKe,NKe))
    A = 0
    do e = 1, NE
        do f = 1, NFe
            if (eprime(f,e) .eq. 0) cycle
            call random_number(tA)
            A(:,:,:,:,f,e) = tA
        end do
    end do
    allocate(storpsi(NE,NKe,NQ))
    call random_number(storpsi)
    !allocate(storpsi1, source = storpsi)

    call PN_MGXS_Lgg &
    (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, storpsi, A)

    source = storpsi
    print *, "HERE"

    ! WRITE TOTAL STORAGE AS FX OF NE, NKe, NQ, AND m!!!!

    NQ = (l_NL+1)**2

    allocate(tempA(NE,NKe,NQ))
    allocate(tempB(NE,NKe,NQ))
    allocate(resid0(NE,NKe,NQ))
    allocate(vec(NE,NKe,NQ,m))
    allocate(Qk(m+1,m+1))
    allocate(local_ident(m+1,m+1))
    allocate(Hess(m+1,m))
    allocate(Rk(m+1,m))
    allocate(Rot(2,2))
    allocate(tempC(2,m+1))
    allocate(gk(m+1))
    allocate(tvec(m+1))
    allocate(e1(m+1))
    allocate(y(m))

    local_ident = identity(m+1)
    e1 = real(0.0,8)
    e1(1) = real(1.0,8)

    tempA = vector

    call PN_MGXS_Lgg &
    (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA, A)

    resid0 = source - tempA

    GMRESm: do
        !ompstart = omp_get_wtime()
        resid = norm2(resid0)

        if (resid .lt. l_convergence) then
            print *, total, resid
            exit GMRESm
        end if

        vec(1:NE,1:NKe,1:NQ,1) = resid0/resid

        Qk = real(local_ident,8)

        !ompfinish = omp_get_wtime()

        !print *, "0:", ompfinish - ompstart, "seconds"

        do j = 1, m-1
            !ompstart = omp_get_wtime()
            !print *, j
            total = total + 1
            ! Gram-Schmidt orthogonalization to construct vec(:,j) vectors and Hess array
            tempA = vec(1:NE,1:NKe,1:NQ,j)

            call PN_MGXS_Lgg &
            (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA, A)
            !ompfinish = omp_get_wtime()

            !print *, j, "1:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            tempB = tempA
            !$OMP parallel do shared(Hess, tempB) private(i)
            do i = 1, j
                Hess(i,j) = sum(tempA*vec(1:NE,1:NKe,1:NQ,i))
                !tempB = tempB - Hess(i,j)*vec(1:NE,1:NKe,1:NQ,i)
            end do
            !$OMP end parallel do
            do i = 1, j
                tempB = tempB - Hess(i,j)*vec(1:NE,1:NKe,1:NQ,i)
            end do
            !ompfinish = omp_get_wtime()

            !print *, j, "2:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            Hess(j+1,j) = norm2(tempB)
            vec(1:NE,1:NKe,1:NQ,j+1) = tempB/Hess(j+1,j)
            !ompfinish = omp_get_wtime()

            !print *, j, "3:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            ! Append latest updated column of Hess to R, with cumulative rotations applied
            Rk(j+1,j) = Hess(j+1,j)
            Rk(1:j,j) = matmul(Qk(1:j,1:j),Hess(1:j,j)) ! Maybe take advantage of Q's structure? I think its Hessian itself.
            !ompfinish = omp_get_wtime()

            !print *, j, "4:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
            denom = hypot(Rk(j,j),Rk(j+1,j))
            Rot(:,1) = [Rk(j,j)/denom, -Rk(j+1,j)/denom]
            Rot(:,2) = [-Rot(2,1), Rot(1,1)]
            !ompfinish = omp_get_wtime()

            !print *, j, "5:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            Rk(j:j+1,j) = [denom, real(0.0,8)]
            ! Update cumulative rotation matrix
            tempC(1:2,1:j+1) = Qk(j:j+1,1:j+1)
            Qk(j:j+1,1:j+1) = matmul(Rot,tempC(1:2,1:j+1))
            !ompfinish = omp_get_wtime()

            !print *, j, "6:", ompfinish - ompstart, "seconds"

            !! Check for convergence
            print *, j, abs(resid*Qk(j+1,1))
            if (abs(resid*Qk(j+1,1)) <= l_convergence) then
                gk(1:j+1) = resid*Qk(1:j+1,1)
                do i = j, 1, -1
                    y(i) = (gk(i) - sum(Rk(i,i+1:j)*y(i+1:j)))/Rk(i,i)
                end do

                do i = 1, j ! THIS STEP MAY BE EASIER IF RESHAPES AND PACKS ARE USED
                    vector = vector + vec(1:NE,1:NKe,1:NQ,i)*y(i)
                end do

                exit GMRESm
            end if
            !if (isnan(resid*Qk(j+1,1))) stop
        end do

        !ompstart = omp_get_wtime()

        total = total + 1

        ! Do mth step differently
        tempA = vec(1:NE,1:NKe,1:NQ,m)

        call PN_MGXS_Lgg &
        (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA, A)

        !ompfinish = omp_get_wtime()

        !print *, "7:", ompfinish - ompstart, "seconds"

        !ompstart = omp_get_wtime()
        !$OMP parallel do shared(Hess) private(i)
        do i = 1, m
            Hess(i,m) = sum(tempA*vec(1:NE,1:NKe,1:NQ,i))
        end do
        !$OMP end parallel do
        Hess(m+1,m) = sqrt(norm2(tempA)**2 - sum(Hess(1:m,m)**2))
        !ompfinish = omp_get_wtime()

        !print *, "8:", ompfinish - ompstart, "seconds"

        !if (any(isnan(Hess))) then
        !    print *, Hess(m+1,m)
        !    print *, norm2(tempA)**2
        !    print *, sum(Hess(1:m,m)**2)
        !    print *, norm2(tempA)**2 - sum(Hess(1:m,m)**2)
        !    stop
        !end if

        ! Append latest updated column of Hess to R, with cumulative rotations applied
        Rk(m+1,m) = Hess(m+1,m)
        Rk(1:m,m) = matmul(Qk(1:m,1:m),Hess(1:m,m)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

        ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
        denom = hypot(Rk(m,m),Rk(m+1,m))
        Rot(:,1) = [Rk(m,m)/denom, -Rk(m+1,m)/denom] ! Should I transpose these? were they supposed to be transposed? (~3/31/23)
        Rot(:,2) = [-Rot(2,1), Rot(1,1)] !

        Rk(m:m+1,m) = [denom, real(0.0,8)]

        ! Update cumulative rotation matrix
        tempC(1:2,1:m+1) = Qk(m:m+1,1:m+1)
        Qk(m:m+1,1:m+1) = matmul(Rot,tempC(1:2,1:m+1))

        ! Construct current solution
        gk(1:m+1) = resid*Qk(1:m+1,1)
        do i = m, 1, -1
            y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
        end do

        do i = 1, m
            vector = vector + vec(1:NE,1:NKe,1:NQ,i)*y(i)
        end do

        tempA = vector
        call PN_MGXS_Lgg &
        (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA, A)

        resid0 = source - tempA

        !! Slower??
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

        if (abs(gk(m+1)) .le. convergence) then
            exit GMRESm
        end if

        if (total >= max_iter) then
            gk(1:m+1) = resid*Qk(1:m+1,1)
            do i = m, 1, -1
                y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
            end do

            do i = 1, m
                vector = vector + vec(1:NE,1:NKe,1:NQ,i)*y(i)
            end do
            print *, "GMRESm: Did not converge within", max_iter, "iterations."
            print *, "Solution being given with:", abs(gk(j+1)), "residual."
            !print *, "iterations:", total
            !print *, "residual 2:", abs(gk(j+1))
            exit GMRESm
            !stop
        end if
    end do GMRESm

    stop
end subroutine PN_MGXS_GMRESm

subroutine PN_MGXS_K &
    (l_NL, Sigma, scatterers, vector)
    implicit none
    integer, intent(in) :: l_NL
    real, dimension(:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:), intent(in) :: scatterers

    real, dimension(:,:), intent(inout) :: vector ! DOESN'T SET IT TO ZERO. MAKE CONSISTEN W/ SN?

    integer :: g, k, l, mat
    integer :: group

    group = size(Sigma,1)

    if (Nmats .eq. 1) then
        !$OMP parallel do shared(vector) private(g, l)
        do l = 0, l_NL
            do g = 1, group
                vector(1:NK,l*(l+1)-l+1:l*(l+1)+l+1) = vector(1:NK,l*(l+1)-l+1:l*(l+1)+l+1) + &
                    (fourpi*Sigma(g,l+1,1)/(2*l+1))*scatterers(1:NK,l*(l+1)-l+1:l*(l+1)+l+1,g)
            end do
        end do
        !$OMP end parallel do
    else
        do l = 0, l_NL
            do g = 1, group
                do mat = 1, Nmats
                    do k = 1, NKmat(mat)
                        vector(nodesinmat(k,mat),l*(l+1)-l+1:l*(l+1)+l+1) = &
                            vector(nodesinmat(k,mat),l*(l+1)-l+1:l*(l+1)+l+1) + &
                            Nelnodes(nodesinmat(k,mat),mat)*(fourpi*Sigma(g,l+1,mat)/(2*l+1))*&
                            scatterers(nodesinmat(k,mat),l*(l+1)-l+1:l*(l+1)+l+1,g)/&
                            Nglobal(nodesinmat(k,mat))
                    end do
                end do
            end do
        end do
    end if
end subroutine PN_MGXS_K

subroutine PN_MGXS_Lgg &
    (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, vector, A)
    implicit none
    ! USES METHOD 1 WHICH SCALES BETTER WITH NE BUT NOT WELL AT ALL WITH NQ
    integer, intent(in) :: l_NL
    integer, dimension(:,:), intent(in) :: eprime ! (f,e)
    integer, dimension(:,:), intent(in) :: Ceff
    integer, dimension(:,:), intent(in) :: sgn
    real, dimension(:,:,:), intent(in) :: Awl ! (q,qp,dir) for Method 1, UNDECIDED for Method 2
    real, dimension(:,:,:), intent(in) :: Aup
    real, dimension(:,:,:), intent(in) :: Adn
    !real, dimension(:,:,:,:), intent(in) :: Aup ! (e,q,f,q) for Method 1, (q,q,f,e) for Method 2
    !real, dimension(:,:,:,:), intent(in) :: Adn ! (e,q,f,q) for Method 1, (q,q,f,e) for Method 2
    real, dimension(:,:,:,:), intent(in) :: I1invI2 ! (e,kp,f,k) for Methods 1 and 2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f ! (e,kp,f,k) for Methods 1 and 2
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector ! (e,kp,dir,k) for Method 1, (dir,e,k,kp) for Method 2
    real, dimension(:,:), intent(in) :: Sigma ! (l+1,mat)
    real, dimension(:), intent(in) :: Sigmat ! (mat)
    real, dimension(:,:,:,:,:,:), intent(in) :: A

    real, dimension(:,:,:), intent(inout) :: vector

    integer :: dir, e, ep, f, k, kp, l, mat, q, qp
    integer :: NQ
    !integer, dimension(:), allocatable :: ell
    real, dimension(:,:,:), allocatable :: Mterm
    real, dimension(:,:,:,:), allocatable :: Amat
    !real :: ompstart, ompfinish
    !real :: sumtest

    NQ = (l_NL+1)**2
    allocate(Mterm(NE,NKe,NQ))
    Mterm = 0
    do e = 1, NE
        do f = 1, NFe
            if (eprime(f,e) .eq. 0) cycle
            !$OMP parallel do collapse(2) shared(Mterm) private(e,k,q,f)
            do k = 1, NKe
                do q = 1, NQ
                    Mterm(e,k,q) = sum(A(:,:,q,k,f,e)*&
                        vector(eprime(f,e),:,:))
                end do
            end do
            !$OMP end parallel do
        end do
    end do

    vector = Mterm

    return
    !allocate(ell(NQ))

    !do q = 1, NQ
    !    ell(q) = ceiling(sqrt(real(q)))-1
    !end do

    !MTERM = 0 ! Remove soon

    ! MTERM
    if (Nmats .eq. 1) then
        !$OMP parallel do shared(Mterm) private(l)
        do l = 0, l_NL
            Mterm(1:NE,1:NKe,l*(l+1)-l+1:l*(l+1)+l+1) = &
                (Sigmat(1) - fourpi*Sigma(l+1,1)/(2*l+1))*vector(1:NE,1:NKe,l*(l+1)-l+1:l*(l+1)+l+1)
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do shared(Mterm) private(e, l)
        do l = 0, l_NL
            do e = 1, NE
                Mterm(e,1:NKe,l*(l+1)-l+1:l*(l+1)+l+1) = &
                    (Sigmat(eltomat(e)) - fourpi*Sigma(l+1,eltomat(e))/(2*l+1))*&
                    vector(e,1:NKe,l*(l+1)-l+1:l*(l+1)+l+1)
            end do
        end do
        !$OMP end parallel do
    end if

    ! GTERM
    !ompstart = omp_get_wtime()
    allocate(Amat(NE,NQ,NKe,3))
    !$OMP parallel do collapse(2) shared(Amat) private(dir, k)
    do dir = 1, 3
        do k = 1, NKe
            Amat(1:NE,1:NQ,k,dir) = matmul(vector(1:NE,k,1:NQ),Awl(1:NQ,1:NQ,dir))
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(2) shared(Mterm) private(dir, k, q) ! Need to do collapse(2) since its a sum
    do q = 1, NQ
        do k = 1, NKe
            do dir = 1, 3
                Mterm(1:NE,k,q) = Mterm(1:NE,k,q) - sum(I1invI3vector(1:NE,1:NKe,dir,k)*Amat(1:NE,q,1:NKe,dir),2)
            end do
        end do
    end do
    !$OMP end parallel do
    deallocate(Amat)

    !!!!!!$OMP parallel do collapse(3) shared(Mterm) private(k,e,dir) ! Use for large NQ?
    !!!!!do k = 1, NKe
    !!!!!    do e = 1, NE
    !!!!!        do dir = 1, 3
    !!!!!            Mterm(e,k,1:NQ) = Mterm(e,k,1:NQ) - &
    !!!!!            matmul(matmul(Awl(1:NQ,1:NQ,dir),transpose(vector(e,1:NKe,1:NQ))),I1invI3vector(e,k,1:NKe,dir))
    !!!!!        end do
    !!!!!    end do
    !!!!!end do
    !!!!!!$OMP end parallel do
    !!!!
    !!!!!!$OMP parallel do collapse(4) shared(Mterm) private(q, dir, kp, k)
    !!!!!do dir = 1, 3
    !!!!!    do q = 1, NQ
    !!!!!        do kp = 1, NKe
    !!!!!            do k = 1, NKe
    !!!!!                Mterm(1:NE,k,q) = Mterm(1:NE,k,q) - I1invI3vector(1:NE,k,kp,dir)*&
    !!!!!                    matmul(vector(1:NE,kp,Ja(Ia(q,dir):Ia(q+1,dir)-1,dir)),Awl(Ia(q,dir):Ia(q+1,dir)-1,dir))
    !!!!!            end do
    !!!!!        end do
    !!!!!    end do
    !!!!!end do
    !!!!!!$OMP end parallel do

    ! FTERM
    allocate(Amat(NE,NQ,NFe,NKe))
    Amat = 0
    !$OMP parallel do collapse(3) shared(Amat) private(k, f, q, kp)
    do k = 1, NKe
        do f = 1, NFe
            do q = 1, NQ
                do kp = 1, NKe
                    Amat(1:NE,q,f,k) = Amat(1:NE,q,f,k) + I1invI2(1:NE,kp,f,k)*vector(1:NE,kp,q)
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(2) shared(Mterm) private(q, k, f) ! Does this need to be collapse(2)?? It seems not, but I don't know why not.
    do q = 1, NQ
        do k = 1, NKe ! LOOPING OVER e MIGHT MAKE THIS BETTER FOR PARALLELIZATION ON SUPERCOMPUTER
            do f = 1, NFe
                Mterm(1:NE,k,q) = Mterm(1:NE,k,q) + &
                    sgn(1:NE,f)*sum(Amat(1:NE,1:NQ,f,k)*Aup(Ceff(1:NE,f),1:NQ,q),2) ! sum(Amat(1:NE,1:NQ,f,k)*Aup(1:NE,1:NQ,f,q),2)
            end do
        end do
    end do
    !$OMP end parallel do

    ! SWEPT FTERM
    Amat = 0 ! Because eprime(f,e) may be zero, so it will cycle out of some [f,e] pairs that will retain their values
    !$OMP parallel do collapse(3) shared(Amat) private(k, f, e)
    do k = 1, NKe
        do f = 1, NFe
            do e = 1, NE
                if (eprime(f,e) .eq. 0) cycle
                Amat(e,1:NQ,f,k) = matmul(transpose(vector(eprime(f,e),1:NKe,1:NQ)),I1invI2f(e,1:NKe,f,k)) ! USE A REORDERED eprime AND DO SBR STUFF AGAIN? as in, reorganize eprime(e,f) such that it's lower tri/more easily accesible? (~3/10/23)
            end do
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(2) shared(Mterm) private(q, k, f) ! Does this need to be collapse(2)?? It seems not, but I don't know why not.
    do q = 1, NQ
        do k = 1, NKe
            do f = 1, NFe
                Mterm(1:NE,k,q) = Mterm(1:NE,k,q) + &
                    sgn(1:NE,f)*sum(Amat(1:NE,1:NQ,f,k)*Adn(Ceff(1:NE,f),1:NQ,q),2) ! sum(Amat(1:NE,1:NQ,f,k)*Adn(1:NE,1:NQ,f,q),2)
            end do
        end do
    end do
    !$OMP end parallel do

    vector = Mterm
    ! CONSIDER PROBLEMS POSED BY ARRAY TEMPORARIES

    !stop
end subroutine PN_MGXS_Lgg

! PN - FEXS
!! FULLY UNOPTIMIZED
subroutine PN_FEXS_prepare_source &
    (l_NL, Cekk, Sigma, source, scatterers, e_source)
    implicit none
    integer, intent(in) :: l_NL
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:,:,:), intent(in) :: Sigma ! (np,n,g,l,mat)
    real, dimension(:,:,:), intent(inout) :: source ! (k,q,n)
    real, dimension(:,:,:), intent(in) :: scatterers ! (k,q,g)

    real, dimension(:,:,:,:), intent(inout) :: e_source ! (e,k,q,n)

    integer :: e, k

    call PN_FEXS_K(l_NL, Sigma, scatterers, source)

    do k = 1, NKe
        do e = 1, NE
            e_source(e,k,:,:) = source(Cekk(e,k),:,:)
        end do
    end do
end subroutine PN_FEXS_prepare_source

subroutine PN_FEXS_GMRESm &
    (l_NL, l_convergence, pmax, mmax, eprime, Ceff, sgn, Awl, Aup, Adn, &
    I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, source, total, vector)
    ! DOES IT FOR A SINGLE GROUP ONLY!!! This is a problem because I need it this way for accelerating SN, but as a standalone, it should do all groups. FIX by changing structure of everything
    implicit none
    integer, intent(in) :: l_NL
    real :: l_convergence
    integer, intent(in) :: pmax
    integer, intent(in) :: mmax
    integer, dimension(:,:), intent(in) :: eprime
    integer, dimension(:,:), intent(in) :: Ceff
    integer, dimension(:,:), intent(in) :: sgn
    real, dimension(:,:,:), intent(in) :: Awl
    real, dimension(:,:,:), intent(in) :: Aup
    real, dimension(:,:,:), intent(in) :: Adn
    real, dimension(:,:,:,:), intent(in) :: I1invI2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:), intent(in) :: Sigmat
    real, dimension(:,:,:,:), intent(in) :: source ! (e,k,q,n)

    integer, intent(inout) :: total
    real, dimension(:,:,:,:), intent(inout) :: vector ! (e,k,q,n)

    integer :: i, j
    integer :: m
    integer :: max_iter
    integer :: NQ
    real, dimension(:,:,:,:), allocatable :: tempA
    real, dimension(:,:,:,:), allocatable :: tempB
    real, dimension(:,:,:,:), allocatable :: resid0
    real :: resid
    real, dimension(:,:,:,:,:), allocatable :: vec
    real, dimension(:,:), allocatable :: Qk
    integer, dimension(:,:), allocatable :: local_ident
    real, dimension(:,:), allocatable :: Hess
    real, dimension(:,:), allocatable :: Rk
    real :: denom
    real, dimension(:,:), allocatable :: Rot
    real, dimension(:,:), allocatable :: tempC
    real, dimension(:), allocatable :: gk
    real, dimension(:), allocatable :: tvec
    real, dimension(:), allocatable :: e1
    real, dimension(:), allocatable :: y
    real :: ompstart, ompfinish

    ! WRITE TOTAL STORAGE AS FX OF NE, NKe, NQ, AND m!!!!

    NQ = (l_NL+1)**2
    m = pmax
    max_iter = mmax

    allocate(tempA(NE,NKe,NQ,2))
    allocate(tempB(NE,NKe,NQ,2))
    allocate(resid0(NE,NKe,NQ,2))
    allocate(vec(NE,NKe,NQ,2,m))
    allocate(Qk(m+1,m+1))
    allocate(local_ident(m+1,m+1))
    allocate(Hess(m+1,m))
    allocate(Rk(m+1,m))
    allocate(Rot(2,2))
    allocate(tempC(2,m+1))
    allocate(gk(m+1))
    allocate(tvec(m+1))
    allocate(e1(m+1))
    allocate(y(m))

    local_ident = identity(m+1)
    e1 = 0
    e1(1) = 1

    tempA = vector

    if (any(isnan(tempA))) then
        print *, "TEMPA 1"
        stop
    end if
    call PN_FEXS_Lgg &
    (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA)

    resid0 = source - tempA
    !if (any(isnan(resid0))) then
    !    print *, "RESID0"
    !    print *, any(isnan(source))
    !    print *, any(isnan(tempA))
    !    stop
    !end if
    GMRESm: do
        !ompstart = omp_get_wtime()
        resid = norm2(resid0)

        if (resid .lt. l_convergence) then
            print *, total, resid
            exit GMRESm
        end if

        vec(1:NE,1:NKe,1:NQ,1:2,1) = resid0/resid

        Qk = local_ident

        !ompfinish = omp_get_wtime()

        !print *, "0:", ompfinish - ompstart, "seconds"

        do j = 1, m-1
            !ompstart = omp_get_wtime()
            !print *, j
            total = total + 1
            ! Gram-Schmidt orthogonalization to construct vec(:,j) vectors and Hess array
            tempA = vec(1:NE,1:NKe,1:NQ,1:2,j)
            !if (any(isnan(tempA))) then
            !    print *, "TEMPA 2"
            !    stop
            !end if
            call PN_FEXS_Lgg &
            (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA)
            !ompfinish = omp_get_wtime()
            !if (any(isnan(tempA))) then
            !    print *, "TEMPA 3"
            !    stop
            !end if
            !print *, j, "1:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            tempB = tempA
            !$OMP parallel do shared(Hess, tempB) private(i)
            do i = 1, j
                Hess(i,j) = sum(tempA*vec(1:NE,1:NKe,1:NQ,1:2,i))
                !tempB = tempB - Hess(i,j)*vec(1:NE,1:NKe,1:NQ,i)
            end do
            !$OMP end parallel do
            do i = 1, j
                tempB = tempB - Hess(i,j)*vec(1:NE,1:NKe,1:NQ,1:2,i)
            end do
            !ompfinish = omp_get_wtime()

            !print *, j, "2:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            Hess(j+1,j) = norm2(tempB)
            vec(1:NE,1:NKe,1:NQ,1:2,j+1) = tempB/Hess(j+1,j)
            !ompfinish = omp_get_wtime()

            !print *, j, "3:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            ! Append latest updated column of Hess to R, with cumulative rotations applied
            Rk(j+1,j) = Hess(j+1,j)
            Rk(1:j,j) = matmul(Qk(1:j,1:j),Hess(1:j,j)) ! Maybe take advantage of Q's structure? I think it's Hessian itself.
            !ompfinish = omp_get_wtime()

            !print *, j, "4:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
            denom = hypot(Rk(j,j),Rk(j+1,j))
            Rot(:,1) = [Rk(j,j)/denom, -Rk(j+1,j)/denom]
            Rot(:,2) = [-Rot(2,1), Rot(1,1)]
            !ompfinish = omp_get_wtime()

            !print *, j, "5:", ompfinish - ompstart, "seconds"

            !ompstart = omp_get_wtime()
            Rk(j:j+1,j) = [denom, 0.0]
            ! Update cumulative rotation matrix
            tempC(1:2,1:j+1) = Qk(j:j+1,1:j+1)
            Qk(j:j+1,1:j+1) = matmul(Rot,tempC(1:2,1:j+1))
            !ompfinish = omp_get_wtime()

            !print *, j, "6:", ompfinish - ompstart, "seconds"

            !! Check for convergence
            print *, j, abs(resid*Qk(j+1,1))
            !if (isnan(abs(resid*Qk(j+1,1)))) then
            !    print *, any(isnan(vec))
            !    print *, any(isnan(Rk))
            !    print *, any(isnan(Qk))
            !    print *, any(isnan(Hess))
            !    print *, any(isnan(tempC))
            !    print *, any(isnan(Rot))
            !    stop
            !end if
            if (abs(resid*Qk(j+1,1)) .le. l_convergence) then
                gk(1:j+1) = resid*Qk(1:j+1,1)
                do i = j, 1, -1
                    y(i) = (gk(i) - sum(Rk(i,i+1:j)*y(i+1:j)))/Rk(i,i)
                end do

                do i = 1, j ! THIS STEP MAY BE EASIER IF RESHAPES AND PACKS ARE USED
                    vector = vector + vec(1:NE,1:NKe,1:NQ,1:2,i)*y(i)
                end do

                exit GMRESm
            end if
            !if (isnan(resid*Qk(j+1,1))) stop
        end do

        !ompstart = omp_get_wtime()

        total = total + 1

        ! Do mth step differently
        tempA = vec(1:NE,1:NKe,1:NQ,1:2,m)

        call PN_FEXS_Lgg &
        (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA)

        !ompfinish = omp_get_wtime()

        !print *, "7:", ompfinish - ompstart, "seconds"

        !ompstart = omp_get_wtime()
        !$OMP parallel do shared(Hess) private(i)
        do i = 1, m
            Hess(i,m) = sum(tempA*vec(1:NE,1:NKe,1:NQ,1:2,i))
        end do
        !$OMP end parallel do
        Hess(m+1,m) = sqrt(norm2(tempA)**2 - sum(Hess(1:m,m)**2))
        !ompfinish = omp_get_wtime()

        !print *, "8:", ompfinish - ompstart, "seconds"

        !if (any(isnan(Hess))) then
        !    print *, Hess(m+1,m)
        !    print *, norm2(tempA)**2
        !    print *, sum(Hess(1:m,m)**2)
        !    print *, norm2(tempA)**2 - sum(Hess(1:m,m)**2)
        !    stop
        !end if

        ! Append latest updated column of Hess to R, with cumulative rotations applied
        Rk(m+1,m) = Hess(m+1,m)
        Rk(1:m,m) = matmul(Qk(1:m,1:m),Hess(1:m,m)) ! Maybe take advantage of Q's structure? I think its Hessian itself.

        ! Rotate e_j and e_j+1 in this column so that R(j+1,j) is brought to zero.
        denom = hypot(Rk(m,m),Rk(m+1,m))
        Rot(:,1) = [Rk(m,m)/denom, -Rk(m+1,m)/denom] ! Should I transpose these? were they supposed to be transposed? (~3/31/23)
        Rot(:,2) = [-Rot(2,1), Rot(1,1)] !

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
            vector = vector + vec(1:NE,1:NKe,1:NQ,1:2,i)*y(i)
        end do

        tempA = vector
        call PN_FEXS_Lgg &
        (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, tempA)

        resid0 = source - tempA

        !! Slower??
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

        if (abs(gk(m+1)) <= convergence) then
            exit GMRESm
        end if

        if (total >= max_iter) then
            gk(1:m+1) = resid*Qk(1:m+1,1)
            do i = m, 1, -1
                y(i) = (gk(i) - sum(Rk(i,i+1:m)*y(i+1:m)))/Rk(i,i)
            end do

            do i = 1, m
                vector = vector + vec(1:NE,1:NKe,1:NQ,1:2,i)*y(i)
            end do
            print *, "GMRESm: Did not converge within", max_iter, "iterations."
            print *, "Solution being given with:", abs(gk(j+1)), "residual."
            !print *, "iterations:", total
            !print *, "residual 2:", abs(gk(j+1))
            exit GMRESm
            !stop
        end if
    end do GMRESm
end subroutine PN_FEXS_GMRESm

subroutine PN_FEXS_K &
    (l_NL, Sigma, scatterers, vector)
    implicit none
    integer, intent(in) :: l_NL
    real, dimension(:,:,:,:,:), intent(in) :: Sigma ! (np,n,gpr,l,mat)
    real, dimension(:,:,:), intent(in) :: scatterers ! (k,q,gpr)

    real, dimension(:,:,:), intent(inout) :: vector ! (k,q,n) ! DOESN'T SET IT TO ZERO. MAKE CONSISTEN W/ SN?

    integer :: gpr, k, l, mat, n, np
    integer :: group

    group = size(Sigma,3)


    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(2) shared(vector) private(gpr, l, n, np)
        do n = 1, 2
            do l = 0, l_NL
                do gpr = 1, group
                    do np = 1, 2
                        vector(1:NK,l*(l+1)-l+1:l*(l+1)+l+1,n) = vector(1:NK,l*(l+1)-l+1:l*(l+1)+l+1,n) + &
                            (fourpi*Sigma(np,n,gpr,l+1,1)/(2*l+1))*&
                            scatterers(1:NK,l*(l+1)-l+1:l*(l+1)+l+1,gpr+np-1)
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        do n = 1, 2
            do l = 0, l_NL
                do gpr = 1, group
                    do np = 1, 2
                        do mat = 1, Nmats
                            do k = 1, NKmat(mat)
                                vector(nodesinmat(k,mat),l*(l+1)-l+1:l*(l+1)+l+1,n) = &
                                    vector(nodesinmat(k,mat),l*(l+1)-l+1:l*(l+1)+l+1,n) + &
                                    Nelnodes(nodesinmat(k,mat),mat)*(fourpi*Sigma(np,n,gpr,l+1,mat)/(2*l+1))*&
                                    scatterers(nodesinmat(k,mat),l*(l+1)-l+1:l*(l+1)+l+1,gpr+np-1)/Nglobal(nodesinmat(k,mat))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
end subroutine PN_FEXS_K

subroutine PN_FEXS_Lgg &
    (l_NL, eprime, Ceff, sgn, Awl, Aup, Adn, I1invI2, I1invI2f, I1invI3vector, Sigma, Sigmat, vector)
    implicit none
    ! USES METHOD 1 WHICH SCALES BETTER WITH NE BUT NOT WELL AT ALL WITH NQ
    integer, intent(in) :: l_NL
    integer, dimension(:,:), intent(in) :: eprime ! (f,e)
    integer, dimension(:,:), intent(in) :: Ceff
    integer, dimension(:,:), intent(in) :: sgn
    real, dimension(:,:,:), intent(in) :: Awl ! (q,q,dir) for Method 1, UNDECIDED for Method 2
    real, dimension(:,:,:), intent(in) :: Aup ! (f,q,q)
    real, dimension(:,:,:), intent(in) :: Adn ! (f,q,q)
    !real, dimension(:,:,:,:), intent(in) :: Aup ! (e,q,f,q) for Method 1, (q,q,f,e) for Method 2
    !real, dimension(:,:,:,:), intent(in) :: Adn ! (e,q,f,q) for Method 1, (q,q,f,e) for Method 2
    real, dimension(:,:,:,:), intent(in) :: I1invI2 ! (e,kp,f,k) for Methods 1 and 2
    real, dimension(:,:,:,:), intent(in) :: I1invI2f ! (e,kp,f,k) for Methods 1 and 2
    real, dimension(:,:,:,:), intent(in) :: I1invI3vector ! (e,kp,dir,k) for Method 1, (dir,e,k,kp) for Method 2
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (np,n,l+1,mat)
    real, dimension(:,:,:), intent(in) :: Sigmat ! (np,n,mat)

    real, dimension(:,:,:,:), intent(inout) :: vector ! (e,k,q,n)

    integer :: dir, e, f, k, kp, l, mat, n, np, q, qp
    integer :: NQ
    !integer, dimension(:), allocatable :: ell
    real, dimension(:,:,:,:), allocatable :: Mterm
    real, dimension(:,:,:,:,:), allocatable :: Amat
    !real :: ompstart, ompfinish
    !real :: sumtest

    NQ = (l_NL+1)**2
    allocate(Mterm(NE,NKe,NQ,2))
    !allocate(ell(NQ))

    !do q = 1, NQ
    !    ell(q) = ceiling(sqrt(real(q)))-1
    !end do

    !if (any(isnan(vector))) then
    !    print *, "VECTOR"
    !    stop
    !end if

    ! MTERM - Not too bad
    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(2) shared(Mterm) private(l, n, np)
        do n = 1, 2
            do l = 0, l_NL
                do np = 1, 2
                    Mterm(1:NE,1:NKe,l*(l+1)-l+1:l*(l+1)+l+1,n) = Mterm(1:NE,1:NKe,l*(l+1)-l+1:l*(l+1)+l+1,n) + &
                        (Sigmat(np,n,1) - fourpi*Sigma(np,n,l+1,1)/(2*l+1))*&
                        vector(1:NE,1:NKe,l*(l+1)-l+1:l*(l+1)+l+1,np)
                end do
            end do
        end do
        !$OMP end parallel do
    else
        do mat = 1, Nmats
            !!$OMP parallel do collapse(2) shared(Mterm) private(l, mat, n, np)
            do n = 1, 2
                do l = 0, l_NL
                    do np = 1, 2
                        Mterm(elsinmat(1:NEmat(mat),mat),1:NKe,l*(l+1)-l+1:l*(l+1)+l+1,n) = &
                            Mterm(elsinmat(1:NEmat(mat),mat),1:NKe,l*(l+1)-l+1:l*(l+1)+l+1,n) + &
                            (Sigmat(np,n,mat) - fourpi*Sigma(np,n,l+1,mat)/(2*l+1))*&
                            vector(elsinmat(1:NEmat(mat),mat),1:NKe,l*(l+1)-l+1:l*(l+1)+l+1,np)
                    end do
                end do
            end do
            !!$OMP end parallel do
        end do
    end if

    !if (any(isnan(Mterm))) then
    !    print *, "MTERM"
    !    stop
    !end if

    ! GTERM - BAD
    !ompstart = omp_get_wtime()
    allocate(Amat(NE,NQ,NKe,3,2))
    !$OMP parallel do collapse(3) shared(Amat) private(dir, k, n)
    do n = 1, 2
        do dir = 1, 3
            do k = 1, NKe
                Amat(1:NE,1:NQ,k,dir,n) = matmul(vector(1:NE,k,1:NQ,n),Awl(1:NQ,1:NQ,dir))
            end do
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(3) shared(Mterm) private(dir, k, n, q)
    do n = 1, 2
        do q = 1, NQ
            do k = 1, NKe
                do dir = 1, 3
                    Mterm(1:NE,k,q,n) = Mterm(1:NE,k,q,n) - &
                        sum(I1invI3vector(1:NE,1:NKe,dir,k)*Amat(1:NE,q,1:NKe,dir,n),2)
                end do
            end do
        end do
    end do
    !$OMP end parallel do
    deallocate(Amat)

    !if (any(isnan(Mterm))) then
    !    print *, "GTERM"
    !    stop
    !end if
    !!!!!!$OMP parallel do collapse(3) shared(Mterm) private(k,e,dir) ! Use for large NQ?
    !!!!!do k = 1, NKe
    !!!!!    do e = 1, NE
    !!!!!        do dir = 1, 3
    !!!!!            Mterm(e,k,1:NQ) = Mterm(e,k,1:NQ) - &
    !!!!!            matmul(matmul(Awl(1:NQ,1:NQ,dir),transpose(vector(e,1:NKe,1:NQ))),I1invI3vector(e,k,1:NKe,dir))
    !!!!!        end do
    !!!!!    end do
    !!!!!end do
    !!!!!!$OMP end parallel do
    !!!!
    !!!!!!$OMP parallel do collapse(4) shared(Mterm) private(q, dir, kp, k)
    !!!!!do dir = 1, 3
    !!!!!    do q = 1, NQ
    !!!!!        do kp = 1, NKe
    !!!!!            do k = 1, NKe
    !!!!!                Mterm(1:NE,k,q) = Mterm(1:NE,k,q) - I1invI3vector(1:NE,k,kp,dir)*&
    !!!!!                    matmul(vector(1:NE,kp,Ja(Ia(q,dir):Ia(q+1,dir)-1,dir)),Awl(Ia(q,dir):Ia(q+1,dir)-1,dir))
    !!!!!            end do
    !!!!!        end do
    !!!!!    end do
    !!!!!end do
    !!!!!!$OMP end parallel do

    ! FTERM - BAD
    allocate(Amat(NE,NQ,NFe,NKe,2))
    Amat = 0
    !$OMP parallel do collapse(4) shared(Amat) private(k, f, q, kp, n)
    do n = 1, 2
        do k = 1, NKe
            do f = 1, NFe
                do q = 1, NQ
                    do kp = 1, NKe
                        Amat(1:NE,q,f,k,n) = Amat(1:NE,q,f,k,n) + &
                            I1invI2(1:NE,kp,f,k)*vector(1:NE,kp,q,n)
                    end do
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(3) shared(Mterm) private(q, k, f, qp, n) ! Does this need to be collapse(2)?? It seems not, but I don't know why not.
    do n = 1, 2
        do q = 1, NQ
            do k = 1, NKe ! LOOPING OVER e MIGHT MAKE THIS BETTER FOR PARALLELIZATION ON SUPERCOMPUTER
                do f = 1, NFe
                    Mterm(1:NE,k,q,n) = Mterm(1:NE,k,q,n) + &
                        sgn(1:NE,f)*sum(Amat(1:NE,1:NQ,f,k,n)*Aup(Ceff(1:NE,f),1:NQ,q),2) ! sum(Amat(1:NE,1:NQ,f,k,n)*Aup(1:NE,1:NQ,f,q),2)
                end do
            end do
        end do
    end do
    !$OMP end parallel do
    !if (any(isnan(Mterm))) then
    !    print *, "FTERM"
    !    stop
    !end if
    ! SWEPT FTERM - BAD
    Amat = 0 ! Because eprime(f,e) may be zero, so it will cycle out of some [f,e] pairs that will retain their values
    !$OMP parallel do collapse(4) shared(Amat) private(k, f, e, n)
    do n = 1, 2
        do k = 1, NKe
            do f = 1, NFe
                do e = 1, NE
                    if (eprime(f,e) .eq. 0) cycle
                    Amat(e,1:NQ,f,k,n) = &
                        matmul(transpose(vector(eprime(f,e),1:NKe,1:NQ,n)),I1invI2f(e,1:NKe,f,k)) ! USE A REORDERED eprime AND DO SBR STUFF AGAIN? as in, reorganize eprime(e,f) such that it's lower tri/more easily accesible? (~3/10/23)
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    !$OMP parallel do collapse(3) shared(Mterm) private(q, k, f, n)
    do n = 1, 2
        do q = 1, NQ
            do k = 1, NKe
                do f = 1, NFe
                    Mterm(1:NE,k,q,n) = Mterm(1:NE,k,q,n) + &
                        sgn(1:NE,f)*sum(Amat(1:NE,1:NQ,f,k,n)*Aup(Ceff(1:NE,f),1:NQ,q),2) ! sum(Amat(1:NE,1:NQ,f,k,n)*Adn(1:NE,1:NQ,f,q),2)
                end do
            end do
        end do
    end do
    !$OMP end parallel do
    !if (any(isnan(Mterm))) then
    !    print *, "SWEPT FTERM"
    !    stop
    !end if
    vector = Mterm
    ! CONSIDER PROBLEMS POSED BY ARRAY TEMPORARIES

    !stop
end subroutine PN_FEXS_Lgg

! MISC
end module iteration