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

module t1_ele_iteration
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

contains
subroutine t1_ele_SN_write_field &
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
end subroutine t1_ele_SN_write_field

subroutine t1_ele_SN_MGXS_EI_SI &
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
    real, dimension(:,:,:,:), allocatable, intent(in) :: source
    real, dimension(:), intent(in) :: fEg
    real, dimension(:,:,:,:), allocatable, intent(in) :: bdysrc

    real, dimension(:,:,:,:), allocatable, intent(inout) :: psi
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

    allocate(psi(NK,Nmu,Nphi,l_G))
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
                do k = 1, NKe
                    do e = 1, NE
                        t_source(e,k,:,:) = source(Cekk(e,k),:,:,1)
                    end do
                end do
            else
                t_source = 0
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(1)*bdysrc
            end if
        else
            allocate(t_Sigma, source = Sigma(1:g-1,g,:,:))
            call t1_ele_SN_MGXS_K &
                (Cekk, wL, cosmat, qPa, qfactorialmat, t_Sigma, psi, t_source)
            if (allocated(source)) then
                do k = 1, NKe
                    do e = 1, NE
                        t_source(e,k,:,:) = t_source(e,k,:,:) + source(Cekk(e,k),:,:,g)
                    end do
                end do
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
        call t1_ele_SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()

        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        call t1_ele_SN_MGXS_SI &
            (convergence, pmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, &
            esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, t_Sigmagg, t_source, t_bdysrc, Tprimeinv, total, t_phiu, t_psi)

        do k = 1, NKe
            do e = 1, NE
                psi(Cekk(e,k),:,:,g) = psi(Cekk(e,k),:,:,g) + &
                    t_psi(e,k,:,:)/Nglobal(Cekk(e,k))
            end do
        end do

        phi(:,:,g) = t_phiu

        ompfinish = omp_get_wtime()

        print *, "CONVERGED WITH:", total, "ITERATES"
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    phi = phi - phiu

    if (particle .eq. 1) then
        call t1_ele_SN_write_field("photon", phi, psi)
    else if (particle .eq. 2) then
        call t1_ele_SN_write_field("electron", phi, psi)
    end if
end subroutine t1_ele_SN_MGXS_EI_SI

subroutine t1_ele_SN_MGXS_construct_inverted_transport_matrix &
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
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine t1_ele_SN_MGXS_construct_inverted_transport_matrix

subroutine t1_ele_SN_MGXS_SI &
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
    real, dimension(:,:), intent(inout) :: phi ! (e,k)
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

    call t1_ele_SN_MGXS_sweep &
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

        call t1_ele_SN_MGXS_Kgg &
            (wL, cosmat, qPa, qfactorialmat, Sigma, source)

        call t1_ele_SN_MGXS_sweep &
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
end subroutine t1_ele_SN_MGXS_SI

subroutine t1_ele_SN_MGXS_K(Cekk, wL, cosmat, qPa, qfactorialmat, Sigma, scatterers, vector)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat ! (j,jp,m+1)
    real, dimension(:,:), intent(in) :: qPa ! (i,q)
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:), intent(in) :: Sigma ! (gp,l+1)
    real, dimension(:,:,:,:), intent(in) :: scatterers ! (k,i,j,g)

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

    group = size(Sigma,1)
    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    allocate(Amat(NE,NKe,Nmu,Nphi,NL+1))
    if (Nmats .eq. 1) then ! THERE MAY BE PARALLELIZATION PROBLEMS HERE!!
        !$OMP parallel do collapse(4) shared(Amat) private(l, jp, ip, k)
        do l = 0, NL
            do jp = 1, Nphi
                do ip = 1, Nmu
                    do k = 1, NKe
                        Amat(1:NE,k,ip,jp,l+1) = twopi*matmul(scatterers(Cekk(1:NE,k),ip,jp,1:group), &
                            Sigma(1:group,l+1,1))/Nphi
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
                        Amat(e,k,ip,1:Nphi,l+1) = twopi*matmul(scatterers(Cekk(e,k),ip,1:Nphi,1:group),&
                            Sigma(1:group,l+1,eltomat(e)))/Nphi
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    end if

    allocate(Bmat(NE,NKe,Nphi,NQ))
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

    if (size(Sigma,2) .eq. NL+2) then ! THERE MAY BE PARALLELIZATION PROBLEMS HERE!!
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(3) shared(vector) private(i, j, k)
            do j = 1, Nphi
                do i = 1, Nmu
                    do k = 1, NKe
                        vector(1:NE,k,i,j) = matmul(scatterers(Cekk(1:NE,k),i,j,1:group),&
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
                        vector(e,k,i,1:Nphi) = matmul(scatterers(Cekk(e,k),i,1:Nphi,1:group),&
                            Sigma(1:group,NL+2,eltomat(e))) ! IDEA: Use mattoel(:,mat), which lists the set of mats for an element? But sizes are different. figure out. ---> do mattoel(1:Nelsinmat(mat),mat).
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
end subroutine t1_ele_SN_MGXS_K

subroutine t1_ele_SN_MGXS_Kgg(wL, cosmat, qPa, qfactorialmat, Sigma, vector)
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
end subroutine t1_ele_SN_MGXS_Kgg

subroutine t1_ele_SN_MGXS_sweep &
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
            ! vector(1:NE,1:NKe,i,j) = transpose(l_vector(1:NKe,l_enew)) !!! TRY THIS
        end do
    end do
    !$OMP end parallel do
end subroutine t1_ele_SN_MGXS_sweep

subroutine t1_ele_SN_MGXS_EI_GMRESm &
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

    allocate(psi(NK,Nmu,Nphi,l_G))
    allocate(phi(NE,NKe,l_G))
    allocate(t_source(NE,NKe,Nmu,Nphi))
    if (allocated(bdysrc)) allocate(t_bdysrc(NKe,size(bdysrc,2),Nmu,Nphi))
    allocate(t_Sigmagg(size(Sigma,3),size(Sigma,4)))
    allocate(t_Sigmat(Nmats))
    allocate(t_psi(NE,NKe,Nmu,Nphi))
    allocate(Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

    psi = 0
    phi = 0
    do g = 1, l_G
        total = 0
        print *, "ENERGY GROUP:", g

        ompstart = omp_get_wtime()
        t_Sigmagg = Sigma(g,g,:,:)
        t_Sigmat = Sigmat(g,:)

        call t1_ele_SN_MGXS_construct_inverted_transport_matrix &
            (khat, I1invI2, I1invI3vector, enew, varsigmaup, t_Sigmat, Tprimeinv)
        ompfinish = omp_get_wtime()
        print *, "TRANSPORT MATRIX INVERTED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        if (g .eq. 1) then
            if (allocated(source)) then
                do k = 1, NKe
                    do e = 1, NE
                        t_source(e,k,:,:) = source(Cekk(e,k),:,:,1)
                    end do
                end do
            else
                t_source = 0
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(1)*bdysrc
            end if
        else
            allocate(t_Sigma, source = Sigma(1:g-1,g,:,:))
            call t1_ele_SN_MGXS_K &
                (Cekk, wL, cosmat, qPa, qfactorialmat, t_Sigma, psi, t_source)
            if (allocated(source)) then
                do k = 1, NKe
                    do e = 1, NE
                        t_source(e,k,:,:) = t_source(e,k,:,:) + source(Cekk(e,k),:,:,g)
                    end do
                end do
            end if
            if (allocated(bdysrc)) then
                t_bdysrc = fEg(g)*bdysrc
            end if
            deallocate(t_Sigma)
        end if

        call t1_ele_SN_MGXS_sweep &
            (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
            varsigmadown, I1invI2f, t_bdysrc, Tprimeinv, t_source)

        ompfinish = omp_get_wtime()
        print *, "SOURCE UPDATED. TIME (s):", ompfinish - ompstart

        ompstart = omp_get_wtime()
        !t_psi = t_source/t_Sigmat(:,1) ! GUESS
        t_psi = t_source
        call t1_ele_SN_MGXS_GMRESm &
            (convergence, pmax, mmax, Cekk, bdyel, wL, cosmat, qPa, qfactorialmat, &
            esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
            I1invI2f, t_Sigmagg, t_source, t_bdysrc, Tprimeinv, total, t_psi)

        do k = 1, NKe
            do e = 1, NE
                psi(Cekk(e,k),:,:,g) = psi(Cekk(e,k),:,:,g) + t_psi(e,k,:,:)/Nglobal(Cekk(e,k))
            end do
        end do

        do k = 1, NKe
            do j = 1, Nphi
                phi(1:NE,k,g) = phi(1:NE,k,g) + twopi*matmul(t_psi(1:NE,k,1:Nmu,j),wL)/Nphi
            end do
        end do

        print *, "CONVERGED WITH:", total, "ITERATES"

        ompfinish = omp_get_wtime()
        print *, "TOTAL TIME (s):", ompfinish - ompstart
    end do

    if (particle .eq. 1) then
        call t1_ele_SN_write_field("photon", phi, psi)
    else if (particle .eq. 2) then
        call t1_ele_SN_write_field("electron", phi, psi)
    end if
end subroutine t1_ele_SN_MGXS_EI_GMRESm

subroutine t1_ele_SN_MGXS_GMRESm &
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

    call t1_ele_SN_MGXS_swept_Lgg &
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

            call t1_ele_SN_MGXS_swept_Lgg &
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

        call t1_ele_SN_MGXS_swept_Lgg &
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

        call t1_ele_SN_MGXS_swept_Lgg &
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
end subroutine t1_ele_SN_MGXS_GMRESm

subroutine t1_ele_SN_MGXS_swept_Lgg &
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

    call t1_ele_SN_MGXS_Kgg(wL, cosmat, qPa, qfactorialmat, Sigma, tempA)

    call t1_ele_SN_MGXS_sweep &
        (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
        varsigmadown, I1invI2f, nobdysrc, Tprimeinv, tempA)

    vector = vector - tempA
end subroutine t1_ele_SN_MGXS_swept_Lgg

end module t1_ele_iteration