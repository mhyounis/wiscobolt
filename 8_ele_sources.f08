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

module ele_sources
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

contains
!! UNOPTIMIZED
subroutine t2_ele_SN_MGXS_construct_uncollided_source &
    (Cekk, Pa, factorialmat, Sigma, Pmlk, cmjk, phiu, source)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:), intent(in) :: Pa
    real, dimension(:,:), intent(in) :: factorialmat
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (gp,g,l,mat)
    real, dimension(:,:,:), intent(in) :: Pmlk
    real, dimension(:,:,:), intent(in) :: cmjk
    real, dimension(:,:,:), intent(in) :: phiu ! (e,k,g)

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: source ! (e,k,i,j,g)

    integer :: e, g, i, j, k, kg, l, m, mat, n, np, q
    integer :: NQ
    integer :: l_G
    integer :: l_NE
    integer :: l_k
    integer, dimension(:), allocatable :: ell
    real, dimension(:,:,:,:), allocatable :: Amat ! (e,k,g,l+1)
    real, dimension(:,:,:), allocatable :: tempBmat ! (e,k,g)
    real, dimension(:,:,:,:,:,:), allocatable :: Cmat ! (e,k,i,g,l+1,m+1)
    real, dimension(:,:,:,:,:), allocatable :: Dmat ! (e,k,i,j,g)
    integer, dimension(:), allocatable :: invnodesinbeam

    allocate(invnodesinbeam(NK))
    invnodesinbeam = 0
    do k = 1, size(nodesinbeam)
        invnodesinbeam(nodesinbeam(k)) = k
    end do

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    l_G = size(Sigma,1)
    l_NE = size(elsinbeam)

    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]

    allocate(Amat(l_NE,NKe,l_G,NL+1))
    if (Nmats .eq. 1) then
        !$OMP parallel do collapse(2) shared(Amat) private(l, k)
        do l = 0, NL
            do k = 1, NKe
                Amat(1:l_NE,k,1:l_G,l+1) = &
                    matmul(phiu(elsinbeam,k,1:l_G),Sigma(1:l_G,1:l_G,l+1,1))
            end do
        end do
        !$OMP end parallel do
    else
        Amat = 0
        !$OMP parallel do collapse(3) shared(Amat) private(e, k, l)
        do l = 0, NL
            do k = 1, NKe
                do e = 1, l_NE
                    Amat(e,k,1:l_G,l+1) = &
                        matmul(transpose(Sigma(1:l_G,1:l_G,l+1,eltomat(elsinbeam(e)))),&
                        phiu(elsinbeam(e),k,1:l_G))
                end do
            end do
        end do
        !$OMP end parallel do
    end if

    allocate(tempBmat(l_NE,NKe,l_G))
    allocate(Cmat(l_NE,NKe,Nmu,l_G,NL+1,NL+1))
    !$OMP parallel do shared(Cmat) private(tempBmat) private(m, l, i) ! CONSIDER USING q SO I CAN PARALLELIZE! ALSO MORE CONVENIENT
    do m = 0, NL
        do l = m, NL
            tempBmat = factorialmat(l+1,m+1)*Amat(:,:,:,l+1)
            do i = 1, Nmu
                Cmat(1:l_NE,1:NKe,i,1:l_G,l+1,m+1) = Pa(i,l+1,m+1)*tempBmat
            end do
        end do
    end do
    !$OMP end parallel do
    deallocate(Amat)
    deallocate(tempBmat)

    allocate(Dmat(l_NE,NKe,Nmu,Nphi,l_G))
    if (beam_angular_dist .eq. "spherical") then
        !$OMP parallel do collapse(4) shared(Dmat) private(e, g, i, k, kg)
        do g = 1, l_G
            do i = 1, Nmu
                do k = 1, NKe
                    do e = 1, l_NE
                        kg = invnodesinbeam(Cekk(elsinbeam(e),k))
                        Dmat(e,k,i,1:Nphi,g) = matmul(cmjk(1:Nphi,1:NL+1,kg),&
                            sum(Pmlk(1:NL+1,1:NL+1,kg)*&
                            Cmat(e,k,i,g,1:NL+1,1:NL+1),1))
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    else
        !$OMP parallel do collapse(4) shared(Dmat) private(e, g, i, k)
        do g = 1, l_G
            do i = 1, Nmu
                do k = 1, NKe
                    do e = 1, l_NE
                        Dmat(e,k,i,1:Nphi,g) = matmul(cmjk(1:Nphi,1:NL+1,1),&
                            sum(Pmlk(1:NL+1,1:NL+1,1)*&
                            Cmat(e,k,i,g,1:NL+1,1:NL+1),1))
                    end do
                end do
            end do
        end do
        !$OMP end parallel do
    end if
    deallocate(Cmat)

    allocate(source(NE,NKe,Nmu,Nphi,l_G))
    !$OMP parallel do shared(source) private(e)
    do e = 1, l_NE
        source(elsinbeam(e),:,:,:,:) = Dmat(e,:,:,:,:) ! This solution may suck. Thoroughly consider fixing this and the next loop. ALSO TRY ADDING ARRAY DIMENSIONS
    end do
    !$OMP end parallel do
end subroutine t2_ele_SN_MGXS_construct_uncollided_source

subroutine ele_SN_to_SN_MGXS_particle_particle_scattering_source &
    (Cekk, wL, cosmat, Pa, factorialmat, qPa, qfactorialmat, Sigma, Pmlk, cmjk, phiu, psi, source)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:), intent(in) :: Pa
    real, dimension(:,:), intent(in) :: factorialmat
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:), intent(in) :: Pmlk
    real, dimension(:,:,:), intent(in) :: cmjk
    real, dimension(:,:,:), allocatable, intent(in) :: phiu ! (e,k,g)
    real, dimension(:,:,:,:,:), intent(in) :: psi ! (e,k,i,j,g)

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: source

    integer :: e, g, gpr, i, ip, jp, k, kg, l, m, mat, q
    integer :: NQ
    integer :: G1
    integer :: G2
    integer :: l_NE
    integer :: l_k
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:), allocatable :: Amat ! (nodes,g,l+1)
    real, dimension(:,:,:), allocatable :: tempBmat ! (nodes,g)
    real, dimension(:,:,:,:,:,:), allocatable :: Cmat ! (nodes,i,g,l+1,m+1)
    real, dimension(:,:,:,:,:), allocatable :: Dmat ! (nodes,i,j,g)
    real, dimension(:,:,:,:,:), allocatable :: g_Amat ! (e,k,i,j,l)
    real, dimension(:,:,:,:), allocatable :: g_Bmat ! (e,k,jp,q)
    real, dimension(:,:,:,:), allocatable :: g_Cmat ! (e,k,j,q)
    integer, dimension(:), allocatable :: invnodesinbeam

    allocate(invnodesinbeam(NK))
    invnodesinbeam = 0
    do k = 1, size(nodesinbeam)
        invnodesinbeam(nodesinbeam(k)) = k
    end do

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    G1 = size(Sigma,2) ! To
    G2 = size(Sigma,1) ! From
    if (G2 .ne. size(psi,5)) then
        print *, "ARRAY SIZE DOESNT MATCH"
        stop
    end if
    l_NE = size(elsinbeam)

    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    if (allocated(phiu)) then
        allocate(Amat(l_NE,NKe,G1,NL+1))
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(2) shared(Amat) private(l, k)
            do l = 0, NL
                do k = 1, NKe
                    Amat(1:l_NE,k,1:G1,l+1) = &
                    matmul(phiu(elsinbeam,k,1:G2),Sigma(1:G2,1:G1,l+1,1))
                end do
            end do
            !$OMP end parallel do
        else
            Amat = 0
            !$OMP parallel do collapse(3) shared(Amat) private(e, k, l)
            do l = 0, NL
                do k = 1, NKe
                    do e = 1, l_NE
                        Amat(e,k,1:G1,l+1) = &
                            matmul(transpose(Sigma(1:G2,1:G1,l+1,eltomat(elsinbeam(e)))),&
                            phiu(elsinbeam(e),k,1:G2))
                    end do
                end do
            end do
            !$OMP end parallel do
        end if

        allocate(tempBmat(l_NE,NKe,G1))
        allocate(Cmat(l_NE,NKe,Nmu,G1,NL+1,NL+1))
        !$OMP parallel do shared(Cmat) private(tempBmat) private(m, l, i) ! CONSIDER USING q SO I CAN PARALLELIZE! ALSO MORE CONVENIENT
        do m = 0, NL
            do l = m, NL
                tempBmat = factorialmat(l+1,m+1)*Amat(:,:,:,l+1)
                do i = 1, Nmu
                    Cmat(1:l_NE,1:NKe,i,1:G1,l+1,m+1) = Pa(i,l+1,m+1)*tempBmat
                end do
            end do
        end do
        !$OMP end parallel do
        deallocate(Amat)
        deallocate(tempBmat)

        allocate(Dmat(l_NE,NKe,Nmu,Nphi,G1))
        if (beam_angular_dist .eq. "spherical") then
            !$OMP parallel do collapse(4) shared(Dmat) private(e, g, i, k, kg)
            do g = 1, G1
                do i = 1, Nmu
                    do k = 1, NKe
                        do e = 1, l_NE
                            kg = invnodesinbeam(Cekk(elsinbeam(e),k))
                            Dmat(e,k,i,:,g) = matmul(cmjk(1:Nphi,1:NL+1,kg),&
                                sum(Pmlk(1:NL+1,1:NL+1,kg)*&
                                Cmat(e,k,i,g,1:NL+1,1:NL+1),1))
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        else
            !$OMP parallel do collapse(4) shared(Dmat) private(e, g, i, k)
            do g = 1, G1
                do i = 1, Nmu
                    do k = 1, NKe
                        do e = 1, l_NE
                            Dmat(e,k,i,:,g) = matmul(cmjk(1:Nphi,1:NL+1,1),&
                                sum(Pmlk(1:NL+1,1:NL+1,1)*&
                                Cmat(e,k,i,g,1:NL+1,1:NL+1),1))
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        end if
        deallocate(Cmat)

        !$OMP parallel do shared(source) private(e)
        do e = 1, l_NE
            source(elsinbeam(e),:,:,:,:) = Dmat(e,:,:,:,:)
        end do
        !$OMP end parallel do

        deallocate(Dmat)
    else
        source = 0
    end if

    allocate(g_Amat(NE,NKe,Nmu,Nphi,NL+1))
    allocate(g_Bmat(NE,NKe,Nphi,NQ))
    allocate(g_Cmat(NE,NKe,Nphi,NQ))

    do g = 1, G1
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(4) shared(g_Amat) private(l, jp, ip, k)
            do l = 0, NL
                do jp = 1, Nphi
                    do ip = 1, Nmu
                        do k = 1, NKe
                            g_Amat(1:NE,k,ip,jp,l+1) = twopi*matmul(psi(1:NE,k,ip,jp,1:G2),&
                                Sigma(1:G2,g,l+1,1))/Nphi
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        else
            !$OMP parallel do collapse(4) shared(g_Amat) private(e, l, ip, k)
            do l = 0, NL
                do ip = 1, Nmu
                    do k = 1, NKe
                        do e = 1, NE
                            g_Amat(e,k,ip,1:Nphi,l+1) = twopi*matmul(psi(e,k,ip,1:Nphi,1:G2),&
                                Sigma(1:G2,g,l+1,eltomat(e)))/Nphi
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        end if

        !$OMP parallel do shared(g_Bmat) private(k, jp, q)
        do q = 1, NQ
            do jp = 1, Nphi
                do k = 1, NKe
                    g_Bmat(1:NE,k,jp,q) = matmul(g_Amat(1:NE,k,1:Nmu,jp,ell(q)+1),&
                        qPa(1:Nmu,q)*wL(1:Nmu))
                end do
            end do
            g_Bmat(1:NE,1:NKe,1:Nphi,q) = qfactorialmat(q)*g_Bmat(1:NE,1:NKe,1:Nphi,q)
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) shared(g_Cmat) private(k, q)
        do q = 1, NQ
            do k = 1, NKe
                g_Cmat(1:NE,k,1:Nphi,q) = matmul(g_Bmat(1:NE,k,1:Nphi,q),&
                    cosmat(1:Nphi,1:Nphi,emm(q)+1))
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do shared(source) private(i, q)
        do i = 1, Nmu
            do q = 1, NQ
                source(1:NE,1:NKe,i,1:Nphi,g) = source(1:NE,1:NKe,i,1:Nphi,g) + &
                    g_Cmat(1:NE,1:NKe,1:Nphi,q)*qPa(i,q)
            end do
        end do
        !$OMP end parallel do
    end do
end subroutine ele_SN_to_SN_MGXS_particle_particle_scattering_source

! wiscoslab stuff
subroutine wiscoslab_MGXS_uncollided_fluence &
    (label, rglobal, Ep, fEspec, fE, Sigmat, phiu)
    implicit none
    character(*), intent(in) :: label
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE
    real, dimension(:,:), intent(in) :: Sigmat

    real, dimension(:,:,:), allocatable, intent(inout) :: phiu

    integer :: e, g, k, mat
    integer :: l_G
    real :: L
    real, dimension(:,:), allocatable :: ells
    real, dimension(:), allocatable :: fEg

    l_G = size(Ep) - 1

    allocate(phiu(NE,NKe,l_G))
    allocate(fEg(l_G))
    allocate(ells(NK,Nmats))

    L = rglobal(1,1)

    ells = 0
    do k = 2, NK
        ! Visit a node
        ! The path length up to this node is initializaed as the path length to the previous node
        ! The path length for the material of the previous element is incremented by the
        ! length of the previous element
        ells(k,:) = ells(k-1,:)
        ells(k,eltomat(k-1)) = ells(k,eltomat(k-1)) + rglobal(1,k-1) - rglobal(1,k)
    end do

    phiu = 0
    if (beam_energy_dist .eq. "polychromatic") then
        do g = 1, l_G
            call trapezoidal_integration_1D &
                (0, fE, fEspec, Ep(g+1), Ep(g), fEg(g))
        end do

        do g = 1, l_G
            phiu(1:NE,1:NKe,g) = fEg(g)
        end do

        do g = 1, l_G
            do k = 1, NKe
                do e = 1, NE
                    mat = eltomat(e)
                    phiu(e,k,g) = phiu(e,k,g)*&
                        exp(-Sigmat(g,mat)*ells(e+k-1,mat)) ! Units: 1/(inc. particle/cm^2)
                end do
            end do
        end do
    else if (beam_energy_dist .eq. "boxcar") then
        do g = 1, boxg
            phiu(1:NE,1:NKe,g) = (Ep(g)-Ep(g+1))/deltaE
        end do

        do g = 1, boxg
            do k = 1, NKe
                do e = 1, NE
                    mat = eltomat(e)
                    phiu(e,k,g) = phiu(e,k,g)*&
                        exp(-Sigmat(g,mat)*ells(e+k-1,mat)) ! Units: 1/(inc. particle/cm^2)
                end do
            end do
        end do
    else
        phiu(1:NE,1:NKe,1) = 1

        do k = 1, NKe
            do e = 1, NE
                mat = eltomat(e)
                phiu(e,k,1) = phiu(e,k,1)*&
                    exp(-Sigmat(1,mat)*ells(e+k-1,mat)) ! Units: 1/(inc. particle/cm^2)
            end do
        end do
    end if

    open(1, file = trim(adjustl(output_fname))//"/"//label//"_uncollided_fluence.dat", form = "unformatted")
    write(1) phiu
    close(1)

    allocate(nodesinbeam, source = [(k,k=1,NK)])
    allocate(elsinbeam, source = [(e,e=1,NE)])
end subroutine wiscoslab_MGXS_uncollided_fluence

subroutine wiscoslab_MGXS_beam_quadrature(rglobal, dz, St, N, phiu)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:), intent(in) :: dz
    real, dimension(:,:), intent(in) :: St
    integer, intent(in) :: N

    real, dimension(:,:,:), intent(inout) :: phiu

    integer :: e, mat, k, kp, p
    real :: L
    integer :: G
    real, dimension(:,:), allocatable :: ells
    real, dimension(:,:), allocatable :: I1inv
    real, dimension(:,:), allocatable :: T
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:), allocatable :: rquad
    real, dimension(:,:), allocatable :: phi
    real, dimension(:,:), allocatable :: l_ells
    real, dimension(:,:), allocatable :: test
    real, dimension(:), allocatable :: tau

    G = size(St,1)
    allocate(ells(NK,Nmats))
    allocate(I1inv(NKe,NKe))
    allocate(T(NKe,N))
    allocate(rquad(N))
    allocate(phi(N,G))
    allocate(l_ells(N,Nmats))
    allocate(tau(NKe))

    L = rglobal(1,1)

    call populate_Gauss_Legendre_quadrature(N, xq, wq)

    ells = 0
    do k = 2, NK
        ells(k,:) = ells(k-1,:)
        ells(k,eltomat(k-1)) = ells(k,eltomat(k-1)) + rglobal(1,k-1) - rglobal(1,k)
    end do

    I1inv = 2*(-1+3*identity(NKe))

    ! Node 1 is still e
    ! Node 2 is e + 1 (so smaller value of r)
    do p = 1, N
        T(1,p) = wq(p)*xq(p)
        T(2,p) = wq(p)*(1-xq(p))
    end do

    do e = 1, NE
        do k = 1, NKe
            tau(k) = sum(St(1,:)*ells(e+k-1,:))
        end do
        if (minval(tau) .lt. 0.5 .and. maxval(tau) .gt. 0.5) then
            do p = 1, N
                phi(p,:) = phiu(1,1,:) ! Set it equal to phiu at bdy, which has no attn.
                rquad(p) = 0.5*dz(e)*xq(p) + &
                           0.5*(rglobal(1,e) + rglobal(1,e+1))
                l_ells(p,:) = ells(e,:)
                l_ells(p,eltomat(e)) = l_ells(p,eltomat(e)) + rglobal(1,e) - rquad(p)

                do mat = 1, Nmats
                    phi(p,:) = phi(p,:)*exp(-St(:,mat)*l_ells(p,mat))
                end do
            end do

            !print *, St(1,1)
            !print *, rglobal(1,e), rglobal(1,e+1)
            !print *, 0.5*dz(e)*matmul(T,phi(:,1))
            !stop

            ! Must figure out the factors due to shape func. Factor due to I1inv is 1/dz(e)
            ! Factor is determinant. Should be 0.5*dz??
            phiu(e,:,:) = 0.5*matmul(I1inv,matmul(T,phi))
        end if
    end do

    !! Homogenizing seems to cause no problems. Even makes it more accurate, at least in 1D 1 MeV e- on Al.
    !allocate(test(NK,G))
    !test = 0
    !do k = 1, NKe
    !    do e = 1, NE
    !        test(e+k-1,:) = test(e+k-1,:) + phiu(e,k,:)/merge(1,2,e+k-1 .eq. 1 .or. e+k-1 .eq. NK)
    !    end do
    !end do
    !
    !do k = 1, NKe
    !    do e = 1, NE
    !        phiu(e,k,:) = test(e+k-1,:)
    !    end do
    !end do
end subroutine wiscoslab_MGXS_beam_quadrature

subroutine wiscoslab_MGXS_uncollided_fluence_correction &
    (r, vol, Ep, fEspec, fE, Sigmat, phiu)
    implicit none
    real, dimension(:,:,:), intent(in) :: r
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE
    real, dimension(:,:), intent(in) :: Sigmat

    real, dimension(:,:,:), intent(inout) :: phiu

    integer :: e, g, i, k
    integer :: l_G
    real, dimension(:), allocatable :: fEg
    real :: RL
    real :: RT
    real :: RTnew
    real :: alpha

    l_G = size(Sigmat,1)

    allocate(fEg(l_G))
    if (beam_energy_dist .eq. "polychromatic") then
        do g = 1, l_G
            call trapezoidal_integration_1D &
                (0, fE, fEspec, Ep(g+1), Ep(g), fEg(g))
        end do
    else if (beam_energy_dist .eq. "boxcar") then
        do g = 1, boxg
            fEg(g) = (Ep(g)-Ep(g+1))/deltaE
        end do
    else
        fEg = 0
        fEg(1) = 1
    end if

    if (beam_energy_dist .eq. "boxcar") l_G = boxg
    if (beam_energy_dist .eq. "monochromatic") l_G = 1
    do g = 1, l_G
        RL = 0
        RT = 0

        RL = RL - phiu(NE,2,g) ! Basically zero though

        do k = 1, NKe
            do e = 1, NE
                RT = RT + Sigmat(g,eltomat(e))*vol(e)*phiu(e,k,g)/NKe
            end do
        end do

        RTnew = fEg(g) - RL

        alpha = (RTnew - RT)/(RT**2)

        do e = 1, NE
            phiu(e,1:NKe,g) = phiu(e,1:NKe,g)*&
                (1+alpha*vol(e)*Sigmat(g,eltomat(e))*phiu(e,1:NKe,g))
        end do
    end do
    !stop
end subroutine wiscoslab_MGXS_uncollided_fluence_correction

subroutine wiscoslab_SN_MGXS_uncollided_source(mu, Sigma, phiu, source)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (gp,g,l,mat)
    real, dimension(:,:,:), intent(in) :: phiu ! (e,k,g)

    real, dimension(:,:,:,:), allocatable, intent(inout) :: source ! (e,k,i,g)

    integer :: e, g, i, k, l
    integer :: l_G
    real, dimension(:,:), allocatable :: P

    l_G = size(phiu,3)

    allocate(source(NE,NKe,Nmu,l_G))
    allocate(P(Nmu,0:NL))

    P = Legendre_polynomials(NL,mu)

    if (Nmats .eq. 1) then
        do g = 1, l_G
            do i = 1, Nmu
                do k = 1, NKe
                    do l = 0, NL
                        source(1:NE,k,i,g) = source(1:NE,k,i,g) + &
                            ((-1)**l)*P(i,l)*matmul(phiu(1:NE,k,1:g),Sigma(1:g,g,l+1,1))
                    end do
                end do
            end do
        end do
    else
        do g = 1, l_G
            do i = 1, Nmu
                do e = 1, NE
                    do l = 0, NL
                        source(e,1:NKe,i,g) = source(e,1:NKe,i,g) + &
                            ((-1)**l)*P(i,l)*matmul(phiu(e,1:NKe,1:g),Sigma(1:g,g,l+1,eltomat(e)))
                    end do
                end do
            end do
        end do
    end if
end subroutine wiscoslab_SN_MGXS_uncollided_source

subroutine wiscoslab_SN_MGXS_boundary_source_term(mu, Fsweep, E, fEspec, fE, bdy, fEg)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:), intent(in) :: E
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE

    real, dimension(:,:,:), allocatable, intent(inout) :: bdy
    real, dimension(:), allocatable :: fEg

    integer :: g, i, l
    integer :: l_G
    real, dimension(:,:), allocatable :: P
    real, dimension(:,:,:), allocatable :: ang_fl

    l_G = size(E)-1

    allocate(bdy(NKe,2,Nmu)) ! (k,bdyels,i)
    allocate(fEg(l_G))
    allocate(P(Nmu,0:NL))
    allocate(ang_fl(Nmu,2,NKe)) ! (i,bdyels,k)

    P = Legendre_polynomials(NL,mu)
    ang_fl = 0
    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            ang_fl(i,2,2) = 0
            !do l = 0, NL
            !    ang_fl(i,2,2) = ang_fl(i,2,2) + &
            !        (2*l+1)*P(i,l)*((-1)**l)/fourpi
            !end do
        else
            do l = 0, NL
                ang_fl(i,1,1) = ang_fl(i,1,1) + &
                    (2*l+1)*P(i,l)*((-1)**l)/fourpi
            end do
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

    if (beam_energy_dist .eq. "polychromatic") then
        do g = 1, l_G
            call trapezoidal_integration_1D &
                (0, fE, fEspec, E(g+1), E(g), fEg(g))
        end do
    else if (beam_energy_dist .eq. "boxcar") then
        fEg = 0
        do g = 1, boxg
            fEg(g) = (E(g)-E(g+1))/deltaE
        end do
    else
        fEg = 0
        fEg(1) = 1
    end if
end subroutine wiscoslab_SN_MGXS_boundary_source_term

subroutine wiscoslab_SN_FEXS_boundary_source_term(mu, Fsweep, E, fEspec, fE, bdy, fEg)
    implicit none
    real, dimension(:), intent(in) :: mu
    real, dimension(:,:,:,:,:), intent(in) :: Fsweep
    real, dimension(:), intent(in) :: E
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE

    real, dimension(:,:,:), allocatable, intent(inout) :: bdy
    real, dimension(:,:), allocatable :: fEg

    integer :: g, i, l, n
    integer :: l_G
    real, dimension(:,:), allocatable :: P
    real, dimension(:,:,:), allocatable :: ang_fl

    l_G = size(E)-1

    allocate(bdy(NKe,2,Nmu)) ! (k,bdyels,i)
    allocate(fEg(2,l_G))
    allocate(P(Nmu,0:NL))
    allocate(ang_fl(Nmu,2,NKe)) ! (i,bdyels,k)

    P = Legendre_polynomials(NL,mu)
    ang_fl = 0
    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            ang_fl(i,2,2) = 0
            !do l = 0, NL
            !    ang_fl(i,2,2) = ang_fl(i,2,2) + &
            !        (2*l+1)*P(i,l)*((-1)**l)/fourpi
            !end do
        else
            do l = 0, NL
                ang_fl(i,1,1) = ang_fl(i,1,1) + &
                    (2*l+1)*P(i,l)*((-1)**l)/fourpi
            end do
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

    if (beam_energy_dist .eq. "polychromatic") then
        do g = 1, l_G
            do n = 1, 2
                fEg(n,g) = interp1d(fEspec,fE,E(g+n-1))
            end do
        end do
    else if (beam_energy_dist .eq. "boxcar") then
        fEg = 0
        do g = 1, boxg
            do n = 1, 2
                fEg(n,g) = 1.0/deltaE
            end do
        end do
    else
        fEg = 0
        fEg(1,1) = 1/(E(1)-E(2))
        fEg(2,1) = 1/(E(1)-E(2))
    end if
end subroutine wiscoslab_SN_FEXS_boundary_source_term

end module ele_sources