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

module angular_discretization
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

contains
subroutine SN_angular_quadrature_abscissae_and_weights(l_Nmu, l_Nphi, xL, wL, phiC)
    implicit none
    integer, intent(in) :: l_Nmu
    integer, intent(in) :: l_Nphi

    real, dimension(:), allocatable, intent(inout) :: xL
    real, dimension(:), allocatable, intent(inout) :: wL
    real, dimension(:), allocatable, intent(inout) :: phiC

    integer :: j

    allocate(phiC(l_Nphi))

    call populate_Gauss_Legendre_quadrature(l_Nmu, xL, wL)

    do j = 1, l_Nphi
        phiC(j) = (2*j-1)*pi/l_Nphi
    end do
end subroutine SN_angular_quadrature_abscissae_and_weights

subroutine construct_discrete_ordinates(xL, phiC, khat)
    implicit none
    real, dimension(:), intent(in) :: xL
    real, dimension(:), intent(in) :: phiC

    real, dimension(:,:,:), allocatable, intent(inout) :: khat

    integer :: i, j
    integer :: l_Nmu
    integer :: l_Nphi

    l_Nmu = size(xL)
    l_Nphi = size(phiC)

    allocate(khat(3,l_Nmu,l_Nphi))

    do j = 1, l_Nphi
        do i = 1, l_Nmu
            khat(1:3,i,j) = [cos(phiC(j))*sqrt(1-xL(i)**2),sin(phiC(j))*sqrt(1-xL(i)**2),xL(i)]
        end do
    end do
end subroutine construct_discrete_ordinates

subroutine SN_scattering_angular_quadrature_terms(xL, phiC, cosmat, Pa, factorialmat, qPa, qfactorialmat)
    implicit none
    real, dimension(:), intent(in) :: xL
    real, dimension(:), intent(in) :: phiC

    real, dimension(:,:,:), allocatable, intent(inout) :: cosmat ! (j,jp,m)
    real, dimension(:,:,:), allocatable, intent(inout) :: Pa ! (i,l,m)
    real, dimension(:,:), allocatable, intent(inout) :: factorialmat ! (l,m)
    real, dimension(:,:), allocatable, intent(inout) :: qPa ! (i,q)
    real, dimension(:), allocatable, intent(inout) :: qfactorialmat ! (q)

    integer :: j, jp, l, m, q
    integer :: NQ
    integer :: l_Nmu
    integer :: l_Nphi
    real, dimension(:,:,:), allocatable :: temp_Pa ! (i,m,l)
    integer, dimension(:), allocatable :: ell

    l_Nmu = size(xL)
    l_Nphi = size(phiC)
    NQ = NL*(NL+3)/2+1
    allocate(cosmat(l_Nphi,l_Nphi,NL+1))
    allocate(Pa(l_Nmu,NL+1,NL+1))
    allocate(factorialmat(NL+1,NL+1))
    allocate(qPa(l_Nmu,NQ))
    allocate(qfactorialmat(NQ))
    allocate(temp_Pa(l_Nmu,NL+1,NL+1))
    allocate(ell(NQ))

    cosmat(1:l_Nphi,1:l_Nphi,1) = 1
    do m = 1, NL
        do jp = 1, l_Nphi
            do j = 1, l_Nphi
                cosmat(j,jp,m+1) = 2*cos(m*(phiC(j)-phiC(jp)))
            end do
        end do
    end do

    do q = 1, NQ
        ell(q) = floor(0.5*sqrt(8.0*q-7.0)-0.5)
    end do

    temp_Pa = assoc_Legendre_polynomials(NL,xL)

    Pa = reshape(temp_Pa, [l_Nmu,NL+1,NL+1], order = [1,3,2])

    factorialmat = 0
    do m = 0, NL
        do l = m, NL
            factorialmat(l+1,m+1) = factorial(l-m)/factorial(l+m)
        end do
    end do

    do l = 0, NL
        do m = 0, l
            qPa(1:Nmu,l*(l+1)/2+m+1) = Pa(1:Nmu,l+1,m+1)
            qfactorialmat(l*(l+1)/2+m+1) = factorialmat(l+1,m+1)
        end do
    end do
end subroutine SN_scattering_angular_quadrature_terms

subroutine associated_Legendre_polynomials(xL, Pa)
    implicit none
    real, dimension(:), intent(in) :: xL

    real, dimension(:,:,:), allocatable, intent(inout) :: Pa

    integer :: i, l, m
    real, dimension(:,:,:), allocatable :: temp_Pa
    integer :: l_Nmu

    l_Nmu = size(xL)
    allocate(Pa(l_Nmu,NL+1,NL+1))
    allocate(temp_Pa(0:NL,0:NL,l_Nmu))

    temp_Pa = assoc_Legendre_polynomials(NL,xL)

    do i = 1, l_Nmu
        do m = 0, NL
            do l = 0, NL
                Pa(i,m+1,l+1) = temp_Pa(l,m,i)
            end do
        end do
    end do
end subroutine associated_Legendre_polynomials

subroutine cosmatrix(phiC, cosmat, factorialmat)
    implicit none
    real, dimension(:), intent(in) :: phiC

    real, dimension(:,:,:), allocatable, intent(inout) :: cosmat
    real, dimension(:,:), allocatable, intent(inout) :: factorialmat

    integer :: j, jp, l, m
    integer :: l_Nphi

    l_Nphi = size(phiC)
    allocate(cosmat(NL+1,l_Nphi,l_Nphi))
    allocate(factorialmat(NL+1,NL+1))

    cosmat(1,:,:) = 1
    do jp = 1, l_Nphi
        do j = 1, l_Nphi
            do m = 1, NL
                cosmat(m+1,j,jp) = 2*cos(m*(phiC(j)-phiC(jp)))
            end do
        end do
    end do

    factorialmat = 0
    do l = 0, NL
        do m = 0, l
            factorialmat(m+1,l+1) = factorial(l-m)/factorial(l+m)
        end do
    end do
end subroutine cosmatrix

subroutine SN_uncollided_angular_quadrature_terms(rglobal, phiC, Pmlk, cmjk)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:), intent(in) :: phiC

    real, dimension(:,:,:), allocatable, intent(inout) :: Pmlk ! (l,m,k)
    real, dimension(:,:,:), allocatable, intent(inout) :: cmjk ! (j,m,k)

    integer :: j, k, m
    real, dimension(3) :: k0
    real, dimension(3) :: rcurs
    real, dimension(:), allocatable :: xvals
    real, dimension(:,:,:), allocatable :: temp_Pmlk ! (k,m,l)
    real :: phival

    k0 = beam_coord_sys(:,3)
    if (beam_angular_dist .eq. "spherical") then
        allocate(Pmlk(NL+1,NL+1,size(nodesinbeam)))
        allocate(cmjk(Nphi,NL+1,size(nodesinbeam)))
        allocate(xvals(size(nodesinbeam)))
        allocate(temp_Pmlk(size(nodesinbeam),NL+1,NL+1))

        do k = 1, size(nodesinbeam)
            rcurs = rglobal(:,nodesinbeam(k))-beam_origin
            xvals(k) = rcurs(3)/norm2(rcurs)
        end do

        temp_Pmlk = assoc_Legendre_polynomials(NL,xvals)

        Pmlk = reshape(temp_Pmlk, [NL+1,NL+1,size(nodesinbeam)], order = [3,2,1])

        cmjk(:,1,:) = 1
        do k = 1, size(nodesinbeam)
            rcurs = rglobal(:,nodesinbeam(k))-beam_origin
            if (rcurs(1) .eq. 0. .and. rcurs(2) .eq. 0.) then
                cmjk(:,2:NL+1,k) = 0. ! If the node has no meaningful phi value (|mu| = 1), the associated
                ! Legendre polynomial Pmlk at this node will kill the whole term anyway, for m > 0.
                ! If m = 0, then cos(m(...)) = 1, regardless.
                cycle
            end if
            do m = 1, NL
                do j = 1, Nphi
                    phival = atan2(rcurs(2)/hypot(rcurs(1),rcurs(2)),&
                                   rcurs(1)/hypot(rcurs(1),rcurs(2)))
                    cmjk(j,m+1,k) = 2*cos(m*(phiC(j)-phival))
                end do
            end do
        end do
    else
        allocate(Pmlk(NL+1,NL+1,1))
        allocate(cmjk(Nphi,NL+1,1))
        allocate(xvals(1))
        allocate(temp_Pmlk(1,NL+1,NL+1))

        xvals = k0(3)

        temp_Pmlk = assoc_Legendre_polynomials(NL,xvals)

        Pmlk = reshape(temp_Pmlk, [NL+1,NL+1,1], order = [3,2,1])

        cmjk(:,1,:) = 1
        phival = atan2(k0(2),k0(1))
        if (abs(atan2(k0(2)/hypot(k0(1),k0(2)),k0(1)/hypot(k0(1),k0(2))) - atan2(k0(2),k0(1))) .ge. 1.0E-8) then
            print *, "SN_uncollided_angular_quadrature_terms"
            print *, "ATAN2 DOES NEED EITHER NORMALIZED OR NOT NORMALIZED VECTOR."
            stop
        end if
        do m = 1, NL
            cmjk(1:Nphi,m+1,1) = 2*cos(m*(phiC-phival))
        end do
    end if
end subroutine SN_uncollided_angular_quadrature_terms

subroutine uncollided_RSH_terms(l_NL, rglobal, RSH_unc)
    implicit none
    integer, intent(in) :: l_NL
    real, dimension(:,:), intent(in) :: rglobal

    real, dimension(:,:), allocatable, intent(inout) :: RSH_unc

    integer :: j, k, m
    integer :: NQ
    real, dimension(3) :: k0
    real, dimension(3) :: rcurs
    real, dimension(:,:), allocatable :: vectors

    NQ = (l_NL+1)**2

    k0 = beam_coord_sys(:,3)
    if (beam_angular_dist .eq. "spherical") then
        allocate(vectors(3,size(nodesinbeam)))

        do k = 1, size(nodesinbeam)
            rcurs = rglobal(:,nodesinbeam(k))-beam_origin
            vectors(:,k) = rcurs/norm2(rcurs)
        end do

        allocate(RSH_unc(size(nodesinbeam),NQ))

        RSH_unc = q_ind_real_spherical_harmonics_vector_input(NQ,vectors)
    else
        allocate(vectors(3,1))
        vectors(:,1) = k0

        allocate(RSH_unc(1,NQ))

        RSH_unc = q_ind_real_spherical_harmonics_vector_input(NQ,vectors)
    end if
end subroutine uncollided_RSH_terms

subroutine PN_angular_integrals(l_NL, intfc, normal, Awl, Aup, Adn)
    ! Must redo from scratch with new derivations
    implicit none
    integer, intent(in) :: l_NL
    integer, dimension(:,:,:), intent(in) :: intfc ! (fg, el or face, prim or sec)
    real, dimension(:,:), intent(in) :: normal ! (dir,f,e)

    real, dimension(:,:,:), allocatable, intent(inout) :: Awl ! (q,qp,x)
    real, dimension(:,:,:), allocatable, intent(inout) :: Aup ! (q,qp,fg)
    real, dimension(:,:,:), allocatable, intent(inout) :: Adn ! (q,qp,fg)

    integer :: fg, i, k, l, lp, m, mp, n, x
    integer :: l_NF
    integer :: NQ
    real, dimension(4) :: kappas
    integer, dimension(4) :: kappas_ind
    integer :: qval1
    integer :: qval2
    integer :: minkl
    real, dimension(:,:,:), allocatable :: LPints
    real, dimension(:,:,:,:), allocatable :: Wigner
    real, dimension(:,:,:), allocatable :: HSints
    real, dimension(:,:,:,:), allocatable :: Hmat
    real :: localstart
    real :: localfinish

    l_NF = size(intfc,1)
    NQ = (l_NL+1)**2

    allocate(Awl(NQ,4,3))
    allocate(Aup(NQ,NQ,l_NF))
    allocate(Adn(NQ,NQ,l_NF))
    allocate(LPints(-l_NL-1:l_NL+1,0:l_NL+1,0:l_NL+1)) ! (mpp,k,lp) Max value of k,lp is increased for completion of Aup recursion relations
    allocate(Wigner(l_NF,-l_NL-1:l_NL+1,0:l_NL+1,-l_NL-1:l_NL+1)) ! (fg,mp,l,m) Max value of mp,l,m increased, see above
    allocate(HSints(l_NF,NQ,(l_NL+2)**2)) ! (fg,qlp,qk) Max value of qk is increased, see above
    allocate(Hmat(NQ,l_NF,(l_NL+2)**2,3)) ! (qlp,fg,qk,x) Max value of qk is increased, see above

    call cpu_time(localstart)
    LPints = 0.
    do l = 0, l_NL+1
        do k = 0, l-1
            if (even(l+k)) cycle
            do m = -k, k
                LPints(m,k,l) = RSH_norm(l,m)*RSH_norm(k,m)*AL_polys_half_integrals(m,k,l)
            end do
            LPints(:,l,k) = LPints(:,k,l)
        end do
    end do
    call cpu_time(localfinish)
    print *, "LOOP 1:", localfinish-localstart

    call cpu_time(localstart)
    Wigner = 0.
    do m = -l_NL-1, l_NL+1
        do l = abs(m), l_NL+1
            do mp = -l, l
                do fg = 1, l_NF
                    Wigner(fg,mp,l,m) = &
                    Wigner_RD_matrix_elements(l,mp,m,0.,-acos(normal(3,fg)),-atan2(normal(2,fg),normal(1,fg)))
                end do
            end do
        end do
    end do
    call cpu_time(localfinish)
    print *, "LOOP 2:", localfinish-localstart

    ! LOOP 2 IS SLOW. NEED TO MAKE IT NOT REDO ENTIRE (l,mp,m)-DEPENDENT TERMS FOR EACH NORMAL VECTOR

    call cpu_time(localstart)
    HSints = 0.
    do k = 0, l_NL+1
        do lp = 0, l_NL
            if (even(lp+k)) then
                if (lp == k) then
                    do mp = -lp, lp
                        HSints(:,lp*(lp+1)+mp+1,lp*(lp+1)+mp+1) = .5
                    end do
                end if
            else
                minkl = min(k,lp)
                do n = -k, k
                    qval2 = k*(k+1)+n+1
                    do mp = -lp, lp
                        qval1 = lp*(lp+1)+mp+1
                        HSints(:,qval1,qval2) = &
                        matmul(Wigner(:,-minkl:minkl,k,n)*Wigner(:,-minkl:minkl,lp,mp),LPints(-minkl:minkl,lp,k))
                    end do
                end do
            end if
        end do
    end do
    call cpu_time(localfinish)
    print *, "LOOP 3:", localfinish-localstart

    call cpu_time(localstart)
    Hmat = 0.
    do x = 1, 3
        do k = 0, l_NL+1
            do n = -k, k
                qval2 = k*(k+1)+n+1
                do lp = 0, l_NL
                    do mp = -lp, lp
                        qval1 = lp*(lp+1)+mp+1
                        Hmat(qval1,:,qval2,x) = HSints(:,qval1,qval2)*normal(x,:)
                    end do
                end do
            end do
        end do
    end do
    call cpu_time(localfinish)
    print *, "LOOP 4:", localfinish-localstart

    call cpu_time(localstart)
    Awl = 0
    Aup = 0
    do l = 0, l_NL
        do m = -l, l
            qval1 = l*(l+1)+m+1
            do x = 1, 2
                kappas = RSH_kappas(l,m,x)
                kappas_ind = RSH_kappas_indices(l,m,x)
                do i = 1, 4
                    if (kappas_ind(i) == 0) cycle
                    Awl(qval1,i,x) = kappas(i)
                    Aup(qval1,:,:) = Aup(qval1,:,:) + kappas(i)*Hmat(:,:,kappas_ind(i),x) ! THIS GIVES fg ON OUTSIDE
                end do
            end do
            kappas(1:2) = RSH_kappas(l,m,3)
            kappas_ind(1:2) = RSH_kappas_indices(l,m,3)
            do i = 1, 2
                if (kappas_ind(i) == 0) cycle
                Awl(qval1,i,3) = kappas(i)
                Aup(qval1,:,:) = Aup(qval1,:,:) + kappas(i)*Hmat(:,:,kappas_ind(i),3)
            end do
        end do
    end do
    call cpu_time(localfinish)
    print *, "LOOP 5:", localfinish-localstart
end subroutine PN_angular_integrals

subroutine PN_angular_integrals_via_quadrature(l_NL, Ceff, globalnormal, sgn, Nint, Awl, Ia, Ja, g_Aup, g_Adn)
    implicit none
    integer, intent(in) :: l_NL
    integer, dimension(:,:), intent(in) :: Ceff
    real, dimension(:,:), intent(in) :: globalnormal
    logical, dimension(:,:), intent(in) :: sgn
    integer, intent(in) :: Nint

    !real, dimension(:,:), allocatable, intent(inout) :: Awl ! (q,qp,dir)
    real, dimension(:,:,:), allocatable, intent(inout) :: Awl
    ! Must store this in sparse storage. Just do (q,i,x), where i = 1, ..., 4, and have some convention for picking the i's. Maybe be able to forego kappas_ind. (~2/13/23)
    integer, dimension(:,:), allocatable, intent(inout) :: Ia
    integer, dimension(:,:), allocatable, intent(inout) :: Ja
    real, dimension(:,:,:), allocatable, intent(inout) :: g_Aup
    real, dimension(:,:,:), allocatable, intent(inout) :: g_Adn
    !real, dimension(:,:,:,:), allocatable, intent(inout) :: Aup ! (q,qp,f,e)
    !real, dimension(:,:,:,:), allocatable, intent(inout) :: Adn ! (q,qp,f,e)
    !real, dimension(:,:,:,:), allocatable, intent(inout) :: A2 ! (q,qp,fg)

    integer :: dir, e, f, i, l, m, n, q, qp
    integer :: NQ
    integer :: NF
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:), allocatable :: dense_Awl ! (q,qp,dir)
    integer, dimension(:,:,:), allocatable :: kappas_ind
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:), allocatable :: xqGJ
    real, dimension(:), allocatable :: wqGJ
    real, dimension(:,:), allocatable :: Legendre ! (i,q)
    real, dimension(:,:), allocatable :: LegendreGJ ! (i,q)
    real, dimension(:), allocatable :: RSH_norms
    real :: phi0
    real :: l_nz
    real :: sin0
    real :: sign0
    real, dimension(:), allocatable :: alphas
    real, dimension(:), allocatable :: betas
    real, dimension(:,:), allocatable :: phivals
    real, dimension(:,:,:), allocatable :: temp_Sm
    real, dimension(:,:,:), allocatable :: temp_Sn
    real, dimension(:,:,:), allocatable :: SmSn
    real, dimension(:,:,:), allocatable :: INTS
    integer :: qval
    integer :: qpval
    integer :: lval
    integer :: mval
    real, dimension(:,:,:), allocatable :: kappas
    integer :: etest, ftest
    real :: mansgn
    real, dimension(:), allocatable :: temp_Awl_1
    real, dimension(:), allocatable :: temp_Awl_2
    real, dimension(:), allocatable :: temp_Awl_3
    integer, dimension(:), allocatable :: temp_Ia_1
    integer, dimension(:), allocatable :: temp_Ia_2
    integer, dimension(:), allocatable :: temp_Ia_3
    integer, dimension(:), allocatable :: temp_Ja_1
    integer, dimension(:), allocatable :: temp_Ja_2
    integer, dimension(:), allocatable :: temp_Ja_3
    integer :: Awl_size
    integer :: max_dim

    NQ = (l_NL+2)**2

    NF = size(globalnormal,2)

    allocate(ell(NQ))
    allocate(emm(NQ))
    allocate(xqGJ(Nint))
    allocate(wqGJ(Nint))
    allocate(Legendre(Nint,NQ))
    allocate(LegendreGJ(Nint,NQ))
    allocate(RSH_norms(NQ))
    allocate(alphas(Nint))
    allocate(betas(Nint))
    allocate(phivals(Nint,Nint))
    allocate(temp_Sm(Nint,Nint,-l_NL-1:l_NL+1))
    allocate(temp_Sn(Nint,Nint,-l_NL-1:l_NL+1))
    allocate(SmSn(Nint,-l_NL-1:l_NL+1,-l_NL-1:l_NL+1))
    allocate(INTS(NQ,NQ,NF))

    call populate_Gauss_Legendre_quadrature(Nint, xq, wq)
    do i = 1, Nint
        xqGJ(i) = cos(i*pi/(Nint+1))
        wqGJ(i) = pi*sin(i*pi/(Nint+1))**2/(Nint+1)
    end do

    Legendre = q_ind_associated_Legendre_polynomials(NQ,xq)
    LegendreGJ = q_ind_associated_Legendre_polynomials(NQ,xqGJ)
    do i = 1, Nint
        LegendreGJ(i,:) = LegendreGJ(i,:)/(1-xqGJ(i)**2)**(0.25)
    end do
    ! LegendreGJ will only be used when LegendreGJ(q)*LegendreGJ(qp) has spare factor of sqrt(1-xqGJ**2).
    ! Therefore LegendreGJ(q)*LegendreGJ(qp) is a pure polynomial.

    do q = 1, NQ
        ell(q) = ceiling(sqrt(real(q)))-1
        emm(q) = q - ell(q)*(ell(q)+1) - 1
        RSH_norms(q) = RSH_norm(ell(q),emm(q))
    end do

    INTS = 0
    !$OMP parallel do shared(INTS) private(f, m, n, qp, q) &
    !$OMP private(phi0, l_nz, sin0, sign0, alphas, betas, phivals, temp_Sm, SmSn)
    do f = 1, NF
        phi0 = atan2(globalnormal(2,f),globalnormal(1,f))
        l_nz = globalnormal(3,f)
        sin0 = sqrt(1-l_nz**2)
        sign0 = sign(1.0,l_nz)
        alphas = phi0 - pi/2
        betas = phi0 + pi/2

        where (xq .le. -sin0)
            alphas = alphas + sign0*pi/2
            betas = betas - sign0*pi/2
        end where
        where (-sin0 .lt. xq .and. xq .lt. sin0) ! Use (abs(xq) <= sin0)? (~3/26/23)
            alphas = alphas - asin(xq*l_nz/(sqrt(1-xq**2)*sin0))
            betas = betas + asin(xq*l_nz/(sqrt(1-xq**2)*sin0))
        end where
        where (sin0 .le. xq)
            alphas = alphas - sign0*pi/2
            betas = betas + sign0*pi/2
        end where

        do i = 1, Nint
            phivals(1:Nint,i) = 0.5*(betas(i)-alphas(i))*xq + 0.5*(betas(i)+alphas(i))
        end do

        do m = -l_NL-1, -1
            temp_Sm(1:Nint,1:Nint,m) = sin(abs(m)*transpose(phivals))
        end do
        do m = 0, l_NL+1
            temp_Sm(1:Nint,1:Nint,m) = cos(m*transpose(phivals))
        end do
        do m = -l_NL-1, l_NL+1
            do n = m, l_NL+1
                SmSn(1:Nint,n,m) = matmul(temp_Sm(1:Nint,1:Nint,m)*temp_Sm(1:Nint,1:Nint,n),wq)
                SmSn(1:Nint,m,n) = SmSn(1:Nint,n,m)
            end do
        end do

        do qp = 1, NQ ! MAKE THIS SKIP EVEN l+k PAIRS
            do q = qp, NQ
                if (even(ell(q) + ell(qp))) then
                    INTS(q,qp,f) = 0.5*Kron_delta(q,qp)
                    INTS(qp,q,f) = INTS(q,qp,f)
                else
                ! FINISH IMPLEMENTING! Problem is that SmSn integrals are eval'd at GL nodes, NOT GJ nodes. Would need to do SmSn for both to make this work.
                ! Is it even worth it? Maybe not...
                !if (odd(emm(q)+emm(qp))) then
                !    INTS(q,qp,f) = 0.5*RSH_norms(q)*RSH_norms(qp)*&
                !        sum(wq*SmSn(1:Nint,emm(q),emm(qp))*(betas-alphas)*Legendre(1:Nint,q)*Legendre(1:Nint,qp))
                !else
                    INTS(q,qp,f) = 0.5*RSH_norms(q)*RSH_norms(qp)*&
                        sum(wq*SmSn(1:Nint,emm(q),emm(qp))*(betas-alphas)*Legendre(1:Nint,q)*Legendre(1:Nint,qp))
                    INTS(qp,q,f) = INTS(q,qp,f)
                end if
            end do
        end do

        if (l_nz .eq. 1) then
            do q = 1, NQ
                do qp = 1, NQ
                    print *, INTS(q,qp,f)
                end do
            end do
            stop
        end if
    end do
    !$OMP end parallel do

    !allocate(dense_Awl(NQ,NQ,3))
    allocate(Awl(NQ,NQ,3))

    NQ = (l_NL+1)**2

    allocate(kappas_ind(4,NQ,3))
    allocate(kappas(4,NQ,3))
    !allocate(A2(NQ,NQ,3,3))

    Awl = 0
    kappas_ind = 0
    kappas = 0
    do l = 0, l_NL
        do m = -l, l
            qval = l*(l+1)+m+1
            do dir = 1, 2
                kappas_ind(1:4,qval,dir) = RSH_kappas_indices(l,m,dir)
                kappas(1:4,qval,dir) = RSH_kappas(l,m,dir)
                do i = 1, 4
                    if (kappas_ind(i,qval,dir) .eq. 0) cycle
                    Awl(qval,kappas_ind(i,qval,dir),dir) = kappas(i,qval,dir)
                    Awl(kappas_ind(i,qval,dir),qval,dir) = Awl(qval,kappas_ind(i,qval,dir),dir)
                end do
            end do
            kappas_ind(1:2,qval,3) = RSH_kappas_indices(l,m,3)
            kappas(1:2,qval,3) = RSH_kappas(l,m,3)
            do i = 1, 2
                if (kappas_ind(i,qval,3) .eq. 0) cycle
                Awl(qval,kappas_ind(i,qval,3),3) = kappas(i,qval,3)
                Awl(kappas_ind(i,qval,3),qval,3) = Awl(qval,kappas_ind(i,qval,3),3)
            end do
        end do
    end do

    allocate(g_Aup(NF,NQ,NQ))
    allocate(g_Adn(NF,NQ,NQ))
    g_Aup = 0
    g_Adn = 0
    ! MUST SWITCH BACK TO q,q,f IF I AM GOING BACK TO ELEMENTAL AUP AND ADN
    !do dir = 1, 3
    !    do q = 1, NQ
    !        do qp = 1, NQ
    !            g_Aup(q,qp,:) = g_Aup(q,qp,:) + globalnormal(dir,:)*&
    !                matmul(transpose(INTS(:,qp,:)),dense_Awl(q,:,dir))
    !        end do
    !    end do
    !end do
    !$OMP parallel do shared(g_Aup) private(dir, q, qp)
    do q = 1, NQ
        do qp = q, NQ
            do dir = 1, 3
                g_Aup(:,q,qp) = g_Aup(:,q,qp) + globalnormal(dir,:)*&
                    matmul(transpose(INTS(:,qp,:)),Awl(:,q,dir))
            end do
        end do
    end do
    !$OMP end parallel do
    do q = 1, NQ
        do qp = q, NQ
            g_Aup(:,qp,q) = g_Aup(:,q,qp)
        end do
    end do

    !do f = 1, NF
    !    do dir = 1, 3
    !        g_Adn(1:NQ,1:NQ,f) = g_Adn(1:NQ,1:NQ,f) + globalnormal(dir,f)*dense_Awl(1:NQ,1:NQ,dir)
    !    end do
    !end do
    do f = 1, NF
        do dir = 1, 3
            g_Adn(f,1:NQ,1:NQ) = g_Adn(f,1:NQ,1:NQ) + globalnormal(dir,f)*Awl(1:NQ,1:NQ,dir)
        end do
    end do

    g_Adn = g_Adn - g_Aup

    !allocate(Aup(NE,NQ,NFe,NQ))
    !allocate(Adn(NE,NQ,NFe,NQ))
    !!$OMP parallel do collapse(2) shared(Aup, Adn) private(f, e)
    !do f = 1, NFe
    !    do e = 1, NE
    !        !Aup(e,1:NQ,f,1:NQ) = merge(1.0,-1.0,sgn(f,e))*g_Aup(1:NQ,1:NQ,Ceff(e,f))
    !        !Adn(e,1:NQ,f,1:NQ) = merge(1.0,-1.0,sgn(f,e))*g_Adn(1:NQ,1:NQ,Ceff(e,f))
    !        if (sgn(f,e)) then
    !            Aup(e,1:NQ,f,1:NQ) = g_Aup(1:NQ,1:NQ,Ceff(e,f))
    !            Adn(e,1:NQ,f,1:NQ) = g_Adn(1:NQ,1:NQ,Ceff(e,f))
    !        else
    !            Aup(e,1:NQ,f,1:NQ) = -g_Adn(1:NQ,1:NQ,Ceff(e,f))
    !            Adn(e,1:NQ,f,1:NQ) = -g_Aup(1:NQ,1:NQ,Ceff(e,f))
    !        end if
    !    end do
    !end do
    !!$OMP end parallel do

    !do q = 1, NQ
    !    if (any(abs(Awl(q,q,:)) .gt. 0)) then
    !        print *, "PROBLEM:", q
    !        stop
    !    end if
    !end do

    !!q = 63
    !!qp = 38
    !q = 6
    !qp = 11
    !f = 4110
    !print *, "INPUTS"
    !print *, atan2(globalnormal(2,f),globalnormal(1,f))
    !print *, globalnormal(:,f)
    !print *, q, qp
    !print *, "RESULTS"
    !print *, INTS(q,qp,f)
    !print *, Awl(q,qp,1)
    !print *, Awl(q,qp,2)
    !print *, Awl(q,qp,3)
    !print *, dot_product(Awl(q,qp,:),globalnormal(:,f))
    !print *, g_Aup(f,q,qp)
    !print *, g_Adn(f,q,qp)
    !stop

    !do j = 1, 3
    !    do i = 1, 3
    !        A2(1:NQ,1:NQ,i,j) = matmul(Awl(1:NQ,1:NQ,i),Awl(1:NQ,1:NQ,j))
    !    end do
    !end do
!
    !call make_SBR(dense_Awl(1:NQ,1:NQ,1), 1.0E-10, temp_Awl_1, temp_Ia_1, temp_Ja_1)
    !call make_SBR(dense_Awl(1:NQ,1:NQ,2), 1.0E-10, temp_Awl_2, temp_Ia_2, temp_Ja_2)
    !call make_SBR(dense_Awl(1:NQ,1:NQ,3), 1.0E-10, temp_Awl_3, temp_Ia_3, temp_Ja_3)
!
    !Awl_size = maxval([size(temp_Awl_1), size(temp_Awl_2) ,size(temp_Awl_3)])
    !max_dim = maxloc([size(temp_Awl_1), size(temp_Awl_2) ,size(temp_Awl_3)], dim = 1)
!
    !allocate(Awl(Awl_size,3))
    !allocate(Ia(NQ+1,3))
    !allocate(Ja(Awl_size,3))
!
    !Awl = 0
    !Awl(1:size(temp_Awl_1),1) = temp_Awl_1
    !Awl(1:size(temp_Awl_2),2) = temp_Awl_2
    !Awl(1:size(temp_Awl_3),3) = temp_Awl_3
!
    !Ia(1:NQ+1,1) = temp_Ia_1
    !Ia(1:NQ+1,2) = temp_Ia_2
    !Ia(1:NQ+1,3) = temp_Ia_3
!
    !Ja = NQ ! So that the padded 0's are considered to be in column 16. (Should I not do this?)
    !Ja(1:size(temp_Ja_1),1) = temp_Ja_1
    !Ja(1:size(temp_Ja_2),2) = temp_Ja_2
    !Ja(1:size(temp_Ja_3),3) = temp_Ja_3
end subroutine PN_angular_integrals_via_quadrature

subroutine verify_Awl
    implicit none
    integer :: q, qp, i, l, m, n, dir, qval, qpval
    integer :: NQ
    integer :: pNint
    integer :: mNint
    real, dimension(:,:,:), allocatable :: pINTS
    real, dimension(:,:,:), allocatable :: mINTS
    real, dimension(:,:,:), allocatable :: INTS
    real, dimension(:,:,:), allocatable :: SmSn
    real, dimension(:), allocatable :: pxq
    real, dimension(:), allocatable :: pwq
    real, dimension(:,:), allocatable :: Legendre
    real, dimension(:), allocatable :: xq
    real, dimension(:), allocatable :: wq
    real, dimension(:,:,:), allocatable :: Awl
    integer, dimension(:,:,:), allocatable :: kappas_ind
    real, dimension(:,:,:), allocatable :: kappas
    integer, dimension(:), allocatable :: emm
    integer, dimension(:), allocatable :: ell
    real :: MVAL
    real, dimension(:), allocatable :: tempL
    real, dimension(:), allocatable :: RSH_norms

    NQ = (NL+1)**2
    pNint = 2*(NL+1)
    mNint = (NL+1)

    allocate(pINTS(NQ,NQ,3))
    allocate(mINTS(NQ,NQ,3))
    allocate(INTS(NQ,NQ,3))
    allocate(SmSn(pNint,-NL:NL,-NL:NL))
    allocate(Legendre(mNint,NQ))
    allocate(ell(NQ))
    allocate(emm(NQ))
    allocate(tempL(mNint))
    allocate(RSH_norms(NQ))

    do q = 1, NQ
        ell(q) = ceiling(sqrt(real(q)))-1
        emm(q) = q - ell(q)*(ell(q)+1) - 1
    end do

    !call populate_Gauss_Legendre_quadrature(pNint, pxq, pwq)
    allocate(pxq(pNint))
    allocate(pwq(pNint))
    do i = 1, pNint
        pxq(i) = (2*i-1)*pi/pNint
        pwq(i) = twopi/pNint
    end do

    call populate_Gauss_Legendre_quadrature(mNint, xq, wq)

    Legendre = q_ind_associated_Legendre_polynomials(NQ,xq)

    do l = 0, NL
        do m = -l, l
            q = l*(l+1)+m+1
            RSH_norms(q) = RSH_norm(l,m)
        end do
    end do

    do m = -NL, NL
        do n = -NL, NL
            do i = 1, pNint
                !SmSn(i,m,n) = merge(cos(m*(0.5*twopi*pxq(i)+twopi)),&
                !                    sin(abs(m)*(0.5*twopi*pxq(i)+twopi)), m .ge. 0)*&
                !              merge(cos(n*(0.5*twopi*pxq(i)+twopi)),&
                !                    sin(abs(n)*(0.5*twopi*pxq(i)+twopi)), n .ge. 0)
                SmSn(i,m,n) = merge(cos(m*pxq(i)),&
                                    sin(abs(m)*pxq(i)), m .ge. 0)*&
                              merge(cos(n*pxq(i)),&
                                    sin(abs(n)*pxq(i)), n .ge. 0)
            end do
        end do
    end do

    pINTS = 0
    mINTS = 0
    do qp = 1, NQ
        do q = 1, NQ
            do i = 1, pNint
                !pINTS(q,qp,1) = pINTS(q,qp,1) + 0.5*twopi*pwq(i)*SmSn(i,emm(q),emm(qp))*&
                !    cos(0.5*twopi*pxq(i)+twopi)
                !pINTS(q,qp,2) = pINTS(q,qp,2) + 0.5*twopi*pwq(i)*SmSn(i,emm(q),emm(qp))*&
                !    sin(0.5*twopi*pxq(i)+twopi)
                !pINTS(q,qp,3) = pINTS(q,qp,3) + 0.5*twopi*pwq(i)*SmSn(i,emm(q),emm(qp))
                pINTS(q,qp,1) = pINTS(q,qp,1) + twopi*SmSn(i,emm(q),emm(qp))*&
                    cos(pxq(i))/pNint
                pINTS(q,qp,2) = pINTS(q,qp,2) + twopi*SmSn(i,emm(q),emm(qp))*&
                    sin(pxq(i))/pNint
                pINTS(q,qp,3) = pINTS(q,qp,3) + twopi*SmSn(i,emm(q),emm(qp))/pNint
            end do

            tempL = Legendre(:,q)*Legendre(:,qp)

            do i = 1, mNint
                mINTS(q,qp,1) = mINTS(q,qp,1) + wq(i)*sqrt(1-xq(i)**2)*tempL(i)
                mINTS(q,qp,3) = mINTS(q,qp,3) + wq(i)*xq(i)*tempL(i)
            end do
            mINTS(q,qp,2) = mINTS(q,qp,1)

            INTS(q,qp,:) = RSH_norms(q)*RSH_norms(qp)*pINTS(q,qp,:)*mINTS(q,qp,:)
        end do
    end do

    NQ = (NL+2)**2
    allocate(Awl(NQ,NQ,3))
    NQ = (NL+1)**2
    allocate(kappas_ind(4,NQ,3))
    allocate(kappas(4,NQ,3))

    Awl = 0
    kappas_ind = 0
    kappas = 0

    do l = 0, NL
        do m = -l, l
            qval = l*(l+1)+m+1
            do dir = 1, 2
                kappas_ind(1:4,qval,dir) = RSH_kappas_indices(l,m,dir)
                kappas(1:4,qval,dir) = RSH_kappas(l,m,dir)
                do i = 1, 4
                    if (kappas_ind(i,qval,dir) .eq. 0) cycle
                    Awl(qval,kappas_ind(i,qval,dir),dir) = kappas(i,qval,dir)
                    Awl(kappas_ind(i,qval,dir),qval,dir) = kappas(i,qval,dir)
                end do
            end do
            kappas_ind(1:2,qval,3) = RSH_kappas_indices(l,m,3)
            kappas(1:2,qval,3) = RSH_kappas(l,m,3)
            do i = 1, 2
                if (kappas_ind(i,qval,3) .eq. 0) cycle
                Awl(qval,kappas_ind(i,qval,3),3) = kappas(i,qval,3)
                Awl(kappas_ind(i,qval,3),qval,3) = kappas(i,qval,3)
            end do
        end do
    end do

    do q = 1, NQ
        do qp = q, NQ
            do dir = 1, 3
                !if (abs(Awl(q,qp,dir) - INTS(q,qp,dir)) .gt. 1.0E-8) then
                !    print *, [q,qp,dir], Awl(q,qp,dir), INTS(q,qp,dir)
                !end if
                !if (abs(pINTS(q,qp,1) - &
                !    (0.5*(1+Kron_delta(emm(q),0)-Kron_delta(emm(q),-1))*&
                !    Kron_delta(emm(qp),emm(q)+1) + &
                !    0.5*(1-Kron_delta(emm(q),0))*Kron_delta(emm(qp),emm(q)-1))*&
                !    pi*(1+Kron_delta(emm(qp),0))) .gt. 1.0E-8) then
                !    print *, [q,qp], (0.5*(1+Kron_delta(emm(q),0)-Kron_delta(emm(q),-1))*&
                !        Kron_delta(emm(qp),emm(q)+1) + &
                !        0.5*(1-Kron_delta(emm(q),0))*Kron_delta(emm(qp),emm(q)-1))*&
                !        pi*(1+Kron_delta(emm(qp),0)), pINTS(q,qp,1)
                !end if
                if (abs(pINTS(q,qp,3) - &
                    pi*(1+Kron_delta(emm(qp),0))*Kron_delta(emm(q),emm(qp))) .gt. 1.0E-8) then
                    print *, [q,qp], pINTS(q,qp,3), pi*(1+Kron_delta(emm(qp),0))*Kron_delta(emm(q),emm(qp))
                end if
            end do
        end do
    end do

    stop
end subroutine verify_Awl

subroutine populate_real_spherical_harmonics(l_NL, mu, phi, RSH)
    implicit none
    integer, intent(in) :: l_NL
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: phi

    real, dimension(:,:,:), allocatable, intent(inout) :: RSH

    integer :: l_Nmu
    integer :: l_Nphi
    integer :: NQ

    l_Nmu = size(mu)
    l_Nphi = size(phi)
    NQ = (l_NL+1)**2

    allocate(RSH(l_Nmu,l_Nphi,NQ))

    RSH = q_ind_real_spherical_harmonics(NQ,mu,phi)
end subroutine populate_real_spherical_harmonics

subroutine PN_to_SN(PN, RSH, SN)
    implicit none
    real, dimension(:,:), intent(in) :: PN
    real, dimension(:,:,:), intent(in) :: RSH

    real, dimension(:,:,:), intent(inout) :: SN

    integer :: i, j, q
    integer :: l_Nmu
    integer :: l_Nphi
    integer :: NQ

    l_Nmu = size(RSH,1)
    l_Nphi = size(RSH,2)
    NQ = size(RSH,3)

    do j = 1, l_Nphi
        do i = 1, l_Nmu
            SN(1:NK,i,j) = matmul(PN(1:NK,1:NQ),RSH(i,j,1:NQ))
        end do
    end do
end subroutine PN_to_SN

subroutine SN_to_PN(SN, wL, RSH, PN)
    implicit none
    real, dimension(:,:,:), intent(in) :: SN
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: RSH

    real, dimension(:,:), intent(inout) :: PN

    integer :: i, j, q
    integer :: l_Nmu
    integer :: l_Nphi
    integer :: NQ

    l_Nmu = size(RSH,1)
    l_Nphi = size(RSH,2)
    NQ = size(RSH,3)

    if (size(wL) .ne. l_Nmu) then
        print *, "SN_to_PN problem, stopping, message wip"
        stop
    end if

    do q = 1, NQ
        do i = 1, l_Nmu
            PN(1:NK,q) = PN(1:NK,q) + twopi*wL(i)*matmul(SN(1:NK,i,1:l_Nphi),RSH(i,1:l_Nphi,q))/l_Nphi
        end do
    end do
end subroutine SN_to_PN

real function AL_polys_half_integrals(m,k,l) ! COMES WITH FACTOR OF pi FOR m =/= 0 AND 2pi FOR m = 0
    implicit none
    integer, intent(in) :: m
    integer, intent(in) :: k
    integer, intent(in) :: l

    integer :: j, mp, mpp
    real :: integral

    if (even(l+k)) then
        print *, &
        "6_angular_discretization.f08: function AL_polys_half_integrals: Given l and k must not have same parity."
        stop
    end if

    if (m == 0) then
        if (even(l) .and. odd(k)) then
            integral = -(-1)**((l+k+1)/2)*factorial(l)*factorial(k)/&
                (2**(l+k-1)*(l-k)*(l+k+1)*(factorial(l/2)**2)*(factorial((k-1)/2)**2))
        else
            integral = -(-1)**((k+l+1)/2)*factorial(k)*factorial(l)/&
                (2**(k+l-1)*(k-l)*(k+l+1)*(factorial(k/2)**2)*(factorial((l-1)/2)**2))
        end if
        integral = 2.*pi*integral
    else if (m > 0) then
        integral = 0.
        do j = 0, m
            do mp = m, k
                do mpp = m, l
                    integral = integral + &
                        factorial(mp)*factorial(mpp)*binomial_int(k,mp)*binomial_int(l,mpp)*&
                        gen_binom(.5*(k+mp-1),k)*gen_binom(.5*(l+mpp-1),l)*binomial_int(m,j)*&
                        (-1)**(3*j+mp+mpp-2*m)/((2*j+mp+mpp-2*m+1)*factorial(mp-m)*factorial(mpp-m))
                end do
            end do
        end do
        integral = pi*integral*(2**(l+k))
    else
        integral = 0.
        do j = 0, abs(m)
            do mp = abs(m), k
                do mpp = abs(m), l
                    integral = integral + &
                        factorial(mp)*factorial(mpp)*binomial_int(k,mp)*binomial_int(l,mpp)*&
                        gen_binom(.5*(k+mp-1),k)*gen_binom(.5*(l+mpp-1),l)*binomial_int(abs(m),j)*&
                        (-1)**(3*j+mp+mpp-2*abs(m))/((2*j+mp+mpp-2*abs(m)+1)*factorial(mp-abs(m))*factorial(mpp-abs(m)))
                end do
            end do
        end do
        integral = pi*integral*&
            2**(l+k)*factorial(k-abs(m))*factorial(l-abs(m))/(factorial(k+abs(m))*factorial(l+abs(m)))
    end if

    AL_polys_half_integrals = integral
end function AL_polys_half_integrals

end module angular_discretization