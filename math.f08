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

module math
implicit none
real, parameter :: pi = 4.d0*atan(1.d0)
real, parameter :: twopi = 2*pi
real, parameter :: fourpi = 4*pi

! Public variables for functions
real :: MATHA, MATHB

contains
integer function Levi_Civita(i,j,k)
    ! i, j, k are not to preceed 1 nor exceed 3
    implicit none
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: k

    Levi_Civita = (i-j)*(j-k)*(k-i)/2
end function Levi_Civita

logical function even(m)
    implicit none
    integer, intent(in) :: m

    even = mod(m,2) .eq. 0
end function even

logical function odd(m)
    implicit none
    integer, intent(in) :: m

    odd = .not. even(m)
end function odd

integer function Kron_delta(i,j)
    implicit none
    integer, intent(in) :: i
    integer, intent(in) :: j

    Kron_delta = merge(1,0,i .eq. j)
end function Kron_delta

real function factorial(n)
    implicit none
    integer, intent(in) :: n

    factorial = gamma(real(n+1))
end function factorial

real function double_factorial(n)
    implicit none
    integer, intent(in) :: n

    integer :: sig

    integer :: k

    double_factorial = 1.0
    if (n .gt. 1) then
        sig = n/2-(n-1)/2
        double_factorial = sig*2**(n/2)*factorial(n/2) + &
        (1-sig)*factorial(n+1)/(2**((n+1)/2)*factorial((n+1)/2))
    end if
end function double_factorial

real function gen_binom(n,k)
    implicit none
    real, intent(in) :: n
    integer, intent(in) :: k

    if (n-floor(n) .eq. 0 .and. k .gt. n) then
        gen_binom = 0
    else
        gen_binom = gamma(real(n+1))/(gamma(real(k+1))*gamma(real(n-k+1)))
    end if
end function gen_binom

real function binomial_int(n,k) ! Inputs are integers. Outputs are single-precision real.
    ! ( n )
    ! ( k )
    ! n .ge. k
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: k

    integer :: i

    binomial_int = gamma(real(n+1))/(gamma(real(k+1))*gamma(real(n-k+1)))
end function binomial_int

real function dble_binomial_int(n,k) ! Inputs are integers. Outputs are double-precision real.
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: k

    dble_binomial_int = gamma(dble(n+1))/(gamma(dble(k+1))*gamma(dble(n-k+1)))
end function dble_binomial_int

!real function binomial_real(n,k)
!    implicit none
!    real, intent(in) :: n
!    real, intent(in) :: k
!
!    binomial_real = gamma(n+1)/(gamma(k+1)*gamma(n-k+1))
!end function binomial_real

function Legendre_polynomials(L,x) ! Should make this one from 1:L+1 rather than 0:L
    ! Generates Legendre polynomials up to order L evaluated at
    ! x-values in array x
    implicit none
    integer, intent(in) :: L
    real, dimension(:), intent(in) :: x

    real, dimension(:,:), allocatable :: Legendre_polynomials ! (i,l)

    integer :: k

    allocate(Legendre_polynomials(size(x),0:L))

    Legendre_polynomials(:,0) = 1
    if (L .gt. 0) then
        Legendre_polynomials(:,1) = x
        if (L .gt. 1) then
            do k = 2, L
                Legendre_polynomials(:,k) = (real(2*k-1)/k)*x*Legendre_polynomials(:,k-1) - &
                (real(k-1)/k)*Legendre_polynomials(:,k-2)
            end do
        end if
    end if
end function Legendre_polynomials

function single_Legendre_polynomial(L,x)
    ! Generates Legendre polynomials of order L evaluated at x
    implicit none
    integer, intent(in) :: L
    real, intent(in) :: x

    real, dimension(:), allocatable :: single_Legendre_polynomial ! (l)

    integer :: k

    allocate(single_Legendre_polynomial(0:L))

    single_Legendre_polynomial(0) = 1
    if (L .gt. 0) then
        single_Legendre_polynomial(1) = x
        if (L .gt. 1) then
            do k = 2, L
                single_Legendre_polynomial(k) = (real(2*k-1)/k)*x*single_Legendre_polynomial(k-1) - &
                (real(k-1)/k)*single_Legendre_polynomial(k-2)
            end do
        end if
    end if
end function single_Legendre_polynomial

function Legendre_polynomial_coeffs(N)
    ! Generates coefficients of x^l in Legendre polynomials up to order N
    implicit none
    integer, intent(in) :: N

    real, dimension(:,:), allocatable :: Legendre_polynomial_coeffs

    real, dimension(:,:), allocatable :: coeff ! (l,lp)
    integer :: l, k

    allocate(coeff(0:N,0:N))

    coeff = 0
    do l = 0, N
        do k = 0, l/2 ! Note: technically gives floor of l/2
            coeff(l,l-2*k) = (-1)**(k)*binomial_int(l,k)*binomial_int(2*l-2*k,l)/(2**l)
            if (coeff(l,l-2*k) .eq. coeff(l,l-2*k) + 1.0) then
                coeff(l,l-2*k) = (-1)**(k)*&
                real(dble_binomial_int(l,k)*dble_binomial_int(2*l-2*k,l))/(2**l)
                ! NOTE: This becomes relevant for N .ge. 18, maybe even 17 and 16.
            end if
        end do
    end do

    allocate(Legendre_polynomial_coeffs(0:N,0:N))
    Legendre_polynomial_coeffs = coeff
end function Legendre_polynomial_coeffs

function assoc_Legendre_polynomials(L,x)
    ! Generates associated Legendre polynomials up to order L evaluated at
    ! x-values in array x
    ! Importantly, m .ge. 0 always. Use full_assoc_Legendre_polynomials for m < 0
    implicit none
    integer, intent(in) :: L
    real, dimension(:), intent(in) :: x

    real, dimension(:,:,:), allocatable :: assoc_Legendre_polynomials ! (i,m,l)

    integer :: m, k, i, ell
    integer :: N
    real, dimension(:,:,:), allocatable :: P

    N = size(x)

    allocate(assoc_Legendre_polynomials(N,L+1,L+1))
    allocate(P(N,0:L,0:L))

    P = 0.0

    do m = 1, L
        P(:,m,m) = (-1)**(m)*double_factorial(2*m-1)*sqrt(1-x**2)**(m)
        if (m .eq. L) exit
        P(:,m,m+1) = x*(2*m+1)*P(:,m,m)
        do k = 1, L-m-1
            P(:,m,m+k+1) = (x*(2*(m+k+1)-1)*P(:,m,m+k)-(2*m+k)*P(:,m,m+k-1))/(k+1)
        end do
    end do

    P(:,0,:) = Legendre_polynomials(L,x)

    assoc_Legendre_polynomials(1:N,1:L+1,1:L+1) = P(1:N,0:L,0:L)
end function assoc_Legendre_polynomials

function full_assoc_Legendre_polynomials(L,x)
    ! Generates associated Legendre polynomials up to order L evaluated at
    ! x-values in array x
    implicit none
    integer, intent(in) :: L
    real, dimension(:), intent(in) :: x

    real, dimension(:,:,:), allocatable :: full_assoc_Legendre_polynomials ! (i,m,l)

    integer :: N
    real, dimension(:,:,:), allocatable :: P
    integer :: m, k, i

    N = size(x)

    allocate(P(N,0:L,0:L))

    P = 0.
    do m = 1, L
        P(:,m,m) = (-1)**(m)*double_factorial(2*m-1)*sqrt(1-x**2)**(m)
        if (m .eq. L) exit
        P(:,m,m+1) = x*(2*m+1)*P(:,m,m)
        do k = 1, L-m-1
            P(:,m,m+k+1) = (x*(2*(m+k+1)-1)*P(:,m,m+k)-(2*m+k)*P(:,m,m+k-1))/(k+1)
        end do
    end do
    P(:,0,:) = Legendre_polynomials(L,x)

    allocate(full_assoc_Legendre_polynomials(N,-L:L,0:L))

    full_assoc_Legendre_polynomials(:,0:L,0:L) = P(:,0:L,0:L)
    do k = 0, L
        do m = -k, -1
            full_assoc_Legendre_polynomials(:,m,k) = (-1)**(abs(m))*factorial(k-abs(m))*P(:,abs(m),k)/factorial(k+abs(m))
        end do
    end do
end function full_assoc_Legendre_polynomials

function q_ind_associated_Legendre_polynomials(NQ,x)
    ! Generates associated Legendre polynomials up to order NQ evaluated at
    ! x-values in array x
    ! Indexed by q = l*(l+1)+m+1
    implicit none
    integer, intent(in) :: NQ
    real, dimension(:), intent(in) :: x

    real, dimension(:,:), allocatable :: q_ind_associated_Legendre_polynomials ! (i,q)

    integer :: m, k, i
    integer :: N
    integer :: L
    real, dimension(:,:,:), allocatable :: P ! (i,m,l)

    N = size(x)
    L = ceiling(sqrt(real(NQ)))-1

    allocate(P(N,0:L,0:L))

    P = 0.0
    do m = 1, L
        P(:,m,m) = (-1)**(m)*double_factorial(2*m-1)*sqrt(1-x**2)**(m)
        if (m .eq. L) exit
        P(:,m,m+1) = x*(2*m+1)*P(:,m,m)
        do k = 1, L-m-1
            P(:,m,m+k+1) = (x*(2*(m+k+1)-1)*P(:,m,m+k)-(2*m+k)*P(:,m,m+k-1))/(k+1)
        end do
    end do
    P(1:N,0,0:L) = Legendre_polynomials(L,x)

    allocate(q_ind_associated_Legendre_polynomials(N,NQ))

    do k = 0, L
        do m = 0, k
            q_ind_associated_Legendre_polynomials(1:N,k*(k+1)+m+1) = P(1:N,m,k)
        end do
    end do
    do k = 0, L
        do m = -k, -1
            q_ind_associated_Legendre_polynomials(1:N,k*(k+1)+m+1) = &
            (-1)**(abs(m))*factorial(k-abs(m))*P(1:N,abs(m),k)/factorial(k+abs(m))
        end do
    end do
end function q_ind_associated_Legendre_polynomials

function real_spherical_harmonics(L,muvals,phivals)
    ! Gives positive RSH up to order L with positive m evaluated at
    ! mu-values in array muvals and
    ! phi-values in array phivals
    implicit none
    integer, intent(in) :: L
    real, dimension(:), intent(in) :: muvals
    real, dimension(:), intent(in) :: phivals

    real, dimension(:,:,:,:), allocatable :: real_spherical_harmonics

    integer :: m, ell, j
    integer :: Nmu
    integer :: Nphi
    real, dimension(:,:,:), allocatable :: Legendre
    real, dimension(:,:), allocatable :: Sm
    real, dimension(:,:), allocatable :: norms

    Nmu = size(muvals)
    Nphi = size(phivals)

    allocate(real_spherical_harmonics(Nmu,Nphi,0:L,0:L))
    allocate(Legendre(Nmu,0:L,0:L))
    allocate(Sm(Nphi,0:L))
    allocate(norms(0:L,0:L))

    Legendre = assoc_Legendre_polynomials(L,muvals)
    Sm(1:Nphi,0) = 1
    do m = 1, L
        Sm(1:Nphi,m) = cos(m*phivals)
    end do
    do ell = 0, L
        do m = 0, ell
            norms(m,ell) = RSH_norm(ell,m)
        end do
    end do

    do ell = 0, L
        do m = 0, ell
            do j = 1, Nphi
                real_spherical_harmonics(1:Nmu,j,m,ell) = &
                norms(m,ell)*Sm(j,m)*Legendre(1:Nmu,m,ell)
            end do
        end do
    end do
end function real_spherical_harmonics

function q_ind_real_spherical_harmonics(NQ,muvals,phivals)
    implicit none
    integer, intent(in) :: NQ
    real, dimension(:), intent(in) :: muvals
    real, dimension(:), intent(in) :: phivals

    real, dimension(:,:,:), allocatable :: q_ind_real_spherical_harmonics ! (i,j,q)

    integer :: m, ell, j, q
    integer :: L
    integer :: Nmu
    integer :: Nphi
    real, dimension(:,:,:), allocatable :: Legendre
    real, dimension(:,:), allocatable :: Sm
    real, dimension(:,:), allocatable :: norms

    L = ceiling(sqrt(real(NQ)))-1
    Nmu = size(muvals)
    Nphi = size(phivals)

    allocate(q_ind_real_spherical_harmonics(Nmu,Nphi,NQ))
    allocate(Legendre(Nmu,-L:L,0:L))
    allocate(Sm(Nphi,-L:L))
    allocate(norms(-L:L,0:L))

    Legendre = full_assoc_Legendre_polynomials(L,muvals)
    Sm(1:Nphi,0) = 1
    do m = 1, L
        Sm(1:Nphi,m) = cos(m*phivals)
        Sm(1:Nphi,-m) = sin(m*phivals)
    end do

    do ell = 0, L
        do m = -ell, ell
            norms(m,ell) = RSH_norm(ell,m)
        end do
    end do

    do ell = 0, L
        do m = -ell, ell
            do j = 1, Nphi
                q_ind_real_spherical_harmonics(1:Nmu,j,ell*(ell+1)+m+1) = &
                norms(m,ell)*Sm(j,m)*Legendre(1:Nmu,m,ell)
            end do
        end do
    end do
end function q_ind_real_spherical_harmonics

function q_ind_real_spherical_harmonics_vector_input(NQ,vectors)
    implicit none
    integer, intent(in) :: NQ
    real, dimension(:,:), intent(in) :: vectors ! (dir,k) ! MUST BE UNIT VECTORS

    real, dimension(:,:), allocatable :: q_ind_real_spherical_harmonics_vector_input ! (k,q)

    integer :: m, ell, k, q
    integer :: L
    integer :: N
    real, dimension(:), allocatable :: muvals
    real, dimension(:), allocatable :: phivals
    real, dimension(:,:,:), allocatable :: Legendre
    real, dimension(:,:), allocatable :: Sm
    real, dimension(:,:), allocatable :: norms

    L = ceiling(sqrt(real(NQ)))-1
    N = size(vectors,2)

    allocate(q_ind_real_spherical_harmonics_vector_input(N,NQ))
    allocate(Legendre(N,-L:L,0:L))
    allocate(Sm(N,-L:L))
    allocate(norms(-L:L,0:L))
    allocate(muvals(N))
    allocate(phivals(N))

    muvals = vectors(3,1:N)
    phivals = atan2(vectors(2,1:N),vectors(1,1:N))

    Legendre = full_assoc_Legendre_polynomials(L,muvals)
    Sm(1:N,0) = 1
    do m = 1, L
        Sm(1:N,m) = cos(m*phivals)
        Sm(1:N,-m) = sin(m*phivals)
    end do

    do ell = 0, L
        do m = -ell, ell
            norms(m,ell) = RSH_norm(ell,m)
        end do
    end do

    do ell = 0, L
        do m = -ell, ell
            do k = 1, N
                q_ind_real_spherical_harmonics_vector_input(k,ell*(ell+1)+m+1) = &
                norms(m,ell)*Sm(k,m)*Legendre(k,m,ell)
            end do
        end do
    end do
end function q_ind_real_spherical_harmonics_vector_input

real function Wigner_d_matrix_elements(l,mp,m,b)
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: mp
    integer, intent(in) :: m
    real, intent(in) :: b

    integer :: s
    integer :: smin
    integer :: smax
    real :: LC
    real :: cb
    real :: sb
    real :: el

    !if (abs(mp) .gt. l .or. abs(m) .gt. l) then
    !    print *, "PROBLEM"
    !    stop
    !end if

    LC = sqrt(factorial(l+mp)*factorial(l-mp)*factorial(l+m)*factorial(l-m))
    smin = max(0,m-mp)
    smax = min(l+m,l-mp)
    cb = cos(0.5*b)
    sb = sin(0.5*b)

    el = 0
    do s = smin, smax
        el = el + &
            (-1)**(mp-m+s)*(cb**(2*l+m-mp-2*s))*(sb**(mp-m+2*s))/&
            real(factorial(l+m-s)*factorial(s)*factorial(mp-m+s)*factorial(l-mp-s))
    end do
    Wigner_d_matrix_elements = LC*el
end function Wigner_d_matrix_elements

real function Wigner_RD_matrix_elements(l,mp,m,a,b,c)
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: mp
    integer, intent(in) :: m
    real, intent(in) :: a
    real, intent(in) :: b
    real, intent(in) :: c

    real :: el

    if (m .gt. 0) then
        if (mp .gt. 0) then
            el = (-1)**(mp-m)*Wigner_d_matrix_elements(l,mp,m,b)*cos(mp*a+m*c) + &
                 (-1)**(mp)*Wigner_d_matrix_elements(l,mp,-m,b)*cos(mp*a-m*c)
        else if (mp .lt. 0) then
            el = -(-1)**(m)*Wigner_d_matrix_elements(l,mp,m,b)*sin(mp*a+m*c) - &
                 Wigner_d_matrix_elements(l,mp,-m,b)*sin(mp*a-m*c)
        else
            el = (-1)**(m)*sqrt(2.0)*Wigner_d_matrix_elements(l,0,m,b)*cos(m*c)
        end if
    else if (m .lt. 0) then
        if (mp .gt. 0) then
            el = (-1)**(m)*Wigner_d_matrix_elements(l,mp,m,b)*sin(mp*a+m*c) - &
                 (-1)**(mp+m)*Wigner_d_matrix_elements(l,mp,-m,b)*sin(mp*a-m*c)
        else if (mp .lt. 0) then
            el = Wigner_d_matrix_elements(l,mp,m,b)*cos(mp*a+m*c) - &
                 (-1)**(m)*Wigner_d_matrix_elements(l,mp,-m,b)*cos(mp*a-m*c)
        else
            el = sqrt(2.0)*Wigner_d_matrix_elements(l,0,m,b)*sin(m*c)
        end if
    else
        if (mp .gt. 0) then
            el = (-1)**(mp)*sqrt(2.0)*Wigner_d_matrix_elements(l,mp,0,b)*cos(mp*a)
        else if (mp .lt. 0) then
            el = -sqrt(2.0)*Wigner_d_matrix_elements(l,mp,0,b)*sin(mp*a)
        else
            el = Wigner_d_matrix_elements(l,0,0,b)
        end if
    end if

    Wigner_RD_matrix_elements = el
end function Wigner_RD_matrix_elements

real function SH_norm(l,m)
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: m

    !if (abs(m) .gt. l) print *, "PROBLEM"

    if (abs(m) .gt. l) then
        SH_norm = 0.0
    else
        SH_norm = sqrt((2*l+1)*factorial(l-m)/(fourpi*factorial(l+m)))
    end if
end function SH_norm

real function RSH_norm(l,m)
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: m

    RSH_norm = merge((-1)**(m),1,m .ge. 0)*&
        sqrt((2*l+1)*factorial(l-m)/(twopi*(1+Kron_delta(m,0))*factorial(l+m)))
end function RSH_norm

function RSH_kappas(l,m,dir)
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: m
    integer, intent(in) :: dir

    real, dimension(:), allocatable :: RSH_kappas

    real :: C1
    real :: C2
    real :: A
    real :: B

    C1 = RSH_norm(l,m)

    if (dir .eq. 1) then
        ! 1: (l-1,m+1)
        ! 2: (l+1,m+1)
        ! 3: (l+1,m-1)
        ! 4: (l-1,m-1)

        allocate(RSH_kappas(4))
        A = 0.5*(1+Kron_delta(m,0)-Kron_delta(m,-1))
        B = 0.5*(1-Kron_delta(m,0))

        if (l-1 .lt. 0 .or. abs(m+1) .gt. l-1) then
            RSH_kappas(1) = 0
        else
            C2 = RSH_norm(l-1,m+1)
            RSH_kappas(1) = C1*A/((2*l+1)*C2)
        end if

        if (l+1 .lt. 0 .or. abs(m+1) .gt. l+1) then
            RSH_kappas(2) = 0
        else
            C2 = RSH_norm(l+1,m+1)
            RSH_kappas(2) = -C1*A/((2*l+1)*C2)
        end if

        if (l+1 .lt. 0 .or. abs(m-1) .gt. l+1) then
            RSH_kappas(3) = 0
        else
            C2 = RSH_norm(l+1,m-1)
            RSH_kappas(3) = C1*B*(l-m+1)*(l-m+2)/((2*l+1)*C2)
        end if

        if (l-1 .lt. 0 .or. abs(m-1) .gt. l-1) then
            RSH_kappas(4) = 0
        else
            C2 = RSH_norm(l-1,m-1)
            RSH_kappas(4) = -C1*B*(l+m-1)*(l+m)/((2*l+1)*C2)
        end if
    else if (dir .eq. 2) then
        ! 1: (l-1,-(m+1))
        ! 2: (l+1,-(m+1))
        ! 3: (l+1,-(m-1))
        ! 4: (l-1,-(m-1))

        allocate(RSH_kappas(4))
        A = 0.5*(1+Kron_delta(m,0))
        B = -0.5*(1-Kron_delta(m,0)-Kron_delta(m,1))

        if (l-1 .lt. 0 .or. abs(-(m+1)) .gt. l-1) then
            RSH_kappas(1) = 0
        else
            C2 = RSH_norm(l-1,-(m+1))
            RSH_kappas(1) = (-1)**(m-1)*C1*A*factorial(l+m)/&
                ((2*l+1)*C2*factorial(l-m-2))
        end if

        if (l+1 .lt. 0 .or. abs(-(m+1)) .gt. l+1) then
            RSH_kappas(2) = 0
        else
            C2 = RSH_norm(l+1,-(m+1))
            RSH_kappas(2) = -(-1)**(m-1)*C1*A*factorial(l+m+2)/&
                ((2*l+1)*C2*factorial(l-m))
        end if

        if (l+1 .lt. 0 .or. abs(-(m-1)) .gt. l+1) then
            RSH_kappas(3) = 0
        else
            C2 = RSH_norm(l+1,-(m-1))
            RSH_kappas(3) = (-1)**(m-1)*C1*B*(l-m+1)*(l-m+2)*factorial(l+m)/&
                ((2*l+1)*C2*factorial(l-m+2))
        end if

        if (l-1 .lt. 0 .or. abs(-(m-1)) .gt. l-1) then
            RSH_kappas(4) = 0
        else
            C2 = RSH_norm(l-1,-(m-1))
            RSH_kappas(4) = -(-1)**(m-1)*C1*B*(l+m-1)*(l+m)*factorial(l+m-2)/&
                ((2*l+1)*C2*factorial(l-m))
        end if
    else if (dir .eq. 3) then
        ! 1: (l+1,m)
        ! 2: (l-1,m)

        allocate(RSH_kappas(2))

        if (l+1 .lt. 0 .or. abs(m) .gt. l+1) then
            RSH_kappas(1) = 0
        else
            C2 = RSH_norm(l+1,m)
            RSH_kappas(1) = C1*(l-m+1)/((2*l+1)*C2)
        end if

        if (l-1 .lt. 0 .or. abs(m) .gt. l-1) then
            RSH_kappas(2) = 0
        else
            C2 = RSH_norm(l-1,m)
            RSH_kappas(2) = C1*(l+m)/((2*l+1)*C2)
        end if
    else
        print *, "math.f08: function RSH_kappas(l,m,dir): Given direction is outside of [1,3]."
        stop
    end if
end function RSH_kappas

function RSH_kappas_indices(l,m,dir)
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: m
    integer, intent(in) :: dir

    integer, dimension(:), allocatable :: RSH_kappas_indices

    if (abs(m) .gt. l) then
        print *, "m greater than l"
        stop
    end if

    if (dir .eq. 1) then
        ! 1: (l-1,m+1)
        ! 2: (l+1,m+1)
        ! 3: (l+1,m-1)
        ! 4: (l-1,m-1)

        allocate(RSH_kappas_indices(4))

        if (l-1 .lt. 0 .or. abs(m+1) .gt. l-1) then
            RSH_kappas_indices(1) = 0
        else
            RSH_kappas_indices(1) = (l-1)*l+m+2
        end if

        if (l+1 .lt. 0 .or. abs(m+1) .gt. l+1) then
            RSH_kappas_indices(2) = 0
        else
            RSH_kappas_indices(2) = (l+1)*(l+2)+m+2
        end if

        if (l+1 .lt. 0 .or. abs(m-1) .gt. l+1) then
            RSH_kappas_indices(3) = 0
        else
            RSH_kappas_indices(3) = (l+1)*(l+2)+m
        end if

        if (l-1 .lt. 0 .or. abs(m-1) .gt. l-1) then
            RSH_kappas_indices(4) = 0
        else
            RSH_kappas_indices(4) = (l-1)*l+m
        end if
    else if (dir .eq. 2) then
        ! 1: (l-1,-(m+1))
        ! 2: (l+1,-(m+1))
        ! 3: (l+1,-(m-1))
        ! 4: (l-1,-(m-1))

        allocate(RSH_kappas_indices(4))

        if (l-1 .lt. 0 .or. abs(-(m+1)) .gt. l-1) then
            RSH_kappas_indices(1) = 0
        else
            RSH_kappas_indices(1) = (l-1)*l-m
        end if

        if (l+1 .lt. 0 .or. abs(-(m+1)) .gt. l+1) then
            RSH_kappas_indices(2) = 0
        else
            RSH_kappas_indices(2) = (l+1)*(l+2)-m
        end if

        if (l+1 .lt. 0 .or. abs(-(m-1)) .gt. l+1) then
            RSH_kappas_indices(3) = 0
        else
            RSH_kappas_indices(3) = (l+1)*(l+2)-m+2
        end if

        if (l-1 .lt. 0 .or. abs(-(m-1)) .gt. l-1) then
            RSH_kappas_indices(4) = 0
        else
            RSH_kappas_indices(4) = (l-1)*l-m+2
        end if
    else if (dir .eq. 3) then
        ! 1: (l+1,m)
        ! 2: (l-1,m)

        allocate(RSH_kappas_indices(2))

        if (l+1 .lt. 0 .or. abs(m) .gt. l+1) then
            RSH_kappas_indices(1) = 0
        else
            RSH_kappas_indices(1) = (l+1)*(l+2)+m+1
        end if

        if (l-1 .lt. 0 .or. abs(m) .gt. l-1) then
            RSH_kappas_indices(2) = 0
        else
            RSH_kappas_indices(2) = (l-1)*l+m+1
        end if
    else
        print *, "math.f08: function RSH_kappas_indices(l,m,dir): Given dir exceeds [1,3]."
        stop
    end if
end function RSH_kappas_indices

function Euler_zyz_vector_to_north_pole_angles(n)
    ! explanation
    ! (1) is beta, (2) is gamma, alpha is assumed to be zero
    implicit none
    real, dimension(3), intent(in) :: n

    real, dimension(2) :: Euler_zyz_vector_to_north_pole_angles

    real, dimension(2) :: angles

    angles(1) = -acos(n(3))
    angles(2) = -atan2(n(2),n(1)) ! Double check both of these, also determine specific cases.
end function Euler_zyz_vector_to_north_pole_angles

pure real function Euclid(u)
    ! Calculates the Euclidean norm squared of a given vector u
    ! u: array containing the coordinates of the point in some ND Euclidean space
    implicit none
    real, dimension(:), intent(in) :: u

    Euclid = dot_product(u,u) ! norm2(u)**2
    ! Use norm2 squared??
end function Euclid

real function interp2d(x,y,z,x0,y0)
    implicit none
    real, dimension(:), intent(in) :: x
    real, dimension(:), intent(in) :: y
    real, dimension(:,:), intent(in) :: z
    real, intent(in) :: x0
    real, intent(in) :: y0

    integer :: xstart
    integer :: ystart

    if (x0 .eq. x(1)) xstart = 1
    if (y0 .eq. y(1)) ystart = 1
    if (x0 .eq. x(size(x))) xstart = size(x)
    if (y0 .eq. y(size(y))) ystart = size(y)
    xstart = findloc(x0 .ge. x(1:size(x)-1) .and. x0 .lt. x(2:size(x)), .true., dim = 1)
    ystart = findloc(y0 .ge. y(1:size(x)-1) .and. y0 .lt. y(2:size(y)), .true., dim = 1)
    if (xstart .eq. 0 .or. ystart .eq. 0) then
        print *, "function interp2d: input variable is not within domain."
        print *, "x0:", x0
        print *, "Min of x:", minval(x)
        print *, "Max of x:", maxval(x)
        print *, "y0:", y0
        print *, "Min of y:", minval(y)
        print *, "Max of y:", maxval(y)
        stop
    end if

    interp2d = (z(xstart,ystart)*(x(xstart+1)-x0)*(y(ystart+1)-y0) + &
                z(xstart+1,ystart)*(x0-x(xstart))*(y(ystart+1)-y0) + &
                z(xstart,ystart+1)*(x(xstart+1)-x0)*(y0-y(ystart)) + &
                z(xstart+1,ystart+1)*(x0-x(xstart))*(y0-y(ystart)))/ &
                ((x(xstart+1)-x(xstart))*(y(ystart+1)-y(ystart)))
end function interp2d

pure integer function Heaviside(x)
    ! Follows convention:
    ! x .gt. 0: Heaviside = 1
    ! x .le. 0: Heaviside = 0
    ! For x .ge. 0 ~ 1, replace floor with ceiling
    implicit none
    real, intent(in) :: x

    Heaviside = int(floor(0.5+0.5*sign(1.0,x)))
end function Heaviside

pure integer function rect(x,a,b)
    ! Follows convention:
    !
    implicit none
    real, intent(in) :: x
    real, intent(in) :: a
    real, intent(in) :: b

    rect = int(ceiling(0.5*sign(1.0,x-a)+0.5)-floor(0.5*sign(1.0,x-b)+0.5))
end function rect

pure function diag_real(A)
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, dimension(:), allocatable :: diag_real

    integer :: N
    integer :: i

    N = size(A,1)

    allocate(diag_real(N))

    forall(i = 1:N)
        diag_real(i) = A(i,i)
    end forall
end function diag_real

pure function diag_int(A)
    implicit none
    integer, dimension(:,:), intent(in) :: A

    integer, dimension(:), allocatable :: diag_int

    integer :: N
    integer :: i

    N = size(A,1)

    allocate(diag_int(N))

    forall(i = 1:N)
        diag_int(i) = A(i,i)
    end forall
end function diag_int

pure real function trace_real(A)
    implicit none
    real, dimension(:,:), intent(in) :: A

    trace_real = sum(diag_real(A))
end function trace_real

pure integer function trace_int(A)
    implicit none
    integer, dimension(:,:), intent(in) :: A

    trace_int = sum(diag_int(A))
end function trace_int

pure function submatrices_all(A)
    ! Gives all of the submatrices of A, the (i,j) is formed by removing the ith row and jth column
    ! Returns submatrices_all(i,j,m,k), where i and j index the submatrix of A, while
    ! m and k index the elements of the submatrix.
    implicit none
    real, dimension(:,:), intent(in) :: A

    integer :: N
    real, dimension(:,:,:,:), allocatable :: submatrices_all

    integer :: i, j

    N = size(A,1)

    allocate(submatrices_all(N,N,N-1,N-1))

    forall(i = 1:size(A,1), j = 1:size(A,1)) ! REVISIT THIS CONSTRUCT
        submatrices_all(i,j,1:i-1,1:j-1) = A(1:i-1,1:j-1)
        submatrices_all(i,j,1:i-1,j:size(A,1)-1) = A(1:i-1,j+1:size(A,1))
        submatrices_all(i,j,i:size(A,1)-1,1:j-1) = A(i+1:size(A,1),1:j-1)
        submatrices_all(i,j,i:size(A,1)-1,j:size(A,1)-1) = A(i+1:size(A,1),j+1:size(A,1))
    end forall
end function submatrices_all

pure function submatrices_one(A,i,j)
    ! Gives the submatrix of A formed by removing the ith row and jth column
    ! Returns submatrices_one(m,k), where m and k index the elements of the submatrix.
    implicit none
    real, dimension(:,:), intent(in) :: A

    integer, intent(in) :: i
    integer, intent(in) :: j

    integer :: N
    real, dimension(:,:), allocatable :: submatrices_one

    N = size(A,1)

    allocate(submatrices_one(N-1,N-1))

    submatrices_one(1:i-1,1:j-1) = A(1:i-1,1:j-1)
    submatrices_one(1:i-1,j:size(A,1)-1) = A(1:i-1,j+1:size(A,1))
    submatrices_one(i:size(A,1)-1,1:j-1) = A(i+1:size(A,1),1:j-1)
    submatrices_one(i:size(A,1)-1,j:size(A,1)-1) = A(i+1:size(A,1),j+1:size(A,1))
end function submatrices_one

pure real function matdet(A)
    ! Defined only for size(A) = 2, 3, 4
    implicit none
    real, dimension(:,:), intent(in) :: A

    integer :: N
    real, dimension(:,:,:,:), allocatable :: B
    integer :: i

    N = size(A,1)

    if (N .eq. 2) then
        matdet = A(1,1)*A(2,2)-A(1,2)*A(2,1)
    else if (N .eq. 3) then
        matdet = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+&
                 A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+&
                 A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    else if (N .eq. 4) then
        allocate(B(N,N,N-1,N-1))
        B = submatrices_all(A)
        matdet = 0.
        do i = 1, N
            matdet = matdet + ((-1)**(i-1))*A(i,1)*&
            (B(i,1,1,1)*(B(i,1,2,2)*B(i,1,3,3)-B(i,1,2,3)*B(i,1,3,2))+&
             B(i,1,1,2)*(B(i,1,2,3)*B(i,1,3,1)-B(i,1,2,1)*B(i,1,3,3))+&
             B(i,1,1,3)*(B(i,1,2,1)*B(i,1,3,2)-B(i,1,2,2)*B(i,1,3,1)))
        end do
    else
        matdet = 0.
    end if
end function matdet

pure function minors_all(A)
    ! Gives the matrix containing the determinants of the submatrices of A
    ! Returns minor(i,j), where i and j index the determinant of the (i,j) submatrix of A
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, dimension(:,:), allocatable :: minors_all

    integer :: N
    real, dimension(:,:,:,:), allocatable :: submatrices
    integer :: i, j

    N = size(A,1)

    allocate(minors_all(N,N))
    allocate(submatrices(N,N,N-1,N-1))

    submatrices(:,:,:,:) = submatrices_all(A)

    do j = 1, N
        do i = 1, N
            minors_all(i,j) = matdet(submatrices(i,j,:,:))
        end do
    end do
end function minors_all

pure real function minors_one(A,i,j)
    ! Gives the determinant of the (i,j) submatrix of A
    implicit none
    real, dimension(:,:), intent(in) :: A
    integer, intent(in) :: i
    integer, intent(in) :: j

    integer :: N
    real, dimension(:,:), allocatable :: submatrix

    N = size(A,1)

    allocate(submatrix(N-1,N-1))

    submatrix = submatrices_one(A,i,j)

    minors_one = matdet(submatrix)
end function minors_one

pure function cofactor(A) ! CHECK THAT I AM NOT REPEATING ANY OPERATIONS DURING PRACTICAL USE
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, dimension(:,:), allocatable :: cofactor

    integer :: N
    real, dimension(:,:), allocatable :: minors
    integer :: i, j

    N = size(A,1)

    allocate(cofactor(N,N))
    allocate(minors(N,N))

    minors = minors_all(A)

    do j = 1, N
        do i = 1, N
            cofactor(i,j) = ((-1)**(i+j))*minors(i,j)
        end do
    end do
end function cofactor

pure function matinv(A)
    ! Recommend use with N not exceeding 4.
    ! Errors may result from lack of pivoting. (Make pivot inv? Make Gaussian elim.?)
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, dimension(:,:), allocatable :: matinv

    integer :: N

    N = size(A,1)

    allocate(matinv(N,N))

    if (N .eq. 2) then
        matinv(1,1) = A(2,2)
        matinv(1,2) = -A(1,2)
        matinv(2,1) = -A(2,1)
        matinv(2,2) = A(1,1)
        matinv = matinv/(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    else
        matinv = transpose(cofactor(A))/(matdet(A))
    end if
end function matinv

function explicit_2x2_matinv(A)
    implicit none
    real, dimension(2,2), intent(in) :: A

    real, dimension(2,2) :: explicit_2x2_matinv

    real :: detA

    detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    explicit_2x2_matinv(1,1) = A(2,2)/detA
    explicit_2x2_matinv(1,2) = -A(1,2)/detA
    explicit_2x2_matinv(2,1) = -A(2,1)/detA
    explicit_2x2_matinv(2,2) = A(1,1)/detA
end function explicit_2x2_matinv

function explicit_4x4_matinv(A)
    ! Can probably be optimized. Or I probably shouldn't be using it.
    implicit none
    real, dimension(4,4), intent(in) :: A

    real, dimension(4,4) :: explicit_4x4_matinv

    real, dimension(4,4) :: B
    real :: detA
    !real, dimension(4,4) :: ident

    !detA = A(1,1)*A(2,2)*A(3,3)*A(4,4) + A(1,1)*A(2,3)*A(3,4)*A(4,2) + A(1,1)*A(2,4)*A(3,2)*A(4,3) + &
    !       A(1,2)*A(2,1)*A(3,4)*A(4,3) + A(1,2)*A(2,3)*A(3,1)*A(4,4) + A(1,2)*A(2,4)*A(3,3)*A(4,1) + &
    !       A(1,3)*A(2,1)*A(3,2)*A(4,4) + A(1,3)*A(2,2)*A(3,4)*A(4,1) + A(1,3)*A(2,4)*A(3,1)*A(4,2) + &
    !       A(1,4)*A(2,1)*A(3,3)*A(4,2) + A(1,4)*A(2,2)*A(3,1)*A(4,3) + A(1,4)*A(2,3)*A(3,2)*A(4,1) - &
    !       A(1,1)*A(2,2)*A(3,4)*A(4,3) - A(1,1)*A(2,3)*A(3,2)*A(4,4) - A(1,1)*A(2,4)*A(3,3)*A(4,2) - &
    !       A(1,2)*A(2,1)*A(3,3)*A(4,4) - A(1,2)*A(2,3)*A(3,4)*A(4,1) - A(1,2)*A(2,4)*A(3,1)*A(4,3) - &
    !       A(1,3)*A(2,1)*A(3,4)*A(4,2) - A(1,3)*A(2,2)*A(3,1)*A(4,4) - A(1,3)*A(2,4)*A(3,2)*A(4,1) - &
    !       A(1,4)*A(2,1)*A(3,2)*A(4,3) - A(1,4)*A(2,2)*A(3,3)*A(4,1) - A(1,4)*A(2,3)*A(3,1)*A(4,2)

    B(1,1) = A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) - &
             A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4) - A(2,4)*A(3,3)*A(4,2)
    B(2,1) = A(2,1)*A(3,4)*A(4,3) + A(2,3)*A(3,1)*A(4,4) + A(2,4)*A(3,3)*A(4,1) - &
             A(2,1)*A(3,3)*A(4,4) - A(2,3)*A(3,4)*A(4,1) - A(2,4)*A(3,1)*A(4,3)
    B(3,1) = A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) - &
             A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4) - A(2,4)*A(3,2)*A(4,1)
    B(4,1) = A(2,1)*A(3,3)*A(4,2) + A(2,2)*A(3,1)*A(4,3) + A(2,3)*A(3,2)*A(4,1) - &
             A(2,1)*A(3,2)*A(4,3) - A(2,2)*A(3,3)*A(4,1) - A(2,3)*A(3,1)*A(4,2)

    B(1,2) = A(1,2)*A(3,4)*A(4,3) + A(1,3)*A(3,2)*A(4,4) + A(1,4)*A(3,3)*A(4,2) - &
             A(1,2)*A(3,3)*A(4,4) - A(1,3)*A(3,4)*A(4,2) - A(1,4)*A(3,2)*A(4,3)
    B(2,2) = A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,3) - &
             A(1,1)*A(3,4)*A(4,3) - A(1,3)*A(3,1)*A(4,4) - A(1,4)*A(3,3)*A(4,1)
    B(3,2) = A(1,1)*A(3,4)*A(4,2) + A(1,2)*A(3,1)*A(4,4) + A(1,4)*A(3,2)*A(4,1) - &
             A(1,1)*A(3,2)*A(4,4) - A(1,2)*A(3,4)*A(4,1) - A(1,4)*A(3,1)*A(4,2)
    B(4,2) = A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) - &
             A(1,1)*A(3,3)*A(4,2) - A(1,2)*A(3,1)*A(4,3) - A(1,3)*A(3,2)*A(4,1)

    B(1,3) = A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) - &
             A(1,2)*A(2,4)*A(4,3) - A(1,3)*A(2,2)*A(4,4) - A(1,4)*A(2,3)*A(4,2)
    B(2,3) = A(1,1)*A(2,4)*A(4,3) + A(1,3)*A(2,1)*A(4,4) + A(1,4)*A(2,3)*A(4,1) - &
             A(1,1)*A(2,3)*A(4,4) - A(1,3)*A(2,4)*A(4,1) - A(1,4)*A(2,1)*A(4,3)
    B(3,3) = A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) - &
             A(1,1)*A(2,4)*A(4,2) - A(1,2)*A(2,1)*A(4,4) - A(1,4)*A(2,2)*A(4,1)
    B(4,3) = A(1,1)*A(2,3)*A(4,2) + A(1,2)*A(2,1)*A(4,3) + A(1,3)*A(2,2)*A(4,1) - &
             A(1,1)*A(2,2)*A(4,3) - A(1,2)*A(2,3)*A(4,1) - A(1,3)*A(2,1)*A(4,2)

    B(1,4) = A(1,2)*A(2,4)*A(3,3) + A(1,3)*A(2,2)*A(3,4) + A(1,4)*A(2,3)*A(3,2) - &
             A(1,2)*A(2,3)*A(3,4) - A(1,3)*A(2,4)*A(3,2) - A(1,4)*A(2,2)*A(3,3)
    B(2,4) = A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) - &
             A(1,1)*A(2,4)*A(3,3) - A(1,3)*A(2,1)*A(3,4) - A(1,4)*A(2,3)*A(3,1)
    B(3,4) = A(1,1)*A(2,4)*A(3,2) + A(1,2)*A(2,1)*A(3,4) + A(1,4)*A(2,2)*A(3,1) - &
             A(1,1)*A(2,2)*A(3,4) - A(1,2)*A(2,4)*A(3,1) - A(1,4)*A(2,1)*A(3,2)
    B(4,4) = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - &
             A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)

    detA = A(1,1)*B(1,1) + A(2,1)*B(1,2) + A(3,1)*B(1,3) + A(4,1)*B(1,4) + &
           A(1,2)*B(2,1) + A(2,2)*B(2,2) + A(3,2)*B(2,3) + A(4,2)*B(2,4) + &
           A(1,3)*B(3,1) + A(2,3)*B(3,2) + A(3,3)*B(3,3) + A(4,3)*B(3,4) + &
           A(1,4)*B(4,1) + A(2,4)*B(4,2) + A(3,4)*B(4,3) + A(4,4)*B(4,4)

    explicit_4x4_matinv = 4*B/detA

    !! DEBUGGING
    !
    !ident = real(identity(4))
    !
    !if (any(abs(matmul(A,explicit_4x4_matinv) - ident) .gt. 1E-3)) then
    !    print *, transpose(A)
    !    print *, transpose(explicit_4x4_matinv)
    !    print *, transpose(B)
    !    print *, detA
    !    print *, transpose(matmul(A,explicit_4x4_matinv))
    !    print *, "BREAK"
    !    !stop
    !end if
end function explicit_4x4_matinv

pure function matinv4D(A) ! MUST DOUBLE CHECK
    implicit none
    real, dimension(4,4), intent(in) :: A

    real, dimension(4,4) :: matinv4D

    integer :: i, j
    real, dimension(2,2) :: Ablock
    real, dimension(2,2) :: Bblock
    real, dimension(2,2) :: Cblock
    real, dimension(2,2) :: Dblock
    real, dimension(2,2) :: Ablockinv
    real, dimension(2,2) :: AinvB
    real, dimension(2,2) :: CAinv
    real, dimension(2,2) :: Dminus
    real, dimension(2,2) :: invAblock
    real, dimension(2,2) :: invBblock
    real, dimension(2,2) :: invCblock
    real, dimension(2,2) :: invDblock

    Ablock(1:2,1:2) = A(1:2,1:2)
    Bblock(1:2,1:2) = A(1:2,3:4)
    Cblock(1:2,1:2) = A(3:4,1:2)
    Dblock(1:2,1:2) = A(3:4,3:4)

    Ablockinv(1,1) = Ablock(2,2)
    Ablockinv(1,2) = -Ablock(1,2)
    Ablockinv(2,1) = -Ablock(2,1)
    Ablockinv(2,2) = Ablock(1,1)
    Ablockinv = (1/(Ablock(1,1)*Ablock(2,2)-Ablock(1,2)*Ablock(2,1)))*Ablockinv

    AinvB = matmul(Ablockinv,Bblock)
    CAinv = matmul(Cblock,Ablockinv)
    Dminus = Dblock - matmul(Cblock,matmul(Ablockinv,Bblock))
    invDblock(1,1) = Dminus(2,2)
    invDblock(1,2) = -Dminus(1,2)
    invDblock(2,1) = -Dminus(2,1)
    invDblock(2,2) = Dminus(1,1)
    invDblock = (1/(Dminus(1,1)*Dminus(2,2)-Dminus(1,2)*Dminus(2,1)))*Dminus

    invAblock = Ablockinv + matmul(AinvB,matmul(invDblock,CAinv))
    invBblock = -matmul(AinvB,invDblock)
    invCblock = -matmul(invDblock,CAinv)

    matinv4D(1:2,1:2) = invAblock(1:2,1:2)
    matinv4D(1:2,3:4) = invBblock(1:2,1:2)
    matinv4D(3:4,1:2) = invCblock(1:2,1:2)
    matinv4D(3:4,3:4) = invDblock(1:2,1:2)
end function matinv4D

pure function cross_product(u,v)
    implicit none
    real, dimension(3), intent(in) :: u
    real, dimension(3), intent(in) :: v

    real, dimension(3) :: cross_product

    cross_product(1) = u(2)*v(3)-u(3)*v(2)
    cross_product(2) = u(3)*v(1)-u(1)*v(3)
    cross_product(3) = u(1)*v(2)-u(2)*v(1)
end function cross_product

pure function outer_product(u,v)
    implicit none
    real, dimension(:), intent(in) :: u
    real, dimension(:), intent(in) :: v

    real, dimension(:,:), allocatable :: outer_product

    integer :: N
    integer :: M
    integer :: i, j

    N = size(u,1)
    M = size(v,1)
    allocate(outer_product(N,M))

    do j = 1, M
        do i = 1, N
            outer_product(i,j) = u(i)*v(j)
        end do
    end do
end function outer_product

pure function identity(N)
    implicit none
    integer, intent(in) :: N

    integer, dimension(:,:), allocatable :: identity

    integer :: i

    allocate(identity(N,N))

    identity(:,:) = 0
    forall(i = 1:N)
        identity(i,i) = 1
    end forall
end function identity

function lower_tri_inv(L)
    ! Defined only for N = 4
    implicit none
    real, dimension(:,:), intent(in) :: L

    real, dimension(:,:), allocatable :: lower_tri_inv

    integer :: N

    !if (size(L,1) .ne. size(L,2)) then
    !    print *, &
    !    "math.f08: function lower_tri_inv: L is not square. Program ended."
    !    stop
    !end if
    !
    !if (size(L,1) .ne. 4) then
    !    print *, &
    !    "math.f08: function lower_tri_inv: Dimension of L is not 4. Program ended."
    !end if

    N = size(L,1)

    allocate(lower_tri_inv(N,N))

    lower_tri_inv(:,1) = [1./L(1,1), 0., 0., 0.]
    lower_tri_inv(:,2) = [-L(2,1)/(L(1,1)*L(2,2)), 1./L(2,2), 0., 0.]
    lower_tri_inv(:,3) = [(1./(L(1,1)*L(3,3)))*(-L(3,1)+(1./L(2,2))*L(2,1)*L(3,2)), -L(3,2)/(L(2,2)*L(3,3)), 1./L(3,3), 0.]
    lower_tri_inv(:,4) = [(1./(L(1,1)*L(4,4)))*(-L(4,1)+(1./L(2,2))*L(2,1)*L(4,2)+(1./L(3,3))*L(3,1)*L(4,3)-&
    (1./(L(2,2)*L(3,3)))*L(2,1)*L(3,2)*L(4,3)), (1./(L(2,2)*L(4,4)))*(-L(4,2)+(1./L(3,3))*L(3,2)*L(4,3)),&
    -L(4,3)/(L(3,3)*L(4,4)), 1./L(4,4)]

    lower_tri_inv = transpose(lower_tri_inv)
end function lower_tri_inv

function upper_tri_inv(U)
    ! Defined only for N = 4
    implicit none
    real, dimension(:,:), intent(in) :: U

    real, dimension(:,:), allocatable :: upper_tri_inv

    integer :: N

    !if (size(U,1) .ne. size(U,2)) then
    !    print *, &
    !    "math.f08: function upper_tri_inv: U is not square. Program ended."
    !    stop
    !end if
    !
    !if (size(U,1) .ne. 4) then
    !    print *, &
    !    "math.f08: function upper_tri_inv: Dimension of U is not 4. Program ended."
    !end if

    N = size(U,1)

    allocate(upper_tri_inv(N,N))

    upper_tri_inv = transpose(lower_tri_inv(transpose(U)))
end function upper_tri_inv

real function lower_tri_det(L)
    implicit none
    real, dimension(:,:), intent(in) :: L

    integer :: N
    real, dimension(:), allocatable :: diagonal

    N = size(L,1)

    diagonal = diag_real(L)

    lower_tri_det = product(diagonal)
end function lower_tri_det

real function upper_tri_det(U)
    implicit none
    real, dimension(:,:), intent(in) :: U

    integer :: N
    real, dimension(:), allocatable :: diagonal

    N = size(U,1)

    diagonal = diag_real(U)

    upper_tri_det = product(diagonal(:))
end function upper_tri_det

function LU_decomposition_inv_func(A)
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, dimension(:,:), allocatable :: LU_decomposition_inv_func

    integer :: N

    N = size(A,1)

    allocate(LU_decomposition_inv_func(N,N))

    call LU_decomposition_inv(A, LU_decomposition_inv_func)
end function LU_decomposition_inv_func

subroutine LU_decomposition(A, P, L, U)
    ! Lower triangular matrix has diagonal of unity
    implicit none
    real, dimension(:,:), intent(in) :: A

    integer, dimension(:,:), allocatable, intent(inout) :: P
    real, dimension(:,:), allocatable, intent(inout) :: L
    real, dimension(:,:), allocatable, intent(inout) :: U

    integer :: N
    real, dimension(:), allocatable :: realStor
    integer, dimension(:), allocatable :: intStor
    real, dimension(:), allocatable :: allocStor
    integer :: i, j, k

    N = size(A,1)

    allocate(P(N,N))
    allocate(L(N,N))
    allocate(U(N,N))
    allocate(realStor(N))
    allocate(intStor(N))

    P = identity(N)
    L = real(identity(N))
    U = A

    i = maxloc(abs(U(:,1)), dim = 1)
    if (i .ne. 1) then
        realStor(:) = U(i,:)
        U(i,:) = U(1,:)
        U(1,:) = realStor(:)
        intStor(:) = P(i,:)
        P(i,:) = P(1,:)
        P(1,:) = intStor(:)
    end if
    L(2:N,1) = U(2:N,1)/U(1,1)
    do j = 2, N
        U(j,:) = U(j,:) - (U(j,1)/U(1,1))*U(1,:)
    end do

    do k = 2, N-1
        i = maxloc(abs(U(k:N,k)), dim = 1) + k-1
        if (i .ne. k) then
            realStor(:) = U(i,:)
            U(i,:) = U(k,:)
            U(k,:) = realStor(:)
            intStor(:) = P(i,:)
            P(i,:) = P(k,:)
            P(k,:) = intStor(:)
            allocate(allocStor(k-1))
            allocStor(:) = L(i,1:k-1)
            L(i,1:k-1) = L(k,1:k-1)
            L(k,1:k-1) = allocStor(:)
            deallocate(allocStor)
        end if
        L(k+1:N,k) = U(k+1:N,k)/U(k,k)
        do j = k+1, N
            U(j,k:N) = U(j,k:N) - (U(j,k)/U(k,k))*U(k,k:N)
        end do
    end do
end subroutine LU_decomposition

subroutine LU_decomposition_inv(A, Ainv)
    ! If determinant is wanted as well, use LU_decomposition_inv_and_det
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, dimension(:,:), intent(inout) :: Ainv
    ! NOTE - I had to remove `allocatable' to use in do loop.

    integer :: N
    integer, dimension(:,:), allocatable :: P
    real, dimension(:,:), allocatable :: L
    real, dimension(:,:), allocatable :: U

    N = size(A,1)

    !allocate(Ainv(N,N))

    call LU_decomposition(A, P, L, U)
    L = lower_tri_inv(L)
    U = upper_tri_inv(U)

    Ainv = matmul(matmul(U,L),P)
end subroutine LU_decomposition_inv

subroutine LU_decomposition_det(A, detA)
    ! If inverse is wanted as well, use LU_decomposition_inv_and_det
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, intent(inout) :: detA

    integer :: N
    integer, dimension(:,:), allocatable :: P
    real, dimension(:,:), allocatable :: L
    real, dimension(:,:), allocatable :: U
    integer :: S

    N = size(A,1)

    call LU_decomposition(A, P, L, U)

    S = N - trace_int(P) - 1

    detA = ((-1)**S)*upper_tri_det(U)
end subroutine LU_decomposition_det

subroutine LU_decomposition_inv_and_det(A, Ainv, detA)
    implicit none
    real, dimension(:,:), intent(in) :: A

    real, dimension(:,:), intent(inout) :: Ainv
    ! NOTE - I had to remove `allocatable' to use in do loop.
    real, intent(inout) :: detA

    integer :: N
    integer, dimension(:,:), allocatable :: P
    real, dimension(:,:), allocatable :: L
    real, dimension(:,:), allocatable :: U
    integer :: S

    N = size(A,1)

    S = 1

    !allocate(Ainv(N,N))

    call LU_decomposition(A, P, L, U)

    S = N - trace_int(P) - 1

    detA = ((-1)**S)*upper_tri_det(U)

    L = lower_tri_inv(L)
    U = upper_tri_inv(U)

    Ainv = matmul(matmul(U,L),real(P))
end subroutine LU_decomposition_inv_and_det

subroutine populate_Gauss_Legendre_quadrature(N, xq, wq)
    implicit none
    integer, intent(in) :: N

    real, dimension(:), allocatable, intent(inout) :: xq
    real, dimension(:), allocatable, intent(inout) :: wq

    integer :: i
    character(50) :: fname
    character(50) :: str

    if (N .eq. 128) then
        allocate(xq(N))
        allocate(wq(N))

        write(str,*) N

        fname = "Quadrature sets/Gauss-Legendre_" // trim(adjustl(str)) // "_abscissae.txt"

        open(1, file = trim(fname))
        read(1,*) xq(1:64)
        close(1)

        xq(65:128) = -xq(1:64)

        fname = "Quadrature sets/Gauss-Legendre_" // trim(adjustl(str)) // "_weights.txt"

        open(1, file = trim(fname))
        read(1,*) wq(1:64)
        close(1)

        wq(65:128) = wq(1:64)
    else
        allocate(xq(N))
        allocate(wq(N))

        write(str,*) N

        fname = "Quadrature sets/Gauss-Legendre_" // trim(adjustl(str)) // "_abscissae.txt"

        open(1, file = trim(fname))
        read(1,*) xq
        close(1)

        fname = "Quadrature sets/Gauss-Legendre_" // trim(adjustl(str)) // "_weights.txt"

        open(1, file = trim(fname))
        read(1,*) wq
        close(1)
    end if
end subroutine populate_Gauss_Legendre_quadrature

subroutine populate_tetrahedral_quadrature(N, cq, wq)
    implicit none
    integer, intent(in) :: N

    real, dimension(:,:), allocatable, intent(inout) :: cq
    real, dimension(:), allocatable, intent(inout) :: wq

    integer :: i
    character(50) :: fname
    character(50) :: str

    allocate(cq(3,N))
    allocate(wq(N))

    write(str,*) N

    fname = "Quadrature sets/tetrahedral_" // trim(adjustl(str)) // "_abscissae.txt"

    open(1, file = trim(fname))
    read(1,*) cq
    close(1)

    fname = "Quadrature sets/tetrahedral_" // trim(adjustl(str)) // "_weights.txt"

    open(1, file = trim(fname))
    read(1,*) wq
    close(1)
end subroutine populate_tetrahedral_quadrature

subroutine Gauss_Legendre_quadrature_1D(func, a, b, xq, wq, integral)
    implicit none
    real, intent(in) :: a
    real, intent(in) :: b
    real, dimension(:), intent(in) :: xq
    real, dimension(:), intent(in) :: wq

    real, intent(inout) :: integral

    integer :: i
    !real, dimension(:), allocatable :: BLOCKS
    interface
        real function func(x)
            implicit none
            real, intent(in) :: x
        end function func
    end interface

    !allocate(BLOCKS(size(xq)))
    !!$OMP parallel do shared(BLOCKS) private(i)
    !do i = 1, size(xq)
    !    BLOCKS(i) = 0.5*(b-a)*xq(i)*func(0.5*(b-a)*xq(i)+0.5*(b+a))
    !end do
    !!$OMP end parallel do
    !integral = sum(BLOCKS)

    integral = 0
    do i = 1, size(xq)
        integral = integral + 0.5*(b-a)*wq(i)*func(0.5*(b-a)*xq(i)+0.5*(b+a))
    end do
end subroutine Gauss_Legendre_quadrature_1D

subroutine Gauss_Legendre_quadrature_2D(func, a, b, c, d, xq, wq, integral)
    implicit none
    real, intent(in) :: a
    real, intent(in) :: b
    real, intent(in) :: c
    real, intent(in) :: d
    real, dimension(:), intent(in) :: xq
    real, dimension(:), intent(in) :: wq

    real, intent(inout) :: integral

    integer :: i, j
    !real, dimension(:,:), allocatable :: BLOCKS
    interface
        real function func(x,y)
            implicit none
            real, intent(in) :: x
            real, intent(in) :: y
        end function func
    end interface

    !allocate(BLOCKS(size(xq),size(xq)))
    !!$OMP parallel do shared(BLOCKS) private(i, j)
    !do j = 1, size(xq)
    !    do i = 1, size(xq)
    !        BLOCKS(i,j) = 0.25*(b-a)*(d-c)*&
    !            wq(i)*wq(j)*func(0.5*(b-a)*xq(i)+0.5*(b+a),0.5*(d-c)*xq(j)+0.5*(d+c))
    !    end do
    !end do
    !!$OMP end parallel do
    !integral = sum(BLOCKS)

    integral = 0
    do i = 1, size(xq)
        do j = 1, size(xq)
            integral = integral + 0.25*(b-a)*(d-c)*&
            wq(i)*wq(j)*func(0.5*(b-a)*xq(i)+0.5*(b+a),0.5*(d-c)*xq(j)+0.5*(d+c))
        end do
    end do
end subroutine Gauss_Legendre_quadrature_2D

recursive subroutine Gauss_Legendre_adaptive_quadrature_1D(func, a, b, xq1, wq1, xq2, wq2, tol, integral)
    implicit none
    real, intent(in) :: a
    real, intent(in) :: b
    real, dimension(5), intent(in) :: xq1
    real, dimension(5), intent(in) :: wq1
    real, dimension(8), intent(in) :: xq2
    real, dimension(8), intent(in) :: wq2
    real, intent(in) :: tol

    real, intent(inout) :: integral

    real :: integral1
    real :: integral2
    real, dimension(2) :: integrallist
    interface
        real function func(x)
            implicit none
            real, intent(in) :: x
        end function func
    end interface

    call Gauss_Legendre_quadrature_1D(func, a, b, xq1, wq1, integral1)
    call Gauss_Legendre_quadrature_1D(func, a, b, xq2, wq2, integral2)

    if (abs(integral2-integral1) .gt. tol) then
        call Gauss_Legendre_adaptive_quadrature_1D(func, a, 0.5*(b+a), xq1, wq1, xq2, wq2, tol, integrallist(1))

        call Gauss_Legendre_adaptive_quadrature_1D(func, 0.5*(b+a), b, xq1, wq1, xq2, wq2, tol, integrallist(2))

        integral = integrallist(1) + integrallist(2)
    else
        integral = integral2
    end if
end subroutine Gauss_Legendre_adaptive_quadrature_1D

recursive subroutine Gauss_Legendre_adaptive_quadrature_2D(func, a, b, c, d, xq1, wq1, xq2, wq2, tol, integral)
    implicit none
    real, intent(in) :: a
    real, intent(in) :: b
    real, intent(in) :: c
    real, intent(in) :: d
    real, dimension(5), intent(in) :: xq1
    real, dimension(5), intent(in) :: wq1
    real, dimension(8), intent(in) :: xq2
    real, dimension(8), intent(in) :: wq2
    real, intent(in) :: tol

    real, intent(inout) :: integral

    real :: integral1
    real :: integral2
    real, dimension(4) :: integrallist
    interface
        real function func(x,y)
            implicit none
            real, intent(in) :: x
            real, intent(in) :: y
        end function func
    end interface

    call Gauss_Legendre_quadrature_2D(func, a, b, c, d, xq1, wq1, integral1)
    call Gauss_Legendre_quadrature_2D(func, a, b, c, d, xq2, wq2, integral2)

    ! auto tol?
    ! autotol = 10.**(floor(log10(integral2))-6)
    ! Would ideally like to set only once
    ! So have a logical input to subroutine which is true for the first call, false for recursive calls below
    ! When the logical is true, a different tol is set like above.
    ! When the logical is false, tol is fed.

    if (abs(integral2-integral1) .gt. tol) then
        ! count = count + 1
        call Gauss_Legendre_adaptive_quadrature_2D &
        (func, 0.5*(b+a), b, 0.5*(d+c), d, xq1, wq1, xq2, wq2, tol, integrallist(1))

        call Gauss_Legendre_adaptive_quadrature_2D &
        (func, a, 0.5*(b+a), 0.5*(d+c), d, xq1, wq1, xq2, wq2, tol, integrallist(2))

        call Gauss_Legendre_adaptive_quadrature_2D &
        (func, a, 0.5*(b+a), c, 0.5*(d+c), xq1, wq1, xq2, wq2, tol, integrallist(3))

        call Gauss_Legendre_adaptive_quadrature_2D &
        (func, 0.5*(b+a), b, c, 0.5*(d+c), xq1, wq1, xq2, wq2, tol, integrallist(4))

        integral = sum(integrallist)
    else
        integral = integral2
    end if
end subroutine Gauss_Legendre_adaptive_quadrature_2D

real function interp1d(x,y,x0)
    implicit none
    real, dimension(:), intent(in) :: x
    real, dimension(:), intent(in) :: y
    real, intent(in) :: x0

    integer :: start

    start = findloc(x0 .ge. x(1:size(x)-1) .and. x0 .lt. x(2:size(x)), .true., dim = 1)
    if (start .eq. 0) then
        if (x0 .lt. x(1)) then
            start = 1
        else
            start = size(x)-1
        end if
    end if

    interp1d = (y(start+1)-y(start))*(x0-x(start))/(x(start+1)-x(start)) + y(start)

    if (isnan(interp1d)) then
        print *, "lin interp: NAN found"
        print *, x
        print *, y
        print *, x0
        stop
    end if
end function interp1d

real function log_interp1d(x,y,x0)
    implicit none
    real, dimension(:), intent(in) :: x
    real, dimension(:), intent(in) :: y
    real, intent(in) :: x0

    integer :: start
    real :: Mval
    real :: Bval

    if (x0 .le. 0) then
        print *, "x0 is negative for log_interp1d"
        stop
    end if

    start = findloc(x0 .ge. x(1:size(x)-1) .and. x0 .lt. x(2:size(x)), .true., dim = 1)

    if (start .eq. 0) then
        if (x0 .lt. x(1)) then
            start = 1
        else
            start = size(x)-1
        end if
    end if

    if (y(start) .le. 0) then ! Do regular linear interpolation
        log_interp1d = (y(start+1)-y(start))*(x0-x(start))/(x(start+1)-x(start))+y(start)
    else
        Mval = log10(y(start+1)/y(start))/log10(x(start+1)/x(start))
        Bval = log10(y(start)) - Mval*log10(x(start))
        log_interp1d = (10**(Bval))*x0**(Mval)
    end if

    if (isnan(log_interp1d)) then
        print *, "log interp: NAN found"
        print *, x
        print *, y
        print *, x0
        stop
    end if
end function log_interp1d

subroutine trapezoidal_integration_1D(t, y, x, a, b, integral)
    ! t = 0, plain log-log interpolation
    ! t = 1, log-log interpolation with factor of (A-x)(B-x)
    ! t = 2, log-log interpolation with factor of (A-x)
    ! Requires manually defining public variables:
    ! MATHA: A in int((A-x)*(B-x)*f(x))
    ! MATHB: B in int((A-x)*(B-x)*f(x))
    implicit none
    integer, intent(in) :: t
    real, dimension(:), intent(in) :: y
    real, dimension(:), intent(in) :: x
    real, intent(in) :: a
    real, intent(in) :: b

    real, intent(inout) :: integral

    integer :: bl, i
    real, dimension(:), allocatable :: l_x
    real, dimension(:), allocatable :: l_y
    real, dimension(:), allocatable :: B_x
    real, dimension(:), allocatable :: B_y
    real :: l_ya
    real :: l_yb
    integer :: start
    integer :: finish
    real :: mval
    real :: bval
    integer :: Nb
    logical, dimension(:), allocatable :: logarray
    integer, dimension(:), allocatable :: binds

    ! quicksort x and y? I think I should

    ! Find or insert a into x
    if (any(a .eq. x)) then
        start = findloc(x, a, dim = 1)
        l_ya = y(start)
    else
        start = findloc(a .ge. x(1:size(x)-1) .and. a .lt. x(2:size(x)), .true., dim = 1)
        if (start .eq. 0) then
            if (a .lt. x(1)) then
                start = 1
            else
                print *, "lin-lin, lower limit is greater than max value of x"
                stop
            end if
        end if
        l_ya = interp1d(x,y,a)
    end if

    ! Find or insert b into x
    if (any(b .eq. x)) then
        finish = findloc(x, b, dim = 1)
        l_yb = y(finish)
        finish = finish - 1 ! Shifted so that I can insert l_yb and it won't be a duplicate
    else
        finish = findloc(b .ge. x(1:size(x)-1) .and. b .lt. x(2:size(x)), .true., dim = 1)
        if (finish .eq. 0) then
            if (b .gt. x(size(x))) then
                finish = size(x)
            else
                print *, "lin-lin, upper limit is lower than min value of x"
                stop
            end if
        end if
        l_yb = interp1d(x,y,b)
    end if

    allocate(B_x, source = [a,x(start+1:finish),b])
    allocate(B_y, source = [l_ya,y(start+1:finish),l_yb])

    call move_alloc(B_x, l_x)
    call move_alloc(B_y, l_y)

    ! Integrate
    if (t .eq. 0) then
        integral = 0
        do i = 1, size(l_x)-1
            integral = integral + 0.5*(l_x(i+1)-l_x(i))*(l_y(i+1)+l_y(i))
        end do
    else if (t .eq. 1) then
        integral = 0
        do i = 1, size(l_x)-1
            if (l_x(i) .eq. l_x(i+1)) cycle
            mval = (l_y(i+1)-l_y(i))/(l_x(i+1)-l_x(i))
            bval = l_y(1) - mval*l_x(1)
            integral = integral + 0.25*mval*(l_x(i+1)**4-l_x(i)**4) + &
                (bval-MATHA*mval-MATHB*mval)*(l_x(i+1)**3-l_x(i)**3)/3 + &
                0.5*(MATHA*MATHB*mval-MATHA*bval-MATHB*bval)*(l_x(i+1)**2-l_x(i)**2) + &
                MATHA*MATHB*bval*(l_x(i+1)-l_x(i))
        end do
    else if (t .eq. 2) then
        integral = 0
        do i = 1, size(l_x)-1
            if (l_x(i) .eq. l_x(i+1)) cycle
            mval = (l_y(i+1)-l_y(i))/(l_x(i+1)-l_x(i))
            bval = l_y(1) - mval*l_x(1)
            integral = integral - mval*(l_x(i+1)**3-l_x(i)**3)/3 + &
                0.5*(MATHA*mval-bval)*(l_x(i+1)**2-l_x(i)**2) + &
                MATHA*bval*(l_x(i+1)-l_x(i))
        end do
    end if

    if (isnan(integral)) then
        print *, "lin-lin NAN FOUND"
        print *, l_x
        print *, l_y
    end if
end subroutine trapezoidal_integration_1D

subroutine log_log_integration_1D(t, y, x, a, b, integral)
    ! t = 0, plain log-log interpolation
    ! t = 1, log-log interpolation with factor of (A-x)(B-x)
    ! t = 2, log-log interpolation with factor of (A-x)
    ! Requires manually defining public variables:
    ! MATHA: A in int((A-x)*(B-x)*f(x))
    ! MATHB: B in int((A-x)*(B-x)*f(x))
    implicit none
    integer, intent(in) :: t
    real, dimension(:), intent(in) :: y
    real, dimension(:), intent(in) :: x
    real, intent(in) :: a
    real, intent(in) :: b

    real, intent(inout) :: integral

    integer :: i
    real, dimension(:), allocatable :: l_x
    real, dimension(:), allocatable :: l_y
    real, dimension(:), allocatable :: B_x
    real, dimension(:), allocatable :: B_y
    real :: l_ya
    real :: l_yb
    integer :: start
    integer :: finish
    real :: Mval
    real :: Bval

    ! quicksort x and y? Make logical intent in for option?

    ! Find or insert a into x
    if (any(a .eq. x)) then
        start = findloc(x, a, dim = 1)
        l_ya = y(start)
    else
        start = findloc(a .ge. x(1:size(x)-1) .and. a .lt. x(2:size(x)), .true., dim = 1)
        if (start .eq. 0) then
            if (a .lt. x(1)) then
                start = 1
            else
                print *, "log-log, lower limit is greater than max value of x"
                print *, x
                print *, y
                print *, a, b
                stop
            end if
        end if
        l_ya = log_interp1d(x,y,a)
    end if

    ! Find or insert b into x
    if (any(b .eq. x)) then
        finish = findloc(x, b, dim = 1)
        l_yb = y(finish)
        finish = finish - 1 ! Shifted so that I can insert l_yb and it won't be a duplicate
    else
        finish = findloc(b .ge. x(1:size(x)-1) .and. b .lt. x(2:size(x)), .true., dim = 1)
        if (finish .eq. 0) then
            if (b .gt. x(size(x))) then
                finish = size(x)
            else
                print *, "log-log, upper limit is lower than min value of x"
            end if
        end if
        l_yb = log_interp1d(x,y,b)
    end if

    allocate(B_x, source = [a,x(start+1:finish),b])
    allocate(B_y, source = [l_ya,y(start+1:finish),l_yb])

    call move_alloc(B_x, l_x)
    call move_alloc(B_y, l_y)

    ! Check if all entries of l_y are zero. If so, integral is zero.
    if (all(l_y .eq. 0)) then
        integral = 0
        return
    end if

    ! Check if l_y is zero after the first entry.
    if (l_y(2) .eq. 0) then
        start = findloc(l_y .gt. 0, .true., dim = 1) - 1

        allocate(B_x, source = l_x(start:size(l_x)))
        allocate(B_y, source = l_y(start:size(l_y)))

        call move_alloc(B_x, l_x)
        call move_alloc(B_y, l_y)
    end if

    !! Check for zeros/negatives in l_x beyond l_x(1)
    !if (any(l_x(2:size(l_x)) .le. 0)) then
    !    print *, "log-log, zero or negative in l_x after l_x(1). How??"
    !    print *, "this is a and b:"
    !    print *, a, b
    !    print *, "l_x:"
    !    print *, l_x
    !    stop
    !end if

    ! If there are negatives in l_y, just use trapezoidal integration.
    ! MAKE TRAPEZOIDAL INTEGRATION RUN ONLY OVER THE NEGATIVE INDICES, THEN FINISH UP WITH LOG
    if (any(l_y .lt. 0)) then
        call trapezoidal_integration_1D(t, l_y, l_x, a, b, integral)
        return
    end if

    ! Initialize integral scalar with trapezoidal term if y starts with zeros
    nanprevent: if (l_y(1) .eq. 0) then
        Mval = (l_y(2)-l_y(1))/(l_x(2)-l_x(1))
        Bval = l_y(1) - Mval*l_x(1)
        if (t .eq. 0) then
            integral = 0.5*(l_x(2)-l_x(1))*(l_y(2)+l_y(1))
        else if (t .eq. 1) then
            integral = 0.25*Mval*(l_x(2)**4-l_x(1)**4) + &
                (Bval - MATHA*Mval - MATHB*Mval)*(l_x(2)**3-l_x(1)**3)/3 + &
                0.5*(MATHA*MATHB*Mval - MATHA*Bval - MATHB*Bval)*(l_x(2)**2-l_x(1)**2) + &
                (MATHA*MATHB*Bval)*(l_x(2)-l_x(1))
        else if (t .eq. 2) then
            integral = -Mval*(l_x(2)**3-l_x(1)**3)/3 + &
                0.5*(MATHA*Mval-Bval)*(l_x(2)**2-l_x(1)**2) + &
                MATHA*Bval*(l_x(2)-l_x(1))
        end if

        allocate(B_x, source = l_x(2:size(l_x)))
        allocate(B_y, source = l_y(2:size(l_y)))

        call move_alloc(B_x, l_x)
        call move_alloc(B_y, l_y)
    else
        integral = 0
    end if nanprevent

    if (t .eq. 0) then
        do i = 1, size(l_x)-1
            if (l_x(i) .eq. l_x(i+1)) cycle
            Mval = log10(l_y(i+1)/l_y(i))/log10(l_x(i+1)/l_x(i))
            Bval = log10(l_y(i)/(l_x(i)**Mval))
            integral = integral + (10**Bval)*(l_x(i+1)**(Mval+1)-l_x(i)**(Mval+1))/(Mval+1)
        end do
    else if (t .eq. 1) then
        do i = 1, size(l_x)-1
            if (l_x(i) .eq. l_x(i+1)) cycle
            Mval = log10(l_y(i+1)/l_y(i))/log10(l_x(i+1)/l_x(i))
            Bval = log10(l_y(i)/(l_x(i)**Mval))
            integral = integral + (10**Bval)*&
                (MATHA*MATHB*(l_x(i+1)**(Mval+1)-l_x(i)**(Mval+1))/(Mval+1) - &
                (MATHA+MATHB)*(l_x(i+1)**(Mval+2)-l_x(i)**(Mval+2))/(Mval+2) + &
                (l_x(i+1)**(Mval+3)-l_x(i)**(Mval+3))/(Mval+3))
        end do
    else if (t .eq. 2) then
        do i = 1, size(l_x)-1
            if (l_x(i) .eq. l_x(i+1)) cycle
            Mval = log10(l_y(i+1)/l_y(i))/log10(l_x(i+1)/l_x(i))
            Bval = log10(l_y(i)/(l_x(i)**Mval))
            integral = integral + (10**Bval)*&
                (-(l_x(i+1)**(Mval+2)-l_x(i)**(Mval+2))/(Mval+2) + &
                MATHA*(l_x(i+1)**(Mval+1)-l_x(i)**(Mval+1))/(Mval+1))
        end do
    end if

    if (isnan(integral)) then
        !print *, "log_log_integration_1D found a NaN. FIX THIS BEFORE PUBLISHING!!"
        !print *, t
        !print *, l_x
        !print *, l_y
        !print *, x
        !print *, y
        !print *, a, b
        !print *, l_ya, l_yb

        start = findloc(l_y .gt. 0, .true., dim = 1) - 1

        allocate(B_x, source = l_x(start:size(l_x)))
        allocate(B_y, source = l_y(start:size(l_y)))

        call move_alloc(B_x, l_x)
        call move_alloc(B_y, l_y)

        !print *, start
        stop
    end if
end subroutine log_log_integration_1D

recursive subroutine qsort(v)
    implicit none
    real, dimension(:), intent(inout) :: v

    integer :: i, j
    real :: x1
    real :: x2
    integer :: N

    N = size(v,1)
    i = 1
    j = N
    x1 = v((N+1)/2)

    do
        do while (v(i) .lt. x1)
            i = i + 1
        end do
        do while (v(j) .gt. x1)
            j = j - 1
        end do

        if (i .ge. j) exit

        x2 = v(i)
        v(i) = v(j)
        v(j) = x2
        i = i + 1
        j = j - 1
    end do

    if (i .gt. 2) call qsort(v(1:i-1))

    if (j .lt. N-1) call qsort(v(j+1:N))
end subroutine qsort

recursive subroutine qsort_int(v)
    implicit none
    integer, dimension(:), intent(inout) :: v

    integer :: i, j
    integer :: x1
    integer :: x2
    integer :: N

    N = size(v,1)
    i = 1
    j = N
    x1 = v((N+1)/2)

    do
        do while (v(i) .lt. x1)
            i = i + 1
        end do
        do while (v(j) .gt. x1)
            j = j - 1
        end do

        if (i .ge. j) exit

        x2 = v(i)
        v(i) = v(j)
        v(j) = x2
        i = i + 1
        j = j - 1
    end do

    if (i .gt. 2) call qsort_int(v(1:i-1))

    if (j .lt. N-1) call qsort_int(v(j+1:N))
end subroutine qsort_int

subroutine make_SBR(A, zerotol, SBR_A, Ia, Ja)
    implicit none
    real, dimension(:,:), intent(in) :: A
    real, intent(in) :: zerotol

    real, dimension(:), allocatable, intent(inout) :: SBR_A
    integer, dimension(:), allocatable, intent(inout) :: Ia
    integer, dimension(:), allocatable, intent(inout) :: Ja

    integer :: k, i
    integer :: N
    real, dimension(:,:), allocatable :: l_A
    logical, dimension(:,:), allocatable :: columnmask
    integer :: NonZ

    N = size(A,1)
    allocate(Ia(N+1))
    allocate(l_A(N,N))
    allocate(columnmask(N,N))

    l_A = transpose(A) ! SLOW?

    do k = 1, N
        columnmask(1:N,k) = abs(l_A(1:N,k)) .le. zerotol
    end do

    ! Ia(row) counts the number of nonzero els up to BUT NOT including the row "row," + 1.
    Ia(1) = 1
    do k = 2, N+1
        Ia(k) = Ia(k-1) + count(columnmask(1:N,k-1))
    end do

    NonZ = Ia(N+1) - 1

    allocate(Ja(NonZ))

    do k = 1, N
        Ja(Ia(k):Ia(k+1)-1) = pack([(i,i=1,N)], mask = columnmask(1:N,k))
    end do

    allocate(SBR_A(NonZ))

    do k = 1, N
        SBR_A(Ia(k):Ia(k+1)-1) = pack(A(1:N,k), mask = columnmask(1:N,k))
    end do
end subroutine make_SBR

subroutine SBR_matmul(SBR_A, Ia, Ja, v, w)
    implicit none
    real, dimension(:), intent(in) :: SBR_A
    integer, dimension(:), intent(in) :: Ia
    integer, dimension(:), intent(in) :: Ja
    real, dimension(:), intent(in) :: v

    real, dimension(:), allocatable, intent(inout) :: w

    integer :: k, ki
    integer :: N
    integer :: NonZ

    N = size(v)
    NonZ = size(Ja)

    allocate(w(N))

    do k = 1, N
        w(k) = dot_product(SBR_A(Ia(k):Ia(k+1)-1),v(Ja(Ia(k):Ia(k+1)-1)))
    end do
end subroutine SBR_matmul

end module math