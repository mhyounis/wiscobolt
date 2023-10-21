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

module spatial_inner_products
    use math
    use physics
    use user_input
implicit none

contains
subroutine inner_products_with_I1inv &
(eprime, fprime, kprime, area, vol, bvector, I1invI2, I1invI2f, I1invI3vector)
    implicit none
    integer, dimension(:,:), intent(in) :: eprime
    integer, dimension(:,:), intent(in) :: fprime
    integer, dimension(:,:,:), intent(in) :: kprime
    real, dimension(:,:), intent(in) :: area
    real, dimension(:), intent(in) :: vol
    real, dimension(:,:,:), intent(in) :: bvector

    real, dimension(:,:,:,:), allocatable, intent(inout) :: I1invI2
    real, dimension(:,:,:,:), allocatable, intent(inout) :: I1invI2f
    real, dimension(:,:,:,:), allocatable, intent(inout) :: I1invI3vector ! (k,kp,dir,e)

    integer :: dir, e, f, k, kp
    real, dimension(:,:), allocatable :: I1inv
    real, dimension(:,:,:,:), allocatable :: I2
    real, dimension(:,:,:,:), allocatable :: I2f
    real, dimension(:,:,:,:), allocatable :: I3vector
    real, dimension(:,:,:,:), allocatable :: testtemp

    allocate(I1invI2(NKe,NKe,NFe,NE))
    allocate(I1invI2f(NKe,NKe,NFe,NE))
    allocate(I1invI3vector(NKe,NKe,3,NE))
    allocate(I1inv(NKe,NKe))
    allocate(I2(NKe,NKe,NFe,NE))
    allocate(I2f(NKe,NKe,NFe,NE))
    allocate(I3vector(NKe,NKe,3,NE))

    I1inv = 20*identity(NKe) - 4

    do e = 1, NE
        do f = 1, NFe
            do kp = 1, NKe
                do k = 1, NKe
                    I2(k,kp,f,e) = area(f,e)*&
                    (1-real(Kron_delta(k,f)))*(1-real(Kron_delta(kp,f)))*&
                    (1+real(Kron_delta(k,kp)))/12
                    if (fprime(f,e) .eq. 0) then
                        I2f(k,kp,f,e) = I2(k,kp,f,e)
                    else
                        I2f(k,kp,f,e) = area(f,e)*&
                        (1-real(Kron_delta(kp,fprime(f,e))))*(1-real(Kron_delta(k,f)))*&
                        (1+real(Kron_delta(kp,kprime(k,f,e))))/12
                    end if
                end do
            end do
        end do
    end do

    do dir = 1, 3
        do kp = 1, NKe
            do k = 1, NKe
                I3vector(k,kp,dir,1:NE) = bvector(dir,k,1:NE)/24
            end do
        end do
    end do

    do e = 1, NE
        do f = 1, NFe
            I1invI2(:,:,f,e) = matmul(I1inv,I2(:,:,f,e))/vol(e) ! Just write these matmuls out
            I1invI2f(:,:,f,e) = transpose(matmul(I1inv,I2f(:,:,f,e)))/vol(e)
        end do
        do dir = 1, 3
            I1invI3vector(:,:,dir,e) = matmul(I1inv,I3vector(:,:,dir,e))/vol(e)
        end do
    end do
end subroutine inner_products_with_I1inv

subroutine wiscoslab_inner_products(dz, mu, GpFterm, Fsweep)
    implicit none
    real, dimension(:), intent(in) :: dz
    real, dimension(:), intent(in) :: mu

    real, dimension(:,:,:,:), allocatable, intent(inout) :: GpFterm ! (k,kp,e,i)
    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: Fsweep ! (k,kp,ip,i,f)

    integer :: e, f, i, k, kp
    real, dimension(:,:,:), allocatable :: Fterm
    real, dimension(:,:), allocatable :: I1inv
    real, dimension(:,:), allocatable :: I2f
    real, dimension(:,:), allocatable :: I3
    real, dimension(:,:,:,:), allocatable :: tempB
    integer :: estart
    integer :: efin
    integer :: einc

    allocate(GpFterm(NKe,NKe,NE,Nmu))
    allocate(Fsweep(NKe,NKe,NE,Nmu,NFe))
    allocate(Fterm(NKe,NKe,Nmu))
    allocate(I1inv(NKe,NKe))
    allocate(I2f(NKe,NKe))
    allocate(I3(NKe,NKe))
    allocate(tempB(NKe,NKe,NE,2))

    I1inv = 2*(-1+3*identity(NKe))

    do i = 1, Nmu
        do kp = 1, NKe
            do k = 1, NKe
                Fterm(k,kp,i) = &
                    merge(0.0, mu(i), mu(i) .lt. 0)*Kron_delta(k,1)*Kron_delta(kp,1) - &
                    merge(mu(i), 0.0, mu(i) .lt. 0)*Kron_delta(k,2)*Kron_delta(kp,2)
            end do
        end do

        Fterm(1:NKe,1:NKe,i) = matmul(I1inv,Fterm(1:NKe,1:NKe,i)) ! Convert to elemental w/ factor of 1/dz(e)
    end do

    do k = 1, NKe
        I3(k,1:NKe) = 0.5*(-1)**(k+1)
    end do

    I3 = matmul(I1inv,I3)

    do i = 1, Nmu
        do e = 1, NE
            GpFterm(1:NKe,1:NKe,e,i) = (-mu(i)*I3 + Fterm(1:NKe,1:NKe,i))/dz(e)
        end do
    end do

    do f = 1, NFe
        do i = 1, Nmu
            do kp = 1, NKe
                do k = 1, NKe
                    I2f(k,kp) = &
                        merge(mu(i), 0.0, mu(i) .lt. 0)*Kron_delta(k,1)*Kron_delta(kp,2)*Kron_delta(f,1) - &
                        merge(0.0, mu(i), mu(i) .lt. 0)*Kron_delta(k,2)*Kron_delta(kp,1)*Kron_delta(f,2)
                end do
            end do

            do e = 1, NE
                Fsweep(1:NKe,1:NKe,e,i,f) = matmul(I1inv,I2f)/dz(e)
            end do
        end do
    end do

    ! Makes Fsweep for upstr. bdy = Fterm
    do i = 1, Nmu
        if (mu(i) .gt. 0) then
            do kp = 1, NKe
                do k = 1, NKe
                    I2f(k,kp) = &
                        -mu(i)*Kron_delta(k,2)*Kron_delta(kp,2)
                end do
            end do
            Fsweep(1:NKe,1:NKe,NE,i,2) = matmul(I1inv,I2f)/dz(NE)
        else
            do kp = 1, NKe
                do k = 1, NKe
                    I2f(k,kp) = &
                        mu(i)*Kron_delta(k,1)*Kron_delta(kp,1)
                end do
            end do
            Fsweep(1:NKe,1:NKe,1,i,1) = matmul(I1inv,I2f)/dz(1)
        end if
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
        tempB = Fsweep(1:NKe,1:NKe,1:NE,i,1:2) ! NOTE: a bit odd to have the faces involved, but doesn't matter since they will be zero when not needed
        Fsweep(1:NKe,1:NKe,1:NE,i,1:2) = tempB(1:NKe,1:NKe,estart:efin:einc,1:2)
    end do
end subroutine wiscoslab_inner_products

end module spatial_inner_products