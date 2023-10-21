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

module shape_functions
    use math
    use physics
    use user_input
implicit none

contains
subroutine tetrahedral_geometry_and_shape_functions(r, area, vol, a, bvector)
    implicit none
    real, dimension(:,:,:), allocatable, intent(in) :: r

    real, dimension(:,:), allocatable, intent(inout) :: area
    real, dimension(:), allocatable, intent(inout) :: vol
    real, dimension(:,:), allocatable, intent(inout) :: a
    real, dimension(:,:,:), allocatable, intent(inout) :: bvector

    integer :: dir, e, f, k
    integer, dimension(:,:), allocatable :: cyclematrix
    real :: alpha
    real :: beta
    real :: gamma
    real :: side
    real, dimension(:,:), allocatable :: volmat
    real, dimension(:,:), allocatable :: b
    real, dimension(:,:), allocatable :: c
    real, dimension(:,:), allocatable :: d
    real, dimension(:,:), allocatable :: localcofactor
    real :: sumtest

    allocate(cyclematrix(NFe,3))
    allocate(area(NFe,NE))
    allocate(vol(NE))
    allocate(a(NKe,NE))
    allocate(bvector(3,NKe,NE))
    allocate(volmat(NKe,NKe))
    allocate(b(NKe,NE))
    allocate(c(NKe,NE))
    allocate(d(NKe,NE))
    allocate(localcofactor(NKe,NKe))

    cyclematrix(1,:) = [2,3,4]
    cyclematrix(2,:) = [1,3,4]
    cyclematrix(3,:) = [1,2,4]
    cyclematrix(4,:) = [1,2,3]

    do e = 1, NE
        do f = 1, NFe
            !alpha = norm2(r(:,cyclematrix(f,1),e) - r(:,cyclematrix(f,2),e))
            !beta = norm2(r(:,cyclematrix(f,1),e) - r(:,cyclematrix(f,3),e))
            !gamma = norm2(r(:,cyclematrix(f,2),e) - r(:,cyclematrix(f,3),e))
            !side = (alpha+beta+gamma)/2.0
            !area(f,e) = sqrt(side*(side-alpha)*(side-beta)*(side-gamma))
            area(f,e) = 0.5*norm2(cross_product(r(:,cyclematrix(f,3),e)-r(:,cyclematrix(f,1),e),&
                                                r(:,cyclematrix(f,2),e)-r(:,cyclematrix(f,1),e)))
        end do
        volmat(2:NKe,:) = 0.0
        volmat(1,:) = 1.0
        do k = 1, NKe
            do dir = 1, 3
                volmat(1+dir,k) = r(dir,k,e)
            end do
        end do
        localcofactor = cofactor(volmat)
        do k = 1, NKe
            a(k,e) = localcofactor(1,k)
            b(k,e) = localcofactor(2,k)
            c(k,e) = localcofactor(3,k)
            d(k,e) = localcofactor(4,k)
        end do
        vol(e) = matdet(volmat)/6.0
    end do

    bvector(1,:,:) = b
    bvector(2,:,:) = c
    bvector(3,:,:) = d
end subroutine tetrahedral_geometry_and_shape_functions

end module shape_functions