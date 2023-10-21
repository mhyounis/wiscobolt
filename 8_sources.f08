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

module sources
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

! Whole module involves A LOT of un-optimized operations. Too many do loops.
! Not fully parallelized.
! Especially bad for Nmats > 1
! Make ele source module?? Would actually clean up a lot of this.
! Take probably one day to implement. Just run on super computer. Huge for peace of mind...

contains
subroutine single_material_ray_trace_vector(R, kvec, rglobal, Cfkf, bdyfc, tol, ells)
    implicit none
    real, dimension(:), intent(in) :: R
    real, dimension(:), intent(in) :: kvec
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cfkf
    integer, dimension(:,:), intent(in) :: bdyfc
    real, intent(in) :: tol

    real :: ells

    integer :: f
    real :: ell
    real, dimension(3) :: vert0
    real, dimension(3) :: vert1
    real, dimension(3) :: vert2
    real, dimension(3) :: edge1
    real, dimension(3) :: edge2
    real, dimension(3) :: pvec
    real, dimension(3) :: tvec
    real, dimension(3) :: qvec
    real, dimension(3) :: Normal
    real :: det
    real :: u
    real :: v
    real :: t
    logical :: success

    ells = 0
    ell = norm2(kvec)
    success = .false.
    do f = 1, size(bdyfc,1)
        vert0 = rglobal(:,Cfkf(bdyfc(f,3),1))
        vert1 = rglobal(:,Cfkf(bdyfc(f,3),3)) ! SWAPPED
        vert2 = rglobal(:,Cfkf(bdyfc(f,3),2)) ! SWAPPED

        edge1 = vert1 - vert0
        edge2 = vert2 - vert0

        ! Check if face is backfacing. If so, cycle
        Normal = cross_product(edge1,edge2)
        det = dot_product(Normal, kvec)
        if (det .gt. -tol) cycle ! Change this back to det .gt. tol? This is weeding out parallel faces

        pvec = cross_product(kvec, edge2)

        det = dot_product(edge1, pvec)

        if (det .lt. tol) cycle

        tvec = R - vert0

        u = dot_product(tvec, pvec)

        if (u .lt. -tol .or. u .gt. det + tol) cycle

        qvec = cross_product(tvec, edge1)

        v = dot_product(kvec, qvec)

        if (v .lt. -tol .or. u + v .gt. det + tol) cycle

        t = dot_product(edge2, qvec)/det

        success = .true.
        exit
    end do

    if (success) then
        ells = (1-t)*ell
    else
        print *, "STOP!"
        print *, "MODULE 8, SUBROUTINE single_material_ray_trace_vector:"
        print *, "A node in your mesh was not successfully traced."
        print *, "Check your mesh and beam settings for issues."
        print *, "If the problem persists, you could try turning off"
        print *, "the FCS method in the input file."
        print *, "PROGRAM ENDING."
        stop
    end if
end subroutine single_material_ray_trace_vector

subroutine single_material_ray_trace_mesh(rglobal, Cfkf, bdyfc, tol, ells)
    !!! BIG NOTE: Instead of using bdyfc, use an array that keeps only faces pointing towards source? But I think this might depend on beam angular dist? Still should be singular criterion.
    !!! Similar idea to what I'm doing for heterogeneous mesh.
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cfkf
    integer, dimension(:,:), intent(in) :: bdyfc
    real, intent(in) :: tol

    real, dimension(:,:), allocatable :: ells

    integer :: f, i, k
    real, dimension(3) :: R
    real, dimension(3) :: kvec

    allocate(ells(size(nodesinbeam),1))

    ells = 0
    R = beam_origin
    do k = 1, size(nodesinbeam)
        if (beam_angular_dist .eq. "planar") then
            kvec = dot_product(rglobal(:,nodesinbeam(k))-beam_origin,k0)*k0
            R = rglobal(:,nodesinbeam(k)) - kvec
        else if (beam_angular_dist .eq. "spherical") then
            kvec = rglobal(:,nodesinbeam(k)) - R
        end if

        call single_material_ray_trace_vector(R, kvec, rglobal, Cfkf, bdyfc, tol, ells(k,1))
    end do
end subroutine single_material_ray_trace_mesh

    !subroutine single_material_ray_tracing(rglobal, Cfkf, bdyfc, tol, ells)
    !    implicit none
    !    real, dimension(:,:), intent(in) :: rglobal
    !    integer, dimension(:,:), intent(in) :: Cfkf
    !    integer, dimension(:,:), intent(in) :: bdyfc
    !    real, intent(in) :: tol
    !
    !    real, dimension(:,:), allocatable :: ells
    !
    !    integer :: f, i, k
    !    real, dimension(3) :: R
    !    real, dimension(3) :: kvec
    !    real, dimension(3) :: vert0
    !    real, dimension(3) :: vert1
    !    real, dimension(3) :: vert2
    !    real, dimension(3) :: edge1
    !    real, dimension(3) :: edge2
    !    real, dimension(3) :: pvec
    !    real, dimension(3) :: tvec
    !    real, dimension(3) :: qvec
    !    real, dimension(3) :: Normal
    !    real :: det
    !    real :: u
    !    real :: v
    !    real :: t
    !    real :: ell
    !    real, dimension(:,:), allocatable :: rcurs
    !    integer :: success
    !    real :: ftest
    !    real, dimension(3) :: temp_rcurs
    !    logical :: diditwork
    !
    !    allocate(ells(size(nodesinbeam),1))
    !
    !    ells = 0
    !    R = beam_origin
    !    success = 0
    !
    !    do k = 1, size(nodesinbeam)
    !        diditwork = .false.
    !        if (beam_angular_dist .eq. "planar") then
    !            ell = dot_product(rglobal(:,nodesinbeam(k))-beam_origin,k0)
    !            kvec = ell*k0
    !            R = beam_origin + &
    !                (rglobal(:,nodesinbeam(k))-beam_origin) - &
    !                kvec
    !        else if (beam_angular_dist .eq. "spherical") then
    !            kvec = rglobal(:,nodesinbeam(k)) - R
    !            ell = norm2(kvec)
    !        end if
    !        do f = 1, size(bdyfc,1)
    !            vert0 = rglobal(:,Cfkf(bdyfc(f,3),1))
    !            vert1 = rglobal(:,Cfkf(bdyfc(f,3),3)) ! SWAPPED
    !            vert2 = rglobal(:,Cfkf(bdyfc(f,3),2)) ! SWAPPED
    !
    !            edge1 = vert1 - vert0
    !            edge2 = vert2 - vert0
    !
    !            ! Check if face is parallel or backfacing
    !            Normal = cross_product(edge1,edge2)
    !            det = dot_product(Normal, kvec)
    !            ! ### Criterion for parallel
    !            !if (det)
    !            ! ### Criterion for backfacing
    !            if (det .gt. tol) cycle
    !
    !            pvec = cross_product(kvec, edge2)
    !
    !            det = dot_product(edge1, pvec)
    !
    !            if (det .lt. tol) cycle
    !
    !            tvec = R - vert0
    !
    !            u = dot_product(tvec, pvec)
    !
    !            if (u .lt. 0 .or. u .gt. det) cycle
    !
    !            qvec = cross_product(tvec, edge1)
    !
    !            v = dot_product(kvec, qvec)
    !
    !            if (v .lt. 0 .or. u + v .gt. det) cycle
    !
    !            t = dot_product(edge2, qvec)/det
    !
    !            diditwork = .true.
    !            success = success + 1
    !            exit
    !        end do
    !
    !        if (success .eq. k) then
    !            ells(k,1) = (1-t)*ell
    !            !temp_rcurs = rglobal(:,nodesinbeam(k)) - beam_origin
    !            !print *, k, ells(k,1), (1 + 100.0/temp_rcurs(3))*norm2(temp_rcurs)
    !        end if
    !
    !        if (.not. diditwork) then
    !            print *, nodesinbeam(k), rglobal(:,nodesinbeam(k))
    !        end if
    !    end do
    !
    !    if (success .ne. size(nodesinbeam)) then
    !        print *, "STOP!"
    !        print *, "MODULE 8, SUBROUTINE single_material_ray_tracing:"
    !        print *, "Some nodes were not successfully traced. The total number of nodes to be traced was:"
    !        print *, size(nodesinbeam)
    !        print *, "while the total number of nodes that were not traced was:"
    !        print *, size(nodesinbeam) - success
    !        print *, "Check your mesh and beam for issues."
    !        print *, "PROGRAM ENDING."
    !        stop
    !    end if
    !end subroutine single_material_ray_tracing

subroutine heterogeneous_medium_ray_trace_vector(R, kvec, rglobal, Cfkf, matfc, tol, ells)
    implicit none
    real, dimension(:), intent(in) :: R
    real, dimension(:), intent(in) :: kvec
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cfkf
    integer, dimension(:,:), intent(in) :: matfc
    real, intent(in) :: tol

    real, dimension(:), allocatable :: ells ! Must be pre-allocated with size Nmats

    integer :: f, i, mat
    real :: ell
    real, dimension(3) :: vert0
    real, dimension(3) :: vert1
    real, dimension(3) :: vert2
    real, dimension(3) :: edge1
    real, dimension(3) :: edge2
    real, dimension(3) :: pvec
    real, dimension(3) :: tvec
    real, dimension(3) :: qvec
    real, dimension(3) :: Normal
    real :: det
    real :: u
    real :: v
    real :: tval
    real, dimension(:), allocatable :: t
    real, dimension(:), allocatable :: B_t
    integer :: N
    real :: t0
    real :: tN
    logical :: success
    integer :: sizeofmat

    t0 = 0
    tN = 1

    ells = 0
    ell = norm2(kvec)
    success = .true. ! Not a mistake. Later on taking success = success .and. (something).
    do mat = 1, Nmats
        allocate(t(0))

        sizeofmat = findloc(matfc(:,mat), 0, dim = 1) - 1

        do f = 1, sizeofmat
            vert0 = rglobal(:,Cfkf(matfc(f,mat),1))
            vert1 = rglobal(:,Cfkf(matfc(f,mat),3)) ! SWAPPED
            vert2 = rglobal(:,Cfkf(matfc(f,mat),2)) ! SWAPPED

            edge1 = vert1 - vert0
            edge2 = vert2 - vert0

            Normal = cross_product(edge1,edge2)
            det = dot_product(Normal, kvec)
            ! NOTES: No culling since material interfaces have normals with direction up to +/- 1,
            ! also removing parallel faces.
            if (abs(det) .lt. tol) cycle

            pvec = cross_product(kvec, edge2)

            det = dot_product(edge1, pvec)

            if (det .lt. tol) cycle ! If stuff doesn't work, investigate this step

            tvec = R - vert0

            u = dot_product(tvec, pvec)

            if (u .lt. -tol .or. u .gt. det + tol) cycle

            qvec = cross_product(tvec, edge1)

            v = dot_product(kvec, qvec)

            if (v .lt. -tol .or. u + v .gt. det + tol) cycle

            tval = dot_product(edge2, qvec)/det

            allocate(B_t(size(t)+1))
            B_t(1:size(t)) = t
            B_t(size(t)+1) = tval
            call move_alloc(B_t, t)
        end do

        success = success .and. size(t) .gt. 0

        call qsort(t)

        ! DOUBLE CHECK. I don't know if this is how it should be.
        ! I used stupid conventions again. When not tired, redo conventions.
        N = size(t) + 1

        allocate(B_t(N+2))
        B_t(1) = t0
        B_t(2:N+1) = t
        B_t(N+2) = tN
        call move_alloc(B_t, t)

        ells(mat) = ell
        do i = 1, floor(0.5*(N+1))
            ells(mat) = ells(mat) - (t(2*i-1 + 1)-t(2*i-2 + 1))*ell ! +1 to shift indexing with t(0) stuff
        end do

        deallocate(t)
    end do

    if (.not. success) then
        print *, "STOP!"
        print *, "MODULE 8, SUBROUTINE heterogeneous_medium_ray_trace_vector:"
        print *, "A node was not successfully traced."
        print *, "Check your mesh and beam settings for issues."
        print *, "If the problem persists, you could try turning off"
        print *, "the FCS method in the input file."
        print *, "PROGRAM ENDING."
        stop
    end if
end subroutine heterogeneous_medium_ray_trace_vector

subroutine heterogeneous_medium_ray_trace_mesh(rglobal, Cfkf, matfc, tol, ells)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cfkf
    integer, dimension(:,:), intent(in) :: matfc
    real, intent(in) :: tol

    real, dimension(:,:), allocatable :: ells

    integer :: f, i, k
    real, dimension(3) :: R
    real, dimension(3) :: kvec
    real, dimension(:), allocatable :: l_ells

    allocate(ells(size(nodesinbeam),Nmats))
    allocate(l_ells(Nmats))

    ells = 0
    R = beam_origin
    do k = 1, size(nodesinbeam)
        if (beam_angular_dist .eq. "planar") then
            kvec = dot_product(rglobal(:,nodesinbeam(k))-beam_origin,k0)*k0
            R = rglobal(:,nodesinbeam(k)) - kvec
        else if (beam_angular_dist .eq. "spherical") then
            kvec = rglobal(:,nodesinbeam(k)) - R
        end if

        call heterogeneous_medium_ray_trace_vector &
            (R, kvec, rglobal, matfc, Cfkf, tol, l_ells)
        ells(k,:) = l_ells
    end do
end subroutine heterogeneous_medium_ray_trace_mesh

subroutine TEMP_MANUAL_ray_tracing(rglobal, ells)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal

    real, dimension(:,:), allocatable, intent(inout) :: ells

    integer :: k
    real, dimension(3) :: rcurs
    real, dimension(3) :: rproj
    integer :: NB

    !! Al/Au/Al
    !real :: tmat1 = 0.015917778
    !real :: tmat2 = 0.010611852
    !real :: Lz = 0.18
    !real :: S

    ! Just Al, but with fully open beam
    !real :: Lz = 0.20
    real :: Lz = 1.0
    real :: S = 100.0

    !real :: tmat1 = 1.0
    !real :: tmat2 = 2.0
    !real :: tmat3 = 7.0
    !real :: Lz = 30.0
    !real :: S = 100.0

    real :: rho = 5.0E-6 ! cm
    real :: R = 5.0 ! cm
    !real :: S = 105.0 ! cm
    real :: mu, d, theta1, theta2, m, p1, p2, water_ell, water_entry, H, ell
    real, dimension(2) :: PP1, PP2

    print *, "YOU ARE USING A MANUALLY RAY TRACED MESH. ARE YOU SURE IT'S CORRECT FOR YOUR MESH???"

    NB = Nmats

    allocate(ells(size(nodesinbeam),NB))
    ells = 0

    !! Planar beam, layers
    !S = beam_origin(3) - Lz
    !do k = 1, size(nodesinbeam)
    !    rcurs = rglobal(:,nodesinbeam(k)) - beam_origin
    !    rproj = dot_product(rcurs,k0)*k0
    !    ell = norm2(rproj)
    !    mu = 1
    !    H = 1/mu
    !
    !    if (ell .lt. H*S) then
    !        cycle
    !    else if (ell .ge. H*S .and. ell .lt. H*(S+tmat1)) then
    !        ells(k,1) = ell - H*S
    !    else if (ell .ge. H*(S+tmat1) .and. ell .lt. H*(S+tmat1+tmat2)) then
    !        ells(k,1) = H*tmat1
    !        ells(k,2) = ell - H*(S+tmat1)
    !    else
    !        ells(k,1) = H*tmat1
    !        ells(k,2) = H*tmat2
    !        ells(k,1) = ells(k,1) + ell - H*(S+tmat1+tmat2)
    !    end if
    !
    !    if (any(ells(k,:) .lt. 0)) then
    !        print *, "manual ray tracing: negative path length found."
    !        stop
    !    end if
    !end do

    !! Planar beam, single material
    !S = beam_origin(3) - Lz
    !do k = 1, size(nodesinbeam)
    !    rcurs = rglobal(:,nodesinbeam(k)) - beam_origin
    !    rproj = dot_product(rcurs,k0)*k0
    !    ell = norm2(rproj)
    !    mu = 1
    !    H = 1/mu
    !
    !    ells(k,1) = ell - H*S
    !    where (ells .lt. 0) ells = 0
    !end do

    ! Spherical beam, single material
    do k = 1, size(nodesinbeam)
        rcurs = rglobal(:,nodesinbeam(k)) - beam_origin
        ell = norm2(rcurs)
        mu = -rcurs(3)/ell
        H = 1/mu

        ells(k,1) = ell - H*S

        where (ells .lt. 0) ells = 0
    end do

    !! Spherical beam, layers
    !do k = 1, size(nodesinbeam)
    !    rcurs = rglobal(:,nodesinbeam(k)) - beam_origin
    !    ell = norm2(rcurs)
    !    mu = -rcurs(3)/ell
    !    H = 1/mu
    !
    !    if (ell .lt. H*S) then
    !        cycle
    !    else if (ell .ge. H*S .and. ell .lt. H*(S+tmat1)) then
    !        ells(k,1) = ell - H*S
    !    else if (ell .ge. H*(S+tmat1) .and. ell .lt. H*(S+tmat1+tmat2)) then
    !        ells(k,1) = H*tmat1
    !        ells(k,2) = ell - H*(S+tmat1)
    !    else if (ell .ge. H*(S+tmat1+tmat2) .and. ell .lt. H*(S+tmat1+tmat2+tmat3)) then
    !        ells(k,1) = H*tmat1
    !        ells(k,2) = H*tmat2
    !        ells(k,3) = ell - H*(S+tmat1+tmat2)
    !    else
    !        ells(k,1) = H*tmat1
    !        ells(k,2) = H*tmat2
    !        ells(k,3) = H*tmat3
    !        ells(k,1) = ells(k,1) + ell - H*(S+tmat1+tmat2+tmat3)
    !    end if
    !end do

    !! GNP
    !do k = 1, size(nodesinbeam)
    !    rcurs = rglobal(:,nodesinbeam(k)) - beam_origin
    !    mu = dot_product(rcurs,beam_coord_sys(:,3))/norm2(rcurs)
    !
    !    if (norm2(rcurs) .eq. abs(rcurs(3))) then
    !        water_entry = (S - R)
    !    else
    !        m = -1/sqrt(1/(mu**2)-1)
    !        d = sqrt(1+m**2-(S/R)**2)
    !        theta1 = atan((m*d+S/R)/(-m*S/R+d))
    !        theta2 = atan((-m*d+S/R)/(-m*S/R-d))
    !        PP1(1) = R*cos(theta1)
    !        PP1(2) = R*sin(theta1) - beam_origin(3)
    !        PP2(1) = R*cos(theta2)
    !        PP2(2) = R*sin(theta2) - beam_origin(3)
    !
    !        water_entry = min(norm2(PP1),norm2(PP2))
    !    end if
    !
    !    water_ell = norm2(rcurs) - water_entry
    !
    !    if (abs(water_ell) .lt. rho) then ! rho is somewhat arbitrarily chosen: I assume that the GNP will be smaller in radius than the mesh density at the surface facing the beam.
    !        ells(k,1) = 0
    !        cycle
    !    end if
    !
    !    if (mu .gt. 1/sqrt(1+(rho/S)**2)) then
    !        if (norm2(rcurs) .eq. abs(rcurs(3))) then
    !            p1 = S - rho - water_entry
    !            p2 = S + rho - water_entry
    !        else
    !            m = -1/sqrt(1/(mu**2)-1)
    !            d = sqrt(1+m**2-(S/rho)**2)
    !            theta1 = atan((m*d+S/rho)/(-m*S/rho+d))
    !            theta2 = atan((-m*d+S/rho)/(-m*S/rho-d))
    !            PP1(1) = rho*cos(theta1) ! ASSUMING BEAM LINE IS STRAIGHT DOWN IN LAB Z AXIS. USES XZ-PLANE (SO ALSO NEEDS SYMMETRY). I probably won't ever need to do otherwise but if I do, work in 3D with PP1 and PP2 as 3D.
    !            PP1(2) = rho*sin(theta1) - beam_origin(3)
    !            PP2(1) = rho*cos(theta2)
    !            PP2(2) = rho*sin(theta2) - beam_origin(3)
    !
    !            p1 = min(norm2(PP1),norm2(PP2)) - water_entry
    !            p2 = max(norm2(PP1),norm2(PP2)) - water_entry
    !        end if
    !
    !        if (water_ell .le. p1) then
    !            ells(k,1) = water_ell
    !        else if (water_ell .le. p2) then
    !            ells(k,1) = p1
    !            ells(k,2) = water_ell - p1
    !        else
    !            ells(k,1) = p1
    !            ells(k,2) = p2 - p1
    !            ells(k,1) = ells(k,1) + water_ell - p2
    !        end if
    !    else
    !        ells(k,1) = water_ell
    !    end if
    !end do
end subroutine TEMP_MANUAL_ray_tracing

subroutine MGXS_discretize_energy_spectrum(E, fEspec, fE, fEg)
    implicit none
    real, dimension(:), intent(in) :: E
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE

    real, dimension(:), allocatable, intent(inout) :: fEg

    integer :: g
    integer :: l_G

    l_G = size(E) - 1
    allocate(fEg(l_G))

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

        fEg(1) = 1.0
    end if
end subroutine MGXS_discretize_energy_spectrum

subroutine FEXS_discretize_energy_spectrum(E, fEspec, fE, fEg)
    implicit none
    real, dimension(:), intent(in) :: E
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE

    real, dimension(:,:), allocatable, intent(inout) :: fEg

    integer :: g, n
    integer :: l_G

    l_G = size(E) - 1
    allocate(fEg(2,l_G))

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
end subroutine FEXS_discretize_energy_spectrum

subroutine SN_external_beam_boundary_source &
    (Cekk, r, eprime, bdyel, normal, khat, RSH, RSH_unc, I1invI2, bdysrc)
    ! Can optimize the storage required by the external beam source.
    ! But would require even more indirect addressing in sweep modules.
    ! I can probably take advantage of esweeplist and esweepbounds, etc., for this purpose
    ! So thoroughly consider that.
    !
    ! At the very least I can make the source not depend on g, and just carry the factors
    ! into the EI subroutines to be applied there.
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: eprime
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:,:,:), intent(in) :: normal
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:), intent(in) :: RSH
    real, dimension(:,:), intent(in) :: RSH_unc
    real, dimension(:,:,:,:), intent(in) :: I1invI2

    real, dimension(:,:,:,:), allocatable, intent(inout) :: bdysrc

    integer :: f, i, ip, j, k
    integer :: NBE
    integer :: NQ
    real, dimension(:), allocatable :: singularity
    real :: dotprod
    integer, dimension(:), allocatable :: invnodesinbeam

    allocate(invnodesinbeam(NK))
    invnodesinbeam = 0
    do k = 1, size(nodesinbeam)
        invnodesinbeam(nodesinbeam(k)) = k
    end do

    NBE = size(bdyel)
    NQ = (NL+1)**2

    allocate(bdysrc(NKe,NBE,Nmu,Nphi))
    allocate(singularity(NKe))

    bdysrc = 0 ! This is for all entries which are never visited
    !$OMP parallel do collapse(3) shared(bdysrc) private(singularity, dotprod) private(f, i, ip, j, k)
    do j = 1, Nphi
        do i = 1, Nmu
            do ip = 1, NBE
                if (.not. any(elsinbeam .eq. bdyel(ip))) cycle
                if (beam_angular_dist .eq. "spherical") then
                    do k = 1, NKe
                        singularity(k) = sum(RSH(i,j,:)*&
                            RSH_unc(invnodesinbeam(Cekk(bdyel(ip),k)),:))
                    end do
                else
                    singularity = sum(RSH(i,j,:)*RSH_unc(1,:)) ! repeated over ip when not necessary, that's fine
                end if
                do f = 1, NFe
                    if (eprime(f,bdyel(ip)) .ne. 0) cycle ! If this isn't a boundary FACE, then cycle
                    dotprod = dot_product(normal(:,f,bdyel(ip)),khat(:,i,j))
                    if (dotprod .ge. 0) cycle ! If this boundary face isn't pointing towards source, then cycle
                    bdysrc(:,ip,i,j) = bdysrc(:,ip,i,j) - &
                        dotprod*matmul(I1invI2(:,:,f,bdyel(ip)),singularity) ! Uses I1invI2 NOT I1invI2f, because e' = 0 means use own element
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    bdysrc = bdysrc/aperture
end subroutine SN_external_beam_boundary_source

subroutine MGXS_construct_uncollided_fluence &
    (label, Cekk, r, ells, Sigmat, fEg, phiu)
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:), intent(in) :: r
    real, dimension(:,:), intent(in) :: ells
    real, dimension(:,:), intent(in) :: Sigmat
    real, dimension(:), intent(in) :: fEg

    real, dimension(:,:,:), allocatable, intent(inout) :: phiu

    integer :: e, k, mat
    integer :: G
    integer, dimension(:), allocatable :: invnodesinbeam

    allocate(invnodesinbeam(NK))
    invnodesinbeam = 0
    do k = 1, size(nodesinbeam)
        invnodesinbeam(nodesinbeam(k)) = k
    end do

    G = size(Sigmat,1)

    allocate(phiu(NE,NKe,G))

    phiu = 0
    if (beam_angular_dist .eq. "spherical") then
        do k = 1, NKe
            do e = 1, size(elsinbeam)
                phiu(elsinbeam(e),k,:) = fEg/&
                    (Euclid(r(1:3,k,elsinbeam(e))-beam_origin)*aperture)
            end do
        end do
    else
        do k = 1, G
            phiu(elsinbeam,1:NKe,k) = fEg(k)/aperture
        end do
    end if

    do k = 1, NKe
        do e = 1, size(elsinbeam)
            do mat = 1, Nmats
                phiu(elsinbeam(e),k,:) = phiu(elsinbeam(e),k,:)*&
                    exp(-Sigmat(:,mat)*ells(invnodesinbeam(Cekk(elsinbeam(e),k)),mat))
            end do
        end do
    end do

    if (any(phiu .lt. 0)) then
        print *, "NEGATIVE"
        stop
    end if

    open(1, file = trim(adjustl(output_fname))//"/"//label//"_uncollided_fluence.dat", form = "unformatted")
    write(1) phiu
    close(1)
end subroutine MGXS_construct_uncollided_fluence

subroutine MGXS_beam_quadrature(Cekk, rglobal, Cfkf, bdyfc, matfc, St, ells, fEg, N, phiu)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cfkf
    integer, dimension(:,:), intent(in) :: bdyfc
    integer, dimension(:,:), intent(in) :: matfc
    real, dimension(:,:), intent(in) :: St
    real, dimension(:,:), intent(in) :: ells
    real, dimension(:), intent(in) :: fEg
    integer, intent(in) :: N

    real, dimension(:,:,:), intent(inout) :: phiu

    integer :: e, k, p
    integer :: G
    real, dimension(:,:), allocatable :: cq
    real, dimension(:), allocatable :: wq
    real, dimension(:,:), allocatable :: I1inv
    real, dimension(:,:), allocatable :: T
    real, dimension(:), allocatable :: tau
    real :: u
    real :: v
    real :: w
    real, dimension(:,:), allocatable :: l_phiu
    integer :: gval
    integer :: eval
    real, dimension(:), allocatable :: rval
    real, dimension(3) :: R
    real, dimension(3) :: kvec
    real, dimension(:,:), allocatable :: l_ells
    real, dimension(:), allocatable :: ll_ells
    integer, dimension(:), allocatable :: invnodesinbeam
    integer :: total = 0

    allocate(invnodesinbeam(NK))
    invnodesinbeam = 0
    do k = 1, size(nodesinbeam)
        invnodesinbeam(nodesinbeam(k)) = k
    end do

    G = size(phiu,3)

    allocate(I1inv(NKe,NKe))
    allocate(T(NKe,N))
    allocate(tau(NKe))
    allocate(l_phiu(N,G))
    allocate(rval(3))
    allocate(l_ells(N,Nmats))
    allocate(ll_ells(Nmats))

    call populate_tetrahedral_quadrature(N, cq, wq)

    I1inv = 20*identity(NKe) - 4

    do p = 1, N
        u = cq(1,p)
        v = cq(2,p)
        w = cq(3,p)

        T(1,p) = wq(p)*(1-u-v-w)
        T(2,p) = wq(p)*u
        T(3,p) = wq(p)*v
        T(4,p) = wq(p)*w
    end do

    ! Choose to examine optical path length criterion of the most populated energy group
    gval = maxloc(fEg, dim = 1)

    R = beam_origin

    do e = 1, size(elsinbeam)
        do k = 1, NKe
            tau(k) = sum(St(gval,:)*ells(invnodesinbeam(Cekk(elsinbeam(e),k)),:))
        end do

        if (minval(tau) .lt. 0.5 .and. maxval(tau) .gt. 0.5) then ! COME UP WITH BETTER CRITERION? IS THIS GOOD?
            eval = elsinbeam(e)
            total = total + 1
            do p = 1, N
                u = cq(1,p)
                v = cq(2,p)
                w = cq(3,p)

                rval = (1-u-v-w)*rglobal(:,Cekk(eval,1)) + &
                               u*rglobal(:,Cekk(eval,2)) + &
                               v*rglobal(:,Cekk(eval,3)) + &
                               w*rglobal(:,Cekk(eval,4))

                l_phiu(p,:) = fEg/&
                    merge(aperture*Euclid(rval-beam_origin), &
                          aperture, &
                          beam_angular_dist .eq. "spherical")

                if (beam_angular_dist .eq. "planar") then
                    kvec = dot_product(rval-beam_origin,k0)*k0
                    R = rval - kvec
                else if (beam_angular_dist .eq. "spherical") then
                    kvec = rval - R
                end if
                if (Nmats .eq. 1) then
                    call single_material_ray_trace_vector &
                        (R, kvec, rglobal, Cfkf, bdyfc, 1.0E-8, l_ells(p,1))
                else
                    call heterogeneous_medium_ray_trace_vector &
                        (R, kvec, rglobal, Cfkf, matfc, 1.0E-8, ll_ells)
                    l_ells(p,:) = ll_ells
                end if
            end do

            do p = 1, N
                l_phiu(p,:) = l_phiu(p,:)*exp(-matmul(St,l_ells(p,:)))
            end do

            phiu(eval,:,:) = 6*matmul(I1inv,matmul(T,l_phiu))
        end if
    end do

    open(1, file = trim(adjustl(output_fname))//"/electron_uncollided_fluence.dat", form = "unformatted")
    write(1) phiu
    close(1)
end subroutine MGXS_beam_quadrature

subroutine MGXS_uncollided_fluence_correction &
    (r, bdyfc, globalnormal, area, vol, Ep, fEspec, fE, Sigmat, phiu)
    implicit none
    real, dimension(:,:,:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: bdyfc
    real, dimension(:,:), intent(in) :: globalnormal
    real, dimension(:,:), intent(in) :: area
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE
    real, dimension(:,:), intent(in) :: Sigmat

    real, dimension(:,:,:), intent(inout) :: phiu

    integer :: e, g, f, i, k, kf
    integer :: l_G
    real, dimension(:), allocatable :: fEg
    integer :: fg
    integer :: kg
    real, dimension(3) :: rcurs
    real, dimension(3) :: kunc
    real :: dotprod
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

        ! Calculate leakage
        if (beam_angular_dist .eq. "spherical") then
            do i = 1, size(bdyfc,1)
                e = bdyfc(i,1)
                f = bdyfc(i,2)
                fg = bdyfc(i,3)
                if (all(phiu(e,1:NKe,g) .eq. 0)) cycle
                do k = 1, NKe
                    if (k .eq. f) cycle
                    rcurs = r(:,k,e) - beam_origin
                    kunc = rcurs/norm2(rcurs)
                    dotprod = dot_product(globalnormal(:,fg),kunc)
                    if (dotprod .le. 0) cycle ! If any node on the face doesn't satisfy k0(r)*n > 0, all nodes will not satisfy it, and the face is pointing TOWARDS source.
                    RL = RL + area(f,e)*dotprod*phiu(e,k,g)/3
                end do
            end do
        else
            kunc = beam_coord_sys(:,3)
            do i = 1, size(bdyfc,1)
                e = bdyfc(i,1)
                f = bdyfc(i,2)
                fg = bdyfc(i,3)
                if (all(phiu(e,1:NKe,g) .eq. 0)) cycle
                do k = 1, NKe
                    if (k .eq. f) cycle
                    dotprod = dot_product(globalnormal(:,fg),kunc)
                    if (dotprod .le. 0) cycle ! If any node on the face doesn't satisfy k0(r)*n > 0, all nodes will not satisfy it, and the face is pointing TOWARDS source.
                    RL = RL + area(f,e)*dotprod*phiu(e,k,g)/3
                    ! SHOULD THIS BE NEGATIVE????
                end do
            end do
        end if

        ! Calculate attenuation
        do k = 1, NKe
            do e = 1, NE
                RT = RT + Sigmat(g,eltomat(e))*vol(e)*phiu(e,k,g)/4
            end do
        end do

        ! Calculate difference between total and leakage (should be ~attenuation)
        RTnew = fEg(g) - RL

        ! Calculate correction factor
        alpha = (RTnew - RT)/(RT**2)
        if (isnan(alpha)) cycle

        ! Correct phiu
        do e = 1, NE
            phiu(e,1:NKe,g) = phiu(e,1:NKe,g)*&
                (1+alpha*vol(e)*Sigmat(g,eltomat(e))*phiu(e,1:NKe,g))
        end do
    end do
end subroutine MGXS_uncollided_fluence_correction

!! UNOPTIMIZED
subroutine SN_MGXS_construct_uncollided_source &
    (Cekk, Pa, factorialmat, Sigma, Pmlk, cmjk, phiu, source)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:), intent(in) :: Pa
    real, dimension(:,:), intent(in) :: factorialmat
    real, dimension(:,:,:,:), intent(in) :: Sigma ! (gp,g,l,mat)
    real, dimension(:,:,:), intent(in) :: Pmlk
    real, dimension(:,:,:), intent(in) :: cmjk
    real, dimension(:,:,:), intent(in) :: phiu ! (e,k,g)

    real, dimension(:,:,:,:), allocatable, intent(inout) :: source ! (k,i,j,g)

    integer :: e, g, i, j, k, kg, l, m, mat, n, np, q
    integer :: NQ
    integer :: l_G
    integer :: l_NE
    integer :: l_k
    integer, dimension(:), allocatable :: ell
    real, dimension(:,:,:), allocatable :: Amat ! (e,k,l+1)
    real, dimension(:,:), allocatable :: tempBmat ! (e,k)
    real, dimension(:,:,:,:,:), allocatable :: Cmat ! (e,k,i,l+1,m+1)
    real, dimension(:,:,:,:), allocatable :: Dmat ! (e,k,i,j)
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

    allocate(source(NK,Nmu,Nphi,l_G))

    source = 0
    do g = 1, l_G
        allocate(Amat(l_NE,NKe,NL+1))
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(2) shared(Amat) private(l, k)
            do l = 0, NL
                do k = 1, NKe
                    Amat(1:l_NE,k,l+1) = &
                        matmul(phiu(elsinbeam,k,1:l_G),Sigma(1:l_G,g,l+1,1))
                end do
            end do
            !$OMP end parallel do
        else
            Amat = 0
            !$OMP parallel do collapse(3) shared(Amat) private(e, k, l)
            do l = 0, NL
                do k = 1, NKe
                    do e = 1, l_NE
                        Amat(e,k,l+1) = &
                            sum(phiu(elsinbeam(e),k,1:l_G)*Sigma(1:l_G,g,l+1,eltomat(elsinbeam(e))))
                    end do
                end do
            end do
            !$OMP end parallel do
        end if

        allocate(tempBmat(l_NE,NKe))
        allocate(Cmat(l_NE,NKe,Nmu,NL+1,NL+1))
        !$OMP parallel do shared(Cmat) private(tempBmat) private(m, l, i) ! CONSIDER USING q SO I CAN PARALLELIZE! ALSO MORE CONVENIENT
        do m = 0, NL
            do l = m, NL
                tempBmat = factorialmat(l+1,m+1)*Amat(:,:,l+1)
                do i = 1, Nmu
                    Cmat(1:l_NE,1:NKe,i,l+1,m+1) = Pa(i,l+1,m+1)*tempBmat
                end do
            end do
        end do
        !$OMP end parallel do
        deallocate(Amat)
        deallocate(tempBmat)

        allocate(Dmat(l_NE,NKe,Nmu,Nphi))
        if (beam_angular_dist .eq. "spherical") then
            !$OMP parallel do collapse(3) shared(Dmat) private(e, i, k, kg)
            do i = 1, Nmu
                do k = 1, NKe
                    do e = 1, l_NE
                        kg = invnodesinbeam(Cekk(elsinbeam(e),k))
                        Dmat(e,k,i,1:Nphi) = matmul(cmjk(1:Nphi,1:NL+1,kg),&
                            sum(Pmlk(1:NL+1,1:NL+1,kg)*&
                            Cmat(e,k,i,1:NL+1,1:NL+1),1))
                    end do
                end do
            end do
            !$OMP end parallel do
        else
            !$OMP parallel do collapse(3) shared(Dmat) private(e, g, i, k)
            do i = 1, Nmu
                do k = 1, NKe
                    do e = 1, l_NE
                        Dmat(e,k,i,1:Nphi) = matmul(cmjk(1:Nphi,1:NL+1,1),&
                            sum(Pmlk(1:NL+1,1:NL+1,1)*&
                            Cmat(e,k,i,1:NL+1,1:NL+1),1))
                    end do
                end do
            end do
            !$OMP end parallel do
        end if
        deallocate(Cmat)

        do k = 1, NKe
            do e = 1, l_NE
                source(Cekk(elsinbeam(e),k),:,:,g) = source(Cekk(elsinbeam(e),k),:,:,g) + &
                    Dmat(e,k,:,:)/Nglobal(Cekk(elsinbeam(e),k))
            end do
        end do
        deallocate(Dmat)
    end do
end subroutine SN_MGXS_construct_uncollided_source

!! UNOPTIMIZED
subroutine PN_MGXS_construct_uncollided_source &
    (l_NL, Cekk, wL, RSH_unc, RSH, Sigma, phiu, source)
    implicit none
    integer, intent(in) :: l_NL
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:), intent(in) :: RSH_unc
    real, dimension(:,:,:), intent(in) :: RSH
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:), intent(in) :: phiu

    real, dimension(:,:,:), allocatable, intent(inout) :: source

    integer :: e, g, ip, jp, k, l, mat, q
    integer :: NQ
    integer :: l_G
    integer :: l_NK
    integer :: l_k
    real, dimension(:,:), allocatable :: l_phiu
    integer, dimension(:), allocatable :: ell

    NQ = (l_NL+1)**2
    l_G = size(Sigma,1)
    l_NK = size(nodesinbeam)

    allocate(ell(NQ))
    ell = [(ceiling(sqrt(real(q)))-1,q=1,NQ)]

    allocate(source(NK,NQ,l_G))
    source = 0

    allocate(l_phiu(NK,l_G))
    l_phiu = 0
    do k = 1, NKe
        do e = 1, NE
            l_phiu(Cekk(e,k),:) = l_phiu(Cekk(e,k),:) + &
                phiu(e,k,:)/Nglobal(Cekk(e,k))
        end do
    end do

    if (beam_angular_dist .eq. "spherical") then
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(2) shared(source) private(g, q)
            do g = 1, l_G
                do q = 1, NQ
                    source(nodesinbeam,q,g) = RSH_unc(1:l_NK,q)*matmul(l_phiu(nodesinbeam,1:g),&
                        Sigma(1:g,g,ell(q)+1,1))/(2*ell(q)+1)
                end do
            end do
            !$OMP end parallel do
        else
            do mat = 1, Nmats
                do k = 1, NKmat(mat)
                    l_k = findloc(nodesinbeam, nodesinmat(k,mat), dim = 1)
                    do g = 1, l_G
                        do q = 1, NQ
                            source(l_k,q,g) = source(l_k,q,g) + &
                                RSH_unc(l_k,q)*Nelnodes(nodesinmat(k,mat),mat)*sum(l_phiu(nodesinmat(k,mat),1:g)*&
                                Sigma(1:g,g,ell(q)+1,mat))/(Nglobal(nodesinmat(k,mat))*(2*ell(q)+1))
                        end do
                    end do
                end do
            end do
        end if
    else
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(2) shared(source) private(g, q)
            do g = 1, l_G
                do q = 1, NQ
                    source(nodesinbeam,q,g) = RSH_unc(1,q)*matmul(l_phiu(nodesinbeam,1:g),&
                        Sigma(1:g,g,ell(q)+1,1))/(2*ell(q)+1)
                end do
            end do
            !$OMP end parallel do
        else
            do mat = 1, Nmats
                do k = 1, NKmat(mat)
                    l_k = findloc(nodesinbeam, nodesinmat(k,mat), dim = 1)
                    do g = 1, l_G
                        do q = 1, NQ
                            source(l_k,q,g) = source(l_k,q,g) + &
                                RSH_unc(1,q)*Nelnodes(nodesinmat(k,mat),mat)*sum(l_phiu(nodesinmat(k,mat),1:g)*&
                                Sigma(1:g,g,ell(q)+1,mat))/(Nglobal(nodesinmat(k,mat))*(2*ell(q)+1))
                        end do
                    end do
                end do
            end do
        end if
    end if

    source = fourpi*source
end subroutine PN_MGXS_construct_uncollided_source

!! UNOPTIMIZED
subroutine SN_to_SN_MGXS_particle_particle_scattering_source &
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
    real, dimension(:,:,:), allocatable, intent(in) :: phiu
    real, dimension(:,:,:,:), intent(in) :: psi

    real, dimension(:,:,:,:), allocatable, intent(inout) :: source

    integer :: e, g, i, ip, jp, k, l, m, mat, q
    integer :: NQ
    integer :: G1
    integer :: G2
    integer :: l_NK
    integer :: l_NE
    integer :: l_k
    integer :: kg
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:), allocatable :: Amat ! (nodes,l+1)
    real, dimension(:,:), allocatable :: tempBmat ! (nodes)
    real, dimension(:,:,:,:,:), allocatable :: Cmat ! (nodes,i,l+1,m+1)
    real, dimension(:,:,:,:), allocatable :: Dmat ! (nodes,i,j)
    real, dimension(:,:,:,:), allocatable :: g_Amat ! (j,i,k,l)
    real, dimension(:,:,:), allocatable :: g_Bmat ! (k,jp,q)
    real, dimension(:,:,:), allocatable :: g_Cmat ! (k,j,q)
    integer, dimension(:), allocatable :: invnodesinbeam

    if (allocated(phiu)) then
        allocate(invnodesinbeam(NK))
        invnodesinbeam = 0
        do k = 1, size(nodesinbeam)
            invnodesinbeam(nodesinbeam(k)) = k
        end do
    end if

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    G1 = size(Sigma,2) ! To
    G2 = size(Sigma,1) ! From
    if (G2 .ne. size(psi,4)) then
        print *, G1, G2
        print *, size(psi,4)
        print *, "ARRAY SIZE DOESNT MATCH"
        stop
    end if
    l_NK = size(nodesinbeam)
    l_NE = size(elsinbeam)

    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5), q=1, NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1, q=1, NQ)]

    allocate(source(NK,Nmu,Nphi,G1))

    source = 0
    if (allocated(phiu)) then
        do g = 1, G1
            allocate(Amat(l_NE,NKe,NL+1))
            if (Nmats .eq. 1) then
                !$OMP parallel do collapse(2) shared(Amat) private(l, k)
                do l = 0, NL
                    do k = 1, NKe
                        Amat(1:l_NE,k,l+1) = &
                        matmul(phiu(elsinbeam,k,1:G2),Sigma(1:G2,g,l+1,1))
                    end do
                end do
                !$OMP end parallel do
            else
                Amat = 0
                !$OMP parallel do collapse(2) shared(Amat) private(e, k)
                do k = 1, NKe
                    do e = 1, l_NE
                        Amat(e,k,1:NL+1) = &
                            matmul(transpose(Sigma(1:G2,g,1:NL+1,eltomat(elsinbeam(e)))),&
                            phiu(elsinbeam(e),k,1:G2))
                    end do
                end do
                !$OMP end parallel do
            end if

            allocate(tempBmat(l_NE,NKe))
            allocate(Cmat(l_NE,NKe,Nmu,NL+1,NL+1))
            !$OMP parallel do shared(Cmat) private(tempBmat) private(m, l, i)
            do m = 0, NL
                do l = m, NL
                    tempBmat = factorialmat(l+1,m+1)*Amat(:,:,l+1)
                    do i = 1, Nmu
                        Cmat(1:l_NE,1:NKe,i,l+1,m+1) = Pa(i,l+1,m+1)*tempBmat
                    end do
                end do
            end do
            !$OMP end parallel do
            deallocate(Amat)
            deallocate(tempBmat)

            allocate(Dmat(l_NE,NKe,Nmu,Nphi))
            if (beam_angular_dist .eq. "spherical") then
                !$OMP parallel do collapse(3) shared(Dmat) private(e, i, k, kg)
                do i = 1, Nmu
                    do k = 1, NKe
                        do e = 1, l_NE
                            kg = invnodesinbeam(Cekk(elsinbeam(e),k))
                            Dmat(e,k,i,:) = matmul(cmjk(1:Nphi,1:NL+1,kg),&
                                sum(Pmlk(1:NL+1,1:NL+1,kg)*&
                                Cmat(e,k,i,1:NL+1,1:NL+1),1))
                        end do
                    end do
                end do
                !$OMP end parallel do
            else
                !$OMP parallel do collapse(3) shared(Dmat) private(e, i, k)
                do i = 1, Nmu
                    do k = 1, NKe
                        do e = 1, l_NE
                            Dmat(e,k,i,:) = matmul(cmjk(1:Nphi,1:NL+1,1),&
                                sum(Pmlk(1:NL+1,1:NL+1,1)*&
                                Cmat(e,k,i,1:NL+1,1:NL+1),1))
                        end do
                    end do
                end do
                !$OMP end parallel do
            end if
            deallocate(Cmat)

            !!$OMP parallel do collapse(2) shared(source) private(e, k)
            do k = 1, NKe
                do e = 1, l_NE
                    source(Cekk(elsinbeam(e),k),:,:,g) = source(Cekk(elsinbeam(e),k),:,:,g) + &
                        Dmat(e,k,:,:)/Nglobal(Cekk(elsinbeam(e),k))
                end do
            end do
            !!$OMP end parallel do

            deallocate(Dmat)
        end do
    end if

    allocate(g_Amat(NK,Nmu,Nphi,NL+1))
    allocate(g_Bmat(NK,Nphi,NQ))
    allocate(g_Cmat(NK,Nphi,NQ))

    do g = 1, G1
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(3) shared(g_Amat) private(l, jp, ip)
            do l = 0, NL
                do jp = 1, Nphi
                    do ip = 1, Nmu
                        g_Amat(1:NK,ip,jp,l+1) = &
                        twopi*matmul(psi(1:NK,ip,jp,1:G2),Sigma(1:G2,g,l+1,1))/Nphi
                    end do
                end do
            end do
            !$OMP end parallel do
        else
            g_Amat = 0
            do l = 0, NL
                do jp = 1, Nphi
                    do ip = 1, Nmu
                        do mat = 1, Nmats
                            do k = 1, NKmat(mat)
                                g_Amat(nodesinmat(k,mat),ip,jp,l+1) = g_Amat(nodesinmat(k,mat),ip,jp,l+1) + &
                                    twopi*Nelnodes(nodesinmat(k,mat),mat)*&
                                    sum(psi(nodesinmat(k,mat),ip,jp,1:G2)*Sigma(1:G2,g,l+1,mat))&
                                    /(Nphi*Nglobal(nodesinmat(k,mat)))
                            end do
                        end do
                    end do
                end do
            end do
        end if

        !$OMP parallel do collapse(2) shared(g_Bmat) private(q, jp)
        do q = 1, NQ
            do jp = 1, Nphi
                g_Bmat(1:NK,jp,q) = qfactorialmat(q)*&
                matmul(g_Amat(1:NK,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do shared(g_Cmat) private(q)
        do q = 1, NQ
            g_Cmat(1:NK,1:Nphi,q) = matmul(g_Bmat(1:NK,1:Nphi,q),cosmat(1:Nphi,1:Nphi,emm(q)+1))
        end do
        !$OMP end parallel do

        !$OMP parallel do shared(source) private(i, q) ! Would generally want this collapse(2), maybe with q then i too. But im not sure how summing conflicts with OMP (~3/31/23)
        do i = 1, Nmu
            do q = 1, NQ
                source(1:NK,i,1:Nphi,g) = source(1:NK,i,1:Nphi,g) + g_Cmat(1:NK,1:Nphi,q)*qPa(i,q)
            end do
        end do
        !$OMP end parallel do
    end do
end subroutine SN_to_SN_MGXS_particle_particle_scattering_source

!! UNOPTIMIZED
subroutine SN_to_SN_FEXS_particle_particle_scattering_source &
    (wL, cosmat, Pa, factorialmat, qPa, qfactorialmat, Sigma, Pmlk, cmjk, psi, source)
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:), intent(in) :: cosmat
    real, dimension(:,:), intent(in) :: qPa
    real, dimension(:), intent(in) :: qfactorialmat
    real, dimension(:,:,:), intent(in) :: Pa
    real, dimension(:,:), intent(in) :: factorialmat
    real, dimension(:,:,:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:), intent(in) :: Pmlk
    real, dimension(:,:,:), intent(in) :: cmjk
    real, dimension(:,:,:,:), intent(in) :: psi

    real, dimension(:,:,:,:,:), allocatable, intent(inout) :: source

    integer :: g, gpr, i, ip, jp, k, l, m, mat, n, np, q
    integer :: NQ
    integer :: G1
    integer :: G2
    integer :: l_NK
    integer :: l_k
    real, dimension(:,:,:), allocatable :: l_phiu
    integer, dimension(:), allocatable :: ell
    integer, dimension(:), allocatable :: emm
    real, dimension(:,:,:,:,:), allocatable :: g_Amat ! (k,n,i,j,l)
    real, dimension(:,:,:,:), allocatable :: g_Bmat ! (k,jp,q,n)
    real, dimension(:,:,:,:), allocatable :: g_Cmat ! (k,j,q,n)

    NQ = NL*(NL+3)/2+1
    allocate(ell(NQ))
    allocate(emm(NQ))
    G1 = size(Sigma,4) ! To
    G2 = size(Sigma,3) ! From
    l_NK = size(nodesinbeam)

    ell = [(floor(0.5*sqrt(8.0*q-7.0)-0.5),q=1,NQ)]
    emm = [(q - (ell(q)*(ell(q)+1))/2 - 1,q=1,NQ)]

    allocate(source(NK,Nmu,Nphi,2,G1))
    source = 0

    allocate(g_Amat(NK,2,Nmu,Nphi,NL+1))
    allocate(g_Bmat(NK,Nphi,NQ,2))
    allocate(g_Cmat(NK,Nphi,NQ,2))

    do g = 1, G1
        g_Amat = 0
        if (Nmats .eq. 1) then
            !$OMP parallel do collapse(4) shared(g_Amat) private(l, gpr, jp, ip, n, np)
            do l = 0, NL
                do jp = 1, Nphi
                    do ip = 1, Nmu
                        do n = 1, 2
                            do np = 1, 2
                                do gpr = 1, G2
                                    g_Amat(1:NK,n,ip,jp,l+1) = g_Amat(1:NK,n,ip,jp,l+1) + &
                                    twopi*psi(1:NK,ip,jp,gpr+np-1)*Sigma(np,n,gpr,g,l+1,1)/Nphi
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        else
            !$OMP parallel do collapse(3) shared(g_Amat) private(l, jp, ip, n, gpr, np, mat, k)
            do l = 0, NL
                do jp = 1, Nphi
                    do ip = 1, Nmu
                        do n = 1, 2
                            do gpr = 1, G2
                                do np = 1, 2
                                    do mat = 1, Nmats
                                        do k = 1, NKmat(mat)
                                            g_Amat(nodesinmat(k,mat),n,ip,jp,l+1) = g_Amat(nodesinmat(k,mat),n,ip,jp,l+1) + &
                                                twopi*Nelnodes(nodesinmat(k,mat),mat)*&
                                                psi(nodesinmat(k,mat),ip,jp,gpr+np-1)*Sigma(np,n,gpr,g,l+1,mat)/&
                                                (Nphi*Nglobal(nodesinmat(k,mat)))
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP end parallel do
        end if

        !$OMP parallel do collapse(2) shared(g_Bmat) private(n, q, jp)
        do n = 1, 2
            do q = 1, NQ
                do jp = 1, Nphi
                    g_Bmat(1:NK,jp,q,n) = matmul(g_Amat(1:NK,n,1:Nmu,jp,ell(q)+1),qPa(1:Nmu,q)*wL(1:Nmu))
                end do
                g_Bmat(1:NK,1:Nphi,q,n) = qfactorialmat(q)*g_Bmat(1:NK,1:Nphi,q,n)
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do collapse(2) shared(g_Cmat) private(g, n, q)
        do n = 1, 2
            do q = 1, NQ
                g_Cmat(1:NK,1:Nphi,q,n) = matmul(g_Bmat(1:NK,1:Nphi,q,n),cosmat(1:Nphi,1:Nphi,emm(q)+1))
            end do
        end do
        !$OMP end parallel do

        !$OMP parallel do shared(source) private(i, q)
        do i = 1, Nmu
            do q = 1, NQ
                source(1:NK,i,1:Nphi,1:2,g) = source(1:NK,i,1:Nphi,1:2,g) + &
                g_Cmat(1:NK,1:Nphi,q,1:2)*qPa(i,q)
            end do
        end do
        !$OMP end parallel do
    end do
end subroutine SN_to_SN_FEXS_particle_particle_scattering_source

!! UNOPTIMIZED
subroutine SN_to_PN_MGXS_particle_particle_scattering_source &
    (l_NL, Cekk, wL, RSH_unc, RSH, Sigma, phiu, psi, source)
    implicit none
    integer, intent(in) :: l_NL
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:), intent(in) :: RSH_unc
    real, dimension(:,:,:), intent(in) :: RSH
    real, dimension(:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:), intent(in) :: psi
    real, dimension(:,:,:), allocatable, intent(in) :: phiu

    real, dimension(:,:,:), allocatable, intent(inout) :: source

    integer :: e, g, ip, jp, k, mat, q
    integer :: NQ
    integer :: G1
    integer :: G2
    integer :: l_NK
    integer :: l_k
    real, dimension(:,:), allocatable :: l_phiu
    integer, dimension(:), allocatable :: ell

    NQ = (l_NL+1)**2
    G1 = size(Sigma,2) ! To
    G2 = size(Sigma,1) ! From
    l_NK = size(nodesinbeam)

    allocate(source(NK,NQ,G1))
    allocate(ell(NQ))
    ell = [(ceiling(sqrt(real(q)))-1,q=1,NQ)]

    source = 0

    if (allocated(phiu)) then
        allocate(l_phiu(NK,G2))
        l_phiu = 0
        do k = 1, NKe
            do e = 1, NE
                l_phiu(Cekk(e,k),:) = l_phiu(Cekk(e,k),:) + &
                    phiu(e,k,:)/Nglobal(Cekk(e,k))
            end do
        end do

        if (Nmats .eq. 1) then
            if (beam_angular_dist .eq. "spherical") then
                !$OMP parallel do collapse(2) shared(source) private(g, q)
                do g = 1, G1
                    do q = 1, NQ
                        source(nodesinbeam,q,g) = RSH_unc(1:l_NK,q)*matmul(l_phiu(nodesinbeam,1:G2),&
                            Sigma(1:G2,g,ell(q)+1,1))/(2*ell(q)+1)
                    end do
                end do
                !$OMP end parallel do
            else
                !$OMP parallel do collapse(2) shared(source) private(g, q)
                do g = 1, G1
                    do q = 1, NQ
                        source(nodesinbeam,q,g) = RSH_unc(1,q)*matmul(l_phiu(nodesinbeam,1:G2),&
                            Sigma(1:G2,g,ell(q)+1,1))/(2*ell(q)+1)
                    end do
                end do
                !$OMP end parallel do
            end if
        else
            if (beam_angular_dist .eq. "spherical") then
                do mat = 1, Nmats
                    do k = 1, NKmat(mat)
                        l_k = findloc(nodesinbeam, nodesinmat(k,mat), dim = 1)
                        !$OMP parallel do collapse(2) shared(source) private(g, k, mat, q)
                        do g = 1, G1
                            do q = 1, NQ
                                source(nodesinmat(k,mat),q,g) = source(nodesinmat(k,mat),q,g) + &
                                    Nelnodes(nodesinmat(k,mat),mat)*&
                                    RSH_unc(l_k,q)*sum(l_phiu(nodesinmat(k,mat),1:G2)*&
                                    Sigma(1:G2,g,ell(q)+1,mat))/((2*ell(q)+1)*Nglobal(nodesinmat(k,mat)))
                            end do
                        end do
                        !$OMP end parallel do
                    end do
                end do
            else
                do mat = 1, Nmats
                    do k = 1, NKmat(mat)
                        !$OMP parallel do collapse(2) shared(source) private(g, k, mat, q)
                        do g = 1, G1
                            do q = 1, NQ
                                source(nodesinmat(k,mat),q,g) = source(nodesinmat(k,mat),q,g) + &
                                    Nelnodes(nodesinmat(k,mat),mat)*&
                                    RSH_unc(1,q)*sum(l_phiu(nodesinmat(k,mat),1:G2)*&
                                    Sigma(1:G2,g,ell(q)+1,mat))/((2*ell(q)+1)*Nglobal(nodesinmat(k,mat)))
                            end do
                        end do
                        !$OMP end parallel do
                    end do
                end do
            end if
        end if
    end if

    !$OMP parallel do collapse(2) shared(source) private(g, q, ip, jp)
    do g = 1, G1
        do q = 1, NQ
            do jp = 1, Nphi
                do ip = 1, Nmu ! Can prevent having to do the psi,Sigma matmul over ell(q) multiple times by outer iterating over l then inner over m.
                    source(1:NK,q,g) = source(1:NK,q,g) + &
                        twopi*wL(ip)*RSH(ip,jp,q)*&
                        matmul(psi(1:NK,ip,jp,1:G2),Sigma(1:G2,g,ell(q)+1,1))/(Nphi*(2*ell(q)+1))
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    source = fourpi*source
end subroutine SN_to_PN_MGXS_particle_particle_scattering_source

!! UNOPTIMIZED
subroutine SN_to_PN_FEXS_particle_particle_scattering_source &
    (l_NL, wL, RSH_unc, RSH, Sigma, psi, source)
    implicit none
    integer, intent(in) :: l_NL
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:), intent(in) :: RSH_unc
    real, dimension(:,:,:), intent(in) :: RSH
    real, dimension(:,:,:,:,:,:), intent(in) :: Sigma
    real, dimension(:,:,:,:), intent(in) :: psi

    real, dimension(:,:,:,:), allocatable, intent(inout) :: source

    integer :: g, gpr, ip, jp, k, mat, n, np, q
    integer :: NQ
    integer :: G1
    integer :: G2
    integer :: l_NK
    integer :: l_k
    integer, dimension(:), allocatable :: ell

    NQ = (l_NL+1)**2
    G1 = size(Sigma,4) ! To
    G2 = size(Sigma,3) ! From
    l_NK = size(nodesinbeam)

    allocate(source(NK,NQ,2,G1))
    allocate(ell(NQ))
    ell = [(ceiling(sqrt(real(q)))-1,q=1,NQ)]

    source = 0

    !$OMP parallel do collapse(3) shared(source) private(g, q, ip, jp, np, n, gpr)
    do g = 1, G1
        do n = 1, 2
            do q = 1, NQ
                do gpr = 1, G2
                    do np = 1, 2
                        do jp = 1, Nphi
                            do ip = 1, Nmu
                                source(1:NK,q,n,g) = source(1:NK,q,n,g) + &
                                    twopi*wL(ip)*RSH(ip,jp,q)*&
                                    psi(1:NK,ip,jp,gpr+np-1)*Sigma(np,n,gpr,g,ell(q)+1,1)/(Nphi*(2*ell(q)+1))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    source = fourpi*source
end subroutine SN_to_PN_FEXS_particle_particle_scattering_source

subroutine SN_construct_T_MMS_vector(function, rglobal, mu, phiC, vector) ! Do a MMS solve for the uncollided beam?
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: phiC

    real, dimension(:,:,:,:), allocatable, intent(inout) :: vector ! (k,i,j,g)

    integer :: i, j, k
    real :: x
    real :: y
    real :: z
    real :: l_mu
    real :: phi
    interface
        real function function(fx,fy,fz,fmu,fphi)
            implicit none
            real, intent(in) :: fx
            real, intent(in) :: fy
            real, intent(in) :: fz
            real, intent(in) :: fmu
            real, intent(in) :: fphi
        end function function
    end interface

    allocate(vector(NK,Nmu,Nphi,1))

    !$OMP parallel do collapse(3) shared(vector) private(i, j, k) private(x, y, z, l_mu, phi)
    do j = 1, Nphi
        do i = 1, Nmu
            do k = 1, NK
                x = rglobal(1,k)
                y = rglobal(2,k)
                z = rglobal(3,k)
                l_mu = mu(i)
                phi = phiC(j)

                vector(k,i,j,1) = function(x,y,z,l_mu,phi)
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine SN_construct_T_MMS_vector

subroutine SN_construct_full_MMS_vectors(rglobal, mu, wL, phiC, RSH, Sigma, vector, source)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:), intent(in) :: phiC
    real, dimension(:,:,:), intent(in) :: RSH

    real, dimension(:,:,:,:), allocatable, intent(inout) :: Sigma ! (gpr,g,l,mat)
    real, dimension(:,:,:,:), allocatable, intent(inout) :: vector ! (k,i,j,g)
    real, dimension(:,:,:,:), allocatable, intent(inout) :: source ! (k,i,j,g)

    integer :: i, j, k, l, q
    integer :: NQ
    integer, dimension(:), allocatable :: ell
    real :: x
    real :: y
    real :: z
    real :: l_mu
    real :: phi
    real, dimension(:,:), allocatable :: P ! (i,l)
    real, dimension(:), allocatable :: Sigmai
    real, dimension(:,:,:), allocatable :: q_vector ! (k,q,g)

    NQ = (NL+1)**2
    allocate(Sigma(1,1,NL+1,1))
    allocate(vector(NK,Nmu,Nphi,1))
    allocate(source(NK,Nmu,Nphi,1))
    allocate(ell(NQ))
    allocate(P(Nmu,0:NL))
    allocate(Sigmai(Nmu))
    allocate(q_vector(NK,NQ,1))

    ell = [(ceiling(sqrt(real(q)))-1,q=1,NQ)]

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

    !$OMP parallel do collapse(3) shared(vector) private(i, j, k) private(x, y, z, l_mu, phi)
    do j = 1, Nphi
        do i = 1, Nmu
            do k = 1, NK
                x = rglobal(1,k)
                y = rglobal(2,k)
                z = rglobal(3,k)
                l_mu = mu(i)
                phi = phiC(j)

                vector(k,i,j,1) = full_MMS_solution(x,y,z,l_mu,phi)
                source(k,i,j,1) = full_MMS_non_K_source(x,y,z,l_mu,phi)
            end do
        end do
    end do
    !$OMP end parallel do

    q_vector = 0
    !$OMP parallel do shared(q_vector) private(j, q)
    do q = 1, NQ
        do j = 1, Nphi
            q_vector(1:NK,q,1) = q_vector(1:NK,q,1) + &
                twopi*matmul(vector(1:NK,1:Nmu,j,1),wL*RSH(1:Nmu,j,q))/Nphi
        end do
    end do
    !$OMP end parallel do

    do k = 1, NK
        do q = 1, NQ
            source(k,1:Nmu,1:Nphi,1) = source(k,1:Nmu,1:Nphi,1) - &
                fourpi*Sigma(1,1,ell(q)+1,1)*RSH(1:Nmu,1:Nphi,q)*q_vector(k,q,1)/(2*ell(q)+1)
        end do
    end do
end subroutine SN_construct_full_MMS_vectors

subroutine SN_boundary_source &
    (Cekk, r, eprime, bdyel, normal, mu, phiC, khat, I1invI2, bdysrc)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: eprime
    integer, dimension(:), intent(in) :: bdyel
    real, dimension(:,:,:), intent(in) :: normal
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: phiC
    real, dimension(:,:,:), intent(in) :: khat
    real, dimension(:,:,:,:), intent(in) :: I1invI2

    real, dimension(:,:,:,:), allocatable, intent(inout) :: bdysrc

    integer :: f, i, ip, j, k
    integer :: NBE
    integer :: NQ
    real :: dotprod
    real, dimension(:,:,:,:), allocatable :: func
    real :: x
    real :: y
    real :: z
    real :: l_mu
    real :: l_phi

    NBE = size(bdyel)
    NQ = (NL+1)**2

    allocate(bdysrc(NKe,NBE,Nmu,Nphi))
    allocate(func(NKe,NBE,Nmu,Nphi))

    do j = 1, Nphi
        l_phi = phiC(j)
        do i = 1, Nmu
            l_mu = mu(i)
            do ip = 1, NBE
                do k = 1, NKe
                    x = r(1,k,bdyel(ip))
                    y = r(2,k,bdyel(ip))
                    z = r(3,k,bdyel(ip))
                    func(k,ip,i,j) = T_MMS_solution(x,y,z,l_mu,l_phi)
                end do
            end do
        end do
    end do

    bdysrc = 0 ! This is for all entries which are never visited
    !$OMP parallel do collapse(3) shared(bdysrc) private(dotprod) private(f, i, ip, j, k)
    do j = 1, Nphi
        do i = 1, Nmu
            do ip = 1, NBE
                do f = 1, NFe
                    if (eprime(f,bdyel(ip)) .ne. 0) cycle
                    dotprod = dot_product(normal(:,f,bdyel(ip)),khat(:,i,j))
                    if (dotprod .ge. 0) cycle
                    bdysrc(:,ip,i,j) = bdysrc(:,ip,i,j) - &
                        dotprod*matmul(I1invI2(:,:,f,bdyel(ip)),func(:,ip,i,j)) ! Uses I1invI2 NOT I1invI2f, because e' = 0 means use own element
                end do
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine SN_boundary_source

subroutine PN_construct_T_MMS_vector(function, rglobal, Cekk, mu, wL, phiC, RSH, vector)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:), intent(in) :: phiC
    real, dimension(:,:,:), intent(in) :: RSH

    real, dimension(:,:,:,:), allocatable, intent(inout) :: vector ! (e,k,q,g)

    integer :: e, i, j, k, q
    integer :: NQ
    real, dimension(:,:,:), allocatable :: g_vector ! (k,q,g)
    real :: x
    real :: y
    real :: z
    real :: l_mu
    real :: phi
    real, dimension(:,:,:), allocatable :: t_RSH
    interface
        real function function(fx,fy,fz,fmu,fphi)
            implicit none
            real, intent(in) :: fx
            real, intent(in) :: fy
            real, intent(in) :: fz
            real, intent(in) :: fmu
            real, intent(in) :: fphi
        end function function
    end interface

    NQ = size(RSH,3)

    allocate(g_vector(NK,NQ,1))
    g_vector = 0

    !$OMP parallel do collapse(2) shared(g_vector) private(i, j, k, q) private(x, y, z, l_mu, phi)
    do q = 1, NQ
        do k = 1, NK
            x = rglobal(1,k)
            y = rglobal(2,k)
            z = rglobal(3,k)
            do j = 1, Nphi
                phi = phiC(j)
                do i = 1, Nmu
                    l_mu = mu(i)
                    g_vector(k,q,1) = g_vector(k,q,1) + twopi*wL(i)*RSH(i,j,q)*function(x,y,z,l_mu,phi)/Nphi
                end do
            end do
        end do
    end do
    !$OMP end parallel do

    allocate(vector(NE,NKe,NQ,1))
    !$OMP parallel do collapse(2) shared(vector) private(e, k)
    do k = 1, NKe
        do e = 1, NE
            vector(e,k,1:NQ,1) = g_vector(Cekk(e,k),1:NQ,1)
        end do
    end do
    !$OMP end parallel do
end subroutine PN_construct_T_MMS_vector

end module sources