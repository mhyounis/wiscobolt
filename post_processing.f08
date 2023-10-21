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

module post_processing
    use math
    use physics
    use user_input
implicit none

contains
subroutine SN_construct_fluence(wL, psi, phi)
    ! If going to keep using, need to make phi elemental
    implicit none
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:), intent(in) :: psi

    real, dimension(:,:), allocatable, intent(inout) :: phi

    integer :: g, i, j, k
    integer :: l_G
    integer :: l_Nmu
    integer :: l_Nnu

    l_G = size(psi,4)
    l_Nmu = size(psi,2)
    if (l_Nmu .ne. size(wL)) then
        print *, l_Nmu
        print *, size(wL)
        print *, "STOPPED - SN_construct fluence - message wip"
        stop
    end if
    l_Nnu = size(psi,3)

    allocate(phi(NK,l_G))

    phi = 0

    do g = 1, l_G
        do j = 1, l_Nnu
            phi(1:NK,g) = phi(1:NK,g) + twopi*matmul(psi(1:NK,1:l_Nmu,j,g),wL)/l_Nnu
        end do
    end do
end subroutine SN_construct_fluence

subroutine PN_construct_fluence(psi, phi)
    implicit none
    real, dimension(:,:,:), intent(in) :: psi

    real, dimension(:,:), allocatable, intent(inout) :: phi

    allocate(phi, source = sqrt(fourpi)*psi(:,1,:))
end subroutine PN_construct_fluence

subroutine homogenize_DFEM(Cekk, DFEM, FEM)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:), intent(in) :: DFEM

    real, dimension(:), intent(inout) :: FEM ! Must be pre-allocated to [NK]

    integer :: e, k

    FEM = 0
    do k = 1, NKe
        do e = 1, NE
            FEM(Cekk(e,k)) = FEM(Cekk(e,k)) + &
                DFEM(e,k)/Nglobal(Cekk(e,k))
        end do
    end do
end subroutine homogenize_DFEM

subroutine construct_mesh_slice_overlay &
    (label, dim, slice, rglobal, Cfkf, slicek, slicerglobal, sliceCfkf)
    implicit none
    character(*), intent(in) :: label
    integer, intent(in) :: dim
    real, intent(in) :: slice
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cfkf

    integer, dimension(:), allocatable, intent(inout) :: slicek
    real, dimension(:,:), allocatable, intent(inout) :: slicerglobal
    integer, dimension(:,:), allocatable, intent(inout) :: sliceCfkf

    integer :: f, i, j, k
    logical, dimension(:), allocatable :: nodemask
    integer, dimension(:), allocatable :: slicekinv
    logical, dimension(:), allocatable :: facemask
    integer, dimension(:), allocatable :: keptfaces
    character(50) :: fname

    allocate(nodemask(NK))
    allocate(slicekinv(NK))
    allocate(facemask(size(Cfkf,1)))

    nodemask = abs(rglobal(dim,1:NK) - slice) .le. 1.0E-15

    allocate(slicek, source = pack([(k,k=1,NK)], nodemask))
    allocate(slicerglobal(2,size(slicek)))

    if (dim .eq. 1) then
        slicerglobal(1,1:size(slicek)) = rglobal(2,slicek)
        slicerglobal(2,1:size(slicek)) = rglobal(3,slicek)
    else if (dim .eq. 2) then
        slicerglobal(1,1:size(slicek)) = rglobal(1,slicek)
        slicerglobal(2,1:size(slicek)) = rglobal(3,slicek)
    else
        slicerglobal(1,1:size(slicek)) = rglobal(1,slicek)
        slicerglobal(2,1:size(slicek)) = rglobal(2,slicek)
    end if

    slicekinv = 0
    do k = 1, size(slicek)
        slicekinv(slicek(k)) = k
    end do
    ! slicekinv takes a global k and gives its address in the slice mesh (if there is one)
    ! Zero is returned if there is no corresponding node

    do f = 1, size(Cfkf,1)
        facemask(f) = any(Cfkf(f,1) .eq. slicek) .and. &
                      any(Cfkf(f,2) .eq. slicek) .and. &
                      any(Cfkf(f,3) .eq. slicek)
    end do

    allocate(keptfaces, source = pack([(f,f=1,size(Cfkf,1))], facemask))
    allocate(sliceCfkf(size(keptfaces),3))

    do f = 1, size(keptfaces)
        sliceCfkf(f,1:3) = Cfkf(keptfaces(f),1:3)
        do k = 1, NKf
            sliceCfkf(f,k) = slicekinv(sliceCfkf(f,k))
        end do
    end do

    if (results_as_dat) then
        fname = trim(adjustl(output_fname))//"/nodes_"//label//".dat"

        open(1, form = "unformatted", file = fname, action = "write")
        write(1) slicerglobal
        close(1)

        fname = trim(adjustl(output_fname))//"/faces_"//label//".dat"

        open(1, form = "unformatted", file = fname, action = "write")
        write(1) sliceCfkf
        close(1)
    else
        fname = "Results/nodes_"//label//".txt"

        open(1, file = fname, status = "replace")
        do j = 1, size(slicerglobal,2)
            do i = 1, size(slicerglobal,1)
                write(1,*) slicerglobal(i,j)
            end do
        end do
        close(1)

        fname = "Results/faces_"//label//".txt"
        open(1, file = fname, status = "replace")
        do j = 1, size(sliceCfkf,2)
            do i = 1, size(sliceCfkf,1)
                write(1, '(I0)') sliceCfkf(i,j)
            end do
        end do
        close(1)
    end if

    print *, "----------------------------------------------------------"
    print *, "### POST-PROCESSING. Number of nodes, faces in mesh slice."
    print *, size(slicek), size(keptfaces)
    fname = trim(adjustl(output_fname))//"/nodes_"//label//".dat"
    print *, "Nodes saved as:", fname
    fname = trim(adjustl(output_fname))//"/faces_"//label//".dat"
    print *, "Faces saved as:", fname
    print *, "----------------------------------------------------------"
end subroutine construct_mesh_slice_overlay

subroutine visualize_mesh_slice_fluence(label, Cekk, slicek, phi)
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: slicek
    real, dimension(:,:), intent(in) :: phi

    integer :: i, j, k
    integer :: l_G
    real, dimension(:,:), allocatable :: slicephi
    !integer, dimension(:), allocatable :: slicekinv
    character(80) :: fname

    l_G = size(phi,2)

    allocate(slicephi(l_G,size(slicek)))

    do k = 1, size(slicek)
        slicephi(1:l_G,k) = phi(slicek(k),1:l_G)
    end do

    if (results_as_dat) then
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_fluence.dat"

        open(1, file = fname, form = "unformatted")
        write(1) slicephi
        close(1)
    else
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_fluence.txt"

        open(1, file = fname, status = "replace")
        do j = 1, size(slicek) ! Integers k and f have no meaning
            do i = 1, l_G
                write(1,*) slicephi(i,j)
            end do
        end do
        close(1)
    end if

    print *, "----------------------------------------"
    print *, "### POST-PROCESSING: Mesh slice fluence."
    print *, "Saved as:", fname
    print *, "----------------------------------------"
end subroutine visualize_mesh_slice_fluence

subroutine SN_visualize_mesh_slice_angular_fluence(label, Cekk, slicek, psi)
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: slicek
    real, dimension(:,:,:,:), intent(in) :: psi

    integer :: g, i, j, k
    integer :: l_G
    real, dimension(:,:,:,:), allocatable :: slicepsi
    !integer, dimension(:), allocatable :: slicekinv
    character(80) :: fname

    l_G = size(psi,4)

    allocate(slicepsi(l_G,Nmu,Nphi,size(slicek)))

    do k = 1, size(slicek)
        do g = 1, l_G
            slicepsi(g,1:Nmu,1:Nphi,k) = psi(slicek(k),1:Nmu,1:Nphi,g)
        end do
    end do

    if (results_as_dat) then
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_angular_fluence.dat"

        open(1, file = fname, form = "unformatted")
        write(1) slicepsi
        close(1)
    else
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_angular_fluence.txt"

        open(1, file = fname, status = "replace")
        do k = 1, size(slicek) ! Integers k and f have no meaning
            do j = 1, Nphi
                do i = 1, Nmu
                    do g = 1, l_G
                        write(1,*) slicepsi(g,i,j,k)
                    end do
                end do
            end do
        end do
        close(1)
    end if

    print *, "------------------------------------------------"
    print *, "### POST-PROCESSING: Mesh slice angular fluence."
    print *, "Saved as:", fname
    print *, "------------------------------------------------"
end subroutine SN_visualize_mesh_slice_angular_fluence

subroutine PN_visualize_mesh_slice_angular_fluence(label, Cekk, slicek, psi)
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: slicek
    real, dimension(:,:,:), intent(in) :: psi

    integer :: g, q, k
    integer :: l_G
    integer :: NQ
    real, dimension(:,:,:), allocatable :: slicepsi
    !integer, dimension(:), allocatable :: slicekinv
    character(80) :: fname

    l_G = size(psi,3)
    NQ = size(psi,2)

    allocate(slicepsi(l_G,NQ,size(slicek)))

    do k = 1, size(slicek)
        do g = 1, l_G
            slicepsi(g,1:NQ,k) = psi(slicek(k),1:NQ,g)
        end do
    end do

    if (results_as_dat) then
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_angular_fluence.dat"

        open(1, file = fname, form = "unformatted")
        write(1) slicepsi
        close(1)
    else
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_angular_fluence.txt"

        open(1, file = fname, status = "replace")
        do k = 1, size(slicek) ! Integers k and f have no meaning
            do q = 1, NQ
                do g = 1, l_G
                    write(1,*) slicepsi(g,q,k)
                end do
            end do
        end do
        close(1)
    end if

    print *, "----------------------------------------------------------"
    print *, "### POST-PROCESSING: Mesh slice angular fluence file name."
    print *, fname
    print *, "----------------------------------------------------------"
end subroutine PN_visualize_mesh_slice_angular_fluence

subroutine construct_line(label, const_axes, const_vals, rglobal, psi)
    ! MUST MAKE THIS WORK FOR NOT JUST STRAIGHT DOWN. ROTATE RGLOBAL AND DO STUFF. FOR NOW, KEEP AS IS.
    ! For z axis pointing down: const_axes = [1,2]
    ! const_vals = [0.0,0.0]
    ! So rotate + translate mesh to system where beam axis is always z-axis and beam origin is 0.0, 0.0. EASY
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:), intent(in) :: const_axes
    real, dimension(:), intent(in) :: const_vals
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:), intent(in) :: psi

    integer :: k
    logical, dimension(:), allocatable :: nodemask
    integer, dimension(:), allocatable :: linek
    !real, dimension(:,:), allocatable :: ROTrglobal
    !real, dimension(:,:), allocatable :: Rot
    !real, dimension(:), allocatable :: transl
    character(100) :: printoutname

    !allocate(ROTrglobal())

    allocate(nodemask(NK))

    nodemask = abs(rglobal(const_axes(1),1:NK) - const_vals(1)) < 1.0E-10 &
        .and. abs(rglobal(const_axes(2),1:NK) - const_vals(2)) < 1.0E-10

    allocate(linek, source = pack([(k,k=1,NK)], mask = nodemask))

    print *, "------------------------------------------------"
    print *, "### POST-PROCESSING: Line nodes printout."
    print *, "LABEL:", label
    do k = 1, size(linek)
        print *, rglobal(6-const_axes(2)-const_axes(1),linek(k))
    end do
    print *, "### POST-PROCESSING: Line distribution printout."
    print *, "LABEL:", label
    do k = 1, size(linek)
        print *, psi(linek(k))
    end do
    print *, "------------------------------------------------"
end subroutine construct_line

subroutine DFEM_construct_line(label, const_axes, const_vals, r, psi)
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:), intent(in) :: const_axes
    real, dimension(:), intent(in) :: const_vals
    real, dimension(:,:,:), intent(in) :: r ! (dir,k,e)
    real, dimension(:,:), intent(in) :: psi ! (e,k)

    integer :: counter, e, i, k
    logical, dimension(:,:), allocatable :: nodemask
    integer, dimension(:,:), allocatable :: lineek
    character(100) :: printoutname

    allocate(nodemask(NKe,NE))

    nodemask = abs(r(const_axes(1),1:NKe,1:NE) - const_vals(1)) .lt. 1.0E-10 &
        .and. abs(r(const_axes(2),1:NKe,1:NE) - const_vals(2)) .lt. 1.0E-10

    allocate(lineek(count(nodemask),2))
    counter = 0
    do k = 1, NKe
        do e = 1, NE
            if (nodemask(k,e)) then
                counter = counter + 1
                lineek(counter,1) = k
                lineek(counter,2) = e
            end if
        end do
    end do

    print *, "------------------------------------------------"
    print *, "### POST-PROCESSING: Line nodes printout."
    print *, "LABEL:", label
    do i = 1, size(lineek,1)
        print *, r(6-const_axes(2)-const_axes(1),lineek(i,1),lineek(i,2))
    end do
    print *, "### POST-PROCESSING: Line distribution printout."
    print *, "LABEL:", label
    do i = 1, size(lineek,1)
        print *, psi(lineek(i,2),lineek(i,1))
    end do
    print *, "------------------------------------------------"
end subroutine DFEM_construct_line

subroutine SN_MGXS_polar_angle_distribution(particle_type, Cekk, vol, mu, wL, Ep, psi)
    implicit none
    character(*), intent(in) :: particle_type
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:), intent(in) :: psi

    integer :: e, g, i, k
    integer :: l_G
    real, dimension(:,:,:), allocatable :: int_space
    real, dimension(:,:), allocatable :: int_azimuth
    real, dimension(:), allocatable :: int_energy
    real, dimension(:), allocatable :: invbetas
    character(100) :: printoutname
    real :: norm

    l_G = size(psi,4)
    allocate(int_space(Nmu,Nphi,l_G))
    allocate(int_azimuth(Nmu,l_G))
    allocate(int_energy(Nmu))
    allocate(invbetas(l_G))

    int_space = 0
    do k = 1, NKe
        do e = 1, NE
            int_space = int_space + vol(e)*psi(Cekk(e,k),1:Nmu,1:Nphi,1:l_G)/4
        end do
    end do

    int_azimuth = twopi*sum(int_space,2)/Nphi

    if (particle_type .eq. "electron") then
        do g = 1, l_G
            invbetas(g) = e_mass_E*(sqrt((Ep(g)/e_mass_E)*(Ep(g)/e_mass_E+2)) - &
                sqrt((Ep(g+1)/e_mass_E)*(Ep(g+1)/e_mass_E+2)))/(Ep(g)-Ep(g+1))
        end do
    else
        invbetas = 1
    end if

    int_energy = matmul(int_azimuth,invbetas)/(4*lightspeed)

    norm = sum(wL*int_energy)

    int_energy = int_energy/norm

    if (particle_type .eq. "photon") then
        printoutname = "LABEL: photon"
    else if (particle_type .eq. "electron") then
        printoutname = "LABEL: electron"
    else
        printoutname = "LABEL: "//trim(adjustl(particle_type))
    end if

    print *, "-------------------------------------------------------"
    print *, "### POST-PROCESSING: Polar angle abscissae printout."
    print *, printoutname
    do i = 1, Nmu
        print *, mu(i)
    end do
    print *, "### POST-PROCESSING: Polar angle distribution printout."
    print *, printoutname
    do i = 1, Nmu
        print *, int_energy(i)
    end do
    print *, "-------------------------------------------------------"
end subroutine SN_MGXS_polar_angle_distribution

subroutine SN_FEXS_polar_angle_distribution(particle_type, Cekk, vol, mu, wL, Ep, psi)
    implicit none
    character(*), intent(in) :: particle_type
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: Ep
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:,:,:,:), intent(in) :: psi

    integer :: e, g, i, k
    integer :: l_G
    real, dimension(:,:,:), allocatable :: int_space
    real, dimension(:,:), allocatable :: int_azimuth
    real, dimension(:), allocatable :: int_energy
    real, dimension(:), allocatable :: invbetas
    character(100) :: printoutname
    real :: norm

    l_G = size(psi,4) - 1
    allocate(int_space(Nmu,Nphi,l_G+1))
    allocate(int_azimuth(Nmu,l_G+1))
    allocate(int_energy(Nmu))
    allocate(invbetas(l_G+1))

    int_space = 0
    do k = 1, NKe
        do e = 1, NE
            int_space = int_space + vol(e)*psi(Cekk(e,k),1:Nmu,1:Nphi,1:l_G+1)/4
        end do
    end do

    int_azimuth = twopi*sum(int_space,2)/Nphi

    if (particle_type .eq. "electron") then
        do g = 1, l_G+1
            invbetas(g) = 1/beta(Ep(g))
        end do
    else
        invbetas = 1
    end if

    int_energy = 0
    do g = 1, l_G
        int_energy = int_energy + 0.5*(int_azimuth(:,g)*invbetas(g) + int_azimuth(:,g+1)*invbetas(g+1))*&
        (Ep(g)-Ep(g+1))/(4*lightspeed)
    end do

    norm = sum(wL*int_energy)

    int_energy = int_energy/norm

    if (particle_type .eq. "photon") then
        printoutname = "LABEL: photon"
    else
        printoutname = "LABEL: electron"
    end if

    print *, "-------------------------------------------------------"
    print *, "### POST-PROCESSING: Polar angle abscissae printout."
    print *, printoutname
    do i = 1, Nmu
        print *, mu(i)
    end do
    print *, "### POST-PROCESSING: Polar angle distribution printout."
    print *, printoutname
    do i = 1, Nmu
        print *, int_energy(i)
    end do
    print *, "-------------------------------------------------------"
end subroutine SN_FEXS_polar_angle_distribution

subroutine MGXS_energy_spectrum(particle_type, label, Cekk, vol, Ep, phi)
    implicit none
    character(*), intent(in) :: particle_type
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: Ep
    real, dimension(:,:), intent(in) :: phi

    integer :: e, g, k, n
    integer :: l_G
    real, dimension(:), allocatable :: int_space
    real, dimension(:), allocatable :: invbetas
    real, dimension(:), allocatable :: ints
    real :: norm
    character(100) :: printoutname

    l_G = size(Ep)-1 ! Use size of phi to inform whether FEXS or MGXS????

    allocate(int_space(l_G))
    allocate(invbetas(l_G))
    allocate(ints(l_G))

    int_space = 0
    do k = 1, NKe
        do e = 1, NE
            int_space = int_space + vol(e)*phi(Cekk(e,k),1:l_G)/4
        end do
    end do

    if (particle_type .eq. "electron") then
        do g = 1, l_G
            invbetas(g) = e_mass_E*(sqrt((Ep(g)/e_mass_E)*(Ep(g)/e_mass_E+2)) - &
                sqrt((Ep(g+1)/e_mass_E)*(Ep(g+1)/e_mass_E+2)))/(Ep(g)-Ep(g+1))
        end do
    else
        invbetas = 1
    end if

    do g = 1, l_G
        ints(g) = int_space(g)*invbetas(g)/((Ep(g)-Ep(g+1))*lightspeed)
    end do

    norm = sum(int_space*invbetas/lightspeed) ! This is for when ints(g) is given as AVG in group g. That's what I have, right?? Also, c cancels. That's fine I just don't want to forget it belongs.

    ints = ints/norm

    if (label .eq. "uncollided") then
        if (particle_type .eq. "photon") then
            printoutname = "LABEL: uncollided photon"
        else
            printoutname = "LABEL: uncollided electron"
        end if
    else
        if (particle_type .eq. "photon") then
            printoutname = "LABEL: photon"
        else
            printoutname = "LABEL: electron"
        end if
    end if

    print *, "------------------------------------------------------"
    print *, "### POST-PROCESSING: Energy bin endpoints printout."
    print *, printoutname
    do g = 1, l_G
        do n = 1, 2
            print *, Ep(g+n-1)
        end do
    end do
    print *, "### POST-PROCESSING: Energy bin distribution printout."
    print *, printoutname
    do g = 1, l_G
        do n = 1, 2
            print *, ints(g)
        end do
    end do
    print *, "------------------------------------------------------"
end subroutine MGXS_energy_spectrum

subroutine DFEM_MGXS_energy_spectrum(particle_type, label, vol, Ep, phi)
    implicit none
    character(*), intent(in) :: particle_type
    character(*), intent(in) :: label
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: Ep
    real, dimension(:,:,:), intent(in) :: phi

    integer :: e, g, k, n
    integer :: l_G
    real, dimension(:), allocatable :: int_space
    real, dimension(:), allocatable :: invbetas
    real, dimension(:), allocatable :: ints
    real :: norm
    character(100) :: printoutname

    l_G = size(Ep)-1 ! Use size of phi to inform whether FEXS or MGXS????

    allocate(int_space(l_G))
    allocate(invbetas(l_G))
    allocate(ints(l_G))

    int_space = 0
    do k = 1, NKe
        do e = 1, NE
            int_space = int_space + vol(e)*phi(e,k,1:l_G)/4
        end do
    end do

    if (particle_type .eq. "electron") then
        do g = 1, l_G
            invbetas(g) = e_mass_E*(sqrt((Ep(g)/e_mass_E)*(Ep(g)/e_mass_E+2)) - &
                sqrt((Ep(g+1)/e_mass_E)*(Ep(g+1)/e_mass_E+2)))/(Ep(g)-Ep(g+1))
        end do
    else
        invbetas = 1
    end if

    do g = 1, l_G
        ints(g) = int_space(g)*invbetas(g)/((Ep(g)-Ep(g+1))*lightspeed)
    end do

    norm = sum(int_space*invbetas/lightspeed)

    ints = ints/norm

    if (label .eq. "uncollided") then
        if (particle_type .eq. "photon") then
            printoutname = "LABEL: uncollided photon"
        else
            printoutname = "LABEL: uncollided electron"
        end if
    else
        if (particle_type .eq. "photon") then
            printoutname = "LABEL: photon"
        else
            printoutname = "LABEL: electron"
        end if
    end if

    print *, "------------------------------------------------------"
    print *, "### POST-PROCESSING: Energy bin endpoints printout."
    print *, printoutname
    do g = 1, l_G
        do n = 1, 2
            print *, Ep(g+n-1)
        end do
    end do
    print *, "### POST-PROCESSING: Energy bin distribution printout."
    print *, printoutname
    do g = 1, l_G
        do n = 1, 2
            print *, ints(g)
        end do
    end do
    print *, "------------------------------------------------------"
end subroutine DFEM_MGXS_energy_spectrum

subroutine FEXS_energy_spectrum(particle_type, label, Cekk, vol, Ep, phi)
    implicit none
    character(*), intent(in) :: particle_type
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: Ep
    real, dimension(:,:), intent(in) :: phi

    integer :: e, g, k
    integer :: l_G
    real, dimension(:), allocatable :: int_space
    real, dimension(:), allocatable :: invbetas
    real, dimension(:), allocatable :: ints
    character(100) :: printoutname
    real :: norm

    l_G = size(Ep)-1

    allocate(int_space(l_G+1))
    allocate(invbetas(l_G+1))
    allocate(ints(l_G+1))

    int_space = 0
    do k = 1, NKe
        do e = 1, NE
            int_space = int_space + vol(e)*phi(Cekk(e,k),:)/4
        end do
    end do

    if (particle_type .eq. "electron") then
        do g = 1, l_G+1
            invbetas(g) = 1/beta(Ep(g))
        end do
    else
        invbetas = 1
    end if

    do g = 1, l_G+1
        ints(g) = int_space(g)*invbetas(g)/lightspeed
    end do

    norm = 0
    do g = 1, l_G
        norm = norm + 0.5*(Ep(g)-Ep(g+1))*(ints(g)+ints(g+1))
    end do

    ints = ints/norm

    if (label .eq. "uncollided") then
        if (particle_type .eq. "photon") then
            printoutname = "LABEL: uncollided photon"
        else
            printoutname = "LABEL: uncollided electron"
        end if
    else
        if (particle_type .eq. "photon") then
            printoutname = "LABEL: photon"
        else
            printoutname = "LABEL: electron"
        end if
    end if

    print *, "-------------------------------------------------------"
    print *, "### POST-PROCESSING: Energy node printout."
    print *, printoutname
    do g = 1, l_G+1
        print *, Ep(g)
    end do
    print *, "### POST-PROCESSING: Energy node distribution printout."
    print *, printoutname
    do g = 1, l_G+1
        print *, ints(g)
    end do
    print *, "-------------------------------------------------------"
end subroutine FEXS_energy_spectrum

subroutine MGXS_calculate_deposition(Cekk, DEPXS, phi, DEP)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:), intent(in) :: DEPXS
    real, dimension(:,:), intent(in) :: phi

    real, dimension(:), allocatable, intent(inout) :: DEP

    integer :: e, g, k, mat

    allocate(DEP(NK))

    DEP = 0
    if (Nmats .eq. 1) then
        do k = 1, NK
            do g = 1, Ge
                DEP(k) = DEP(k) + phi(k,g)*DEPXS(g,1)
            end do
        end do
    else
        do k = 1, NKe
            do e = 1, NE
                do g = 1, Ge
                    DEP(Cekk(e,k)) = DEP(Cekk(e,k)) + phi(Cekk(e,k),g)*DEPXS(g,eltomat(e))/Nglobal(Cekk(e,k))
                end do
            end do
        end do
    end if
end subroutine MGXS_calculate_deposition

subroutine DFEM_MGXS_calculate_deposition(Cekk, DEPXS, phi, DEP)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:), intent(in) :: DEPXS
    real, dimension(:,:,:), intent(in) :: phi

    real, dimension(:,:), allocatable, intent(inout) :: DEP

    integer :: e, g, k

    allocate(DEP(NE,NKe))

    DEP = 0
    if (Nmats .eq. 1) then
        do k = 1, NKe
            DEP(1:NE,k) = DEP(1:NE,k) + matmul(phi(1:NE,k,:),DEPXS(:,1))
        end do
    else
        !$OMP parallel do collapse(2) shared(DEP) private(e, k)
        do k = 1, NKe
            do e = 1, NE
                DEP(e,k) = DEP(e,k) + sum(phi(e,k,:)*DEPXS(:,eltomat(e)))
            end do
        end do
        !$OMP end parallel do
    end if
end subroutine DFEM_MGXS_calculate_deposition

subroutine total_deposition(label, Cekk, vol, DEP)
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: DEP

    integer :: e, k
    real :: total

    total = 0

    do k = 1, NKe
        do e = 1, NE
            total = total + vol(e)*DEP(Cekk(e,k))/NKe
        end do
    end do

    if (label .eq. "energy") then
        print *, "---------------------------------------------"
        print *, "### POST-PROCESSING: Total energy deposition."
        print *, total
        print *, "---------------------------------------------"
    else if (label .eq. "dose") then
        print *, "-------------------------------------------"
        print *, "### POST-PROCESSING: Total dose deposition."
        print *, total
        print *, "-------------------------------------------"
    else if (label .eq. "charge") then
        print *, "---------------------------------------------"
        print *, "### POST-PROCESSING: Total charge deposition."
        print *, total
        print *, "---------------------------------------------"
    end if
end subroutine total_deposition

subroutine DFEM_total_deposition(label, vol, DEP)
    implicit none
    character(*), intent(in) :: label
    real, dimension(:), intent(in) :: vol
    real, dimension(:,:), intent(in) :: DEP

    integer :: e, k
    real :: total

    total = sum(matmul(transpose(DEP),vol))/NKe ! May not be general. Works for 1D, and 3D tetrahedra.

    if (label .eq. "energy") then
        print *, "---------------------------------------------"
        print *, "### POST-PROCESSING: Total energy deposition."
        print *, total
        print *, "---------------------------------------------"
    else if (label .eq. "charge") then
        print *, "---------------------------------------------"
        print *, "### POST-PROCESSING: Total charge deposition."
        print *, total
        print *, "---------------------------------------------"
    else if (label .eq. "dose") then
        print *, "-------------------------------------------"
        print *, "### POST-PROCESSING: Total dose deposition."
        print *, total
        print *, "-------------------------------------------"
    end if
end subroutine DFEM_total_deposition

subroutine save_deposition(label, DEP)
    implicit none
    character(*), intent(in) :: label
    real, dimension(:), intent(in) :: DEP

    character(100) :: fname

    fname = trim(adjustl(output_fname))//"/"//label//"_deposition.dat"

    open(1, file = fname, form = "unformatted")
    write(1) DEP
    close(1)

    if (label .eq. "energy") then
        print *, "---------------------------------------"
        print *, "### POST-PROCESSING: Energy deposition."
        print *, "Saved as:", fname
        print *, "---------------------------------------"
    else if (label .eq. "charge") then
        print *, "---------------------------------------"
        print *, "### POST-PROCESSING: Charge deposition."
        print *, "Saved as:", fname
        print *, "---------------------------------------"
    else if (label .eq. "dose") then
        print *, "-------------------------------------"
        print *, "### POST-PROCESSING: Dose deposition."
        print *, "Saved as:", fname
        print *, "-------------------------------------"
    end if
end subroutine save_deposition

subroutine visualize_mesh_slice_deposition(label, Cekk, slicek, DEP)
    implicit none
    character(*), intent(in) :: label
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:), intent(in) :: slicek
    real, dimension(:), intent(in) :: DEP

    integer :: k
    real, dimension(:), allocatable :: slicedep
    integer :: l_G
    character(80) :: fname

    allocate(slicedep(size(slicek)))

    do k = 1, size(slicek)
        slicedep(k) = DEP(slicek(k))
    end do

    if (results_as_dat) then
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_deposition.dat"

        open(1, file = fname, form = "unformatted")
        write(1) slicedep
        close(1)
    else
        fname = trim(adjustl(output_fname))//"/slice_"//label//"_deposition.txt"

        open(1, file = fname, status = "replace")
        do k = 1, size(slicek)
            write(1,*) slicedep(k)
        end do
        close(1)
    end if

    print *, "-------------------------------------------"
    print *, "### POST-PROCESSING: Mesh slice deposition."
    print *, "Saved as:", fname
    print *, "-------------------------------------------"
end subroutine visualize_mesh_slice_deposition

subroutine wiscoslab_DFEM_SN_MGXS_check_energy_conservation &
    (vol, mu, wL, Emid, EDEP, fEg, phiu, phi, psi)
    implicit none
    real, dimension(:), intent(in) :: vol
    real, dimension(:), intent(in) :: mu
    real, dimension(:), intent(in) :: wL
    real, dimension(:), intent(in) :: Emid
    real, dimension(:,:), intent(in) :: EDEP
    real, dimension(:), intent(in) :: fEg
    real, dimension(:,:,:), intent(in) :: phiu
    real, dimension(:,:,:), intent(in) :: phi
    real, dimension(:,:,:,:), intent(in) :: psi

    integer :: e, g, gpr, i, k
    integer :: l_G
    real :: S1
    real :: S2
    real :: S3
    real :: RHS
    real :: LHS

    l_G = size(Emid)

    !!! NOTE: These cross sections are delta-down'd. I think I should be using
    ! non delta-down'd XS, or I should include the delta-down terms. The latter is better.

    RHS = 0
    LHS = 0
    S1 = 0
    S2 = 0
    S3 = 0
    do g = 1, l_G
        ! Leakage
        do i = 1, Nmu
            S1 = S1 - twopi*Emid(g)*mu(i)*wL(i)*psi(NE,2,i,g) ! Just over outgoing??? Maybe not??? Backscatter should be negligible anyway
        end do
        if (FCS) then
            S1 = S1 - (-1)*phiu(NE,2,g)
        end if

        ! Energy deposited
        do k = 1, NKe
            do e = 1, NE
                S2 = S2 + &
                    EDEP(g,eltomat(e))*phi(e,k,g)*vol(e)/NKe
            end do
        end do

        ! Total energy
        if (beam_energy_dist .eq. "monochromatic") then
            S3 = S3 + merge(Emid(1), 0.0, g .eq. 1)
        else if (beam_energy_dist .eq. "polychromatic") then
            S3 = S3 + fEg(g)*Emid(g)
        else
            print *, "BOXCAR ENERGY CONS. NOT DONE"
            stop
        end if
    end do

    LHS = S1 + S2
    RHS = S3

    print *, "-----------------------------------------------"
    print *, "### POST-PROCESSING: Energy conservation check."
    print *, S1, S2, S3
    print *, LHS, RHS
    print *, "-----------------------------------------------"
end subroutine wiscoslab_DFEM_SN_MGXS_check_energy_conservation

subroutine FEXS_calculate_dose(Cekk, Ee, phi, Scol, DOSE)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:), intent(in) :: Ee
    real, dimension(:,:), intent(in) :: phi
    real, dimension(:,:), intent(in) :: Scol

    real, dimension(:), allocatable, intent(inout) :: DOSE

    integer :: e, g, k, mat

    allocate(DOSE(NK))

    print *, "FEXS DOSE NOT READY"
    stop
end subroutine FEXS_calculate_dose

subroutine DFEM_FEXS_calculate_deposition(Cekk, DEPXS, phi, DEP)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:,:), intent(in) :: DEPXS
    real, dimension(:,:,:,:), intent(in) :: phi

    real, dimension(:,:), allocatable, intent(inout) :: DEP

    integer :: e, g, k, n

    allocate(DEP(NE,NKe))

    DEP = 0
    if (Nmats .eq. 1) then
        do k = 1, NKe
            do n = 1, 2
                DEP(1:NE,k) = DEP(1:NE,k) + &
                    matmul(phi(1:NE,k,n,:),DEPXS(n,:,1))
            end do
        end do
    else
        !$OMP parallel do collapse(2) shared(DEP) private(e, k)
        do k = 1, NKe
            do e = 1, NE
                DEP(e,k) = DEP(e,k) + sum(phi(e,k,:,:)*DEPXS(:,:,eltomat(e)))
            end do
        end do
        !$OMP end parallel do
    end if
end subroutine DFEM_FEXS_calculate_deposition

subroutine planar_integrals(r, eprime, vol, a, bvector, tol, d0, dN, N, F)
    ! For arbitrary axes, take a rotation matrix. Or rotate mesh
    ! Currently takes coordinates w.r.t. mesh coord sys. NOT beam axis. Should be beam axis.
    implicit none
    real, dimension(:,:,:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: eprime
    real, dimension(:), intent(in) :: vol
    real, dimension(:,:), intent(in) :: a
    real, dimension(:,:,:), intent(in) :: bvector
    real, intent(in) :: tol
    real, intent(in) :: d0
    real, intent(in) :: dN
    integer, intent(in) :: N
    real, dimension(:,:), intent(in) :: F ! (e,k)

    integer :: dir, e, ep, face, i, v, k, kp, s
    real, dimension(:), allocatable :: Fbar ! (i) i.e., depth nodes. Not sure why I'm calling it Fbar, it's not an average over area. Should I make it? I would need cross sectional area. So maybe no.
    real :: z0
    logical, dimension(:), allocatable :: logic_elslist
    real, dimension(:), allocatable :: maxzs
    real, dimension(:), allocatable :: minzs
    integer, dimension(:), allocatable :: elslist
    logical, dimension(4) :: numline
    logical, dimension(4) :: hit
    integer, dimension(3) :: kpvals
    integer, dimension(4) :: tfvals
    real :: t
    real, dimension(:,:), allocatable :: Rvecs ! (x or y dir, success)
    real :: xbar
    real :: ybar
    real :: area
    real, dimension(3) :: Heron
    real :: HS
    real, dimension(4) :: angle
    real, dimension(4) :: x
    real, dimension(4) :: y
    integer, dimension(4) :: svals
    integer :: ind1
    integer :: ind2
    real :: contrib

    ! Wrote a lot of this down in my notebook. Write it up?

    allocate(Fbar(N))
    allocate(logic_elslist(NE))
    allocate(maxzs(NE))
    allocate(minzs(NE))
    allocate(Rvecs(2,4))

    do e = 1, NE
        maxzs(e) = maxval(r(3,1:NKe,e))
        minzs(e) = minval(r(3,1:NKe,e))
    end do

    Fbar = 0
    !!$OMP parallel do shared(Fbar) private(dir, e, ep, f, i, k, kp, v, s)
    do i = 1, N
        ! Form the z0 value at which the integral is being evaluated
        z0 = d0 + (i-1)*(dN-d0)/(N-1)

        ! Determine which elements contribute to this z0 slice
        logic_elslist = (z0 .gt. minzs - (dN-d0)*tol) .and. (z0 .lt. maxzs + (dN-d0)*tol)

        allocate(elslist, source = pack([(e,e=1,NE)], mask = logic_elslist))

        ! Determine the intersection of the plane at z0 with the tetrahedra
        do ep = 1, size(elslist)
            e = elslist(ep)

            numline = .false.
            hit = .false.
            do k = 1, NKe
                numline(k) = z0 .gt. r(3,k,e) ! Labels true all nodes above z0
                hit(k) = abs(z0 - r(3,k,e)) .lt. (dN-d0)*tol
            end do

            if (count(hit) .eq. 3) then
                ! Face is flat. Can immediately construct needed quantities, and become wary of
                ! double counting
                face = findloc(hit, .false., 1)
                xbar = 0
                ybar = 0

                kpvals = pack([(k,k=1,NKe)], mask = hit)

                do s = 1, 3
                    k = kpvals(s)

                    Rvecs(:,s) = r(1:2,k,e)
                end do

                ! Calculate the average value of x and y in the intersection face
                xbar = sum(Rvecs(1,1:3))/3
                ybar = sum(Rvecs(2,1:3))/3

                ! Calculate the area of the intersection face using Heron's formula
                ! Also, apply factor of 1/2 if this face will be repeated in elslist
                Heron(1) = norm2(Rvecs(:,1) - Rvecs(:,2))
                Heron(2) = norm2(Rvecs(:,1) - Rvecs(:,3))
                Heron(3) = norm2(Rvecs(:,2) - Rvecs(:,3))
                HS = 0.5*sum(Heron)

                area = merge(0.5, 1.0, any(elslist .eq. eprime(face,e)))*&
                    sqrt(HS*(HS-Heron(1))*(HS-Heron(2))*(HS-Heron(3)))
            else if (count(numline) .eq. 3 .or. count(numline) .eq. 1) then
                ! Intersection face is a triangle.
                ! The node whose logical differs from the rest is attached to every edge intersected.

                if (count(numline) .eq. 1) numline = .not. numline

                k = findloc(numline, .false., 1)

                kpvals = pack([(k,k=1,NKe)], mask = numline)

                do s = 1, 3
                    kp = kpvals(s)

                    ! Determine the intersection points by parametrizing the edge
                    t = (z0-r(3,kp,e))/(r(3,k,e)-r(3,kp,e))

                    ! Determine the intersection points themselves
                    do dir = 1, 2
                        Rvecs(dir,s) = t*(r(dir,k,e) - r(dir,kp,e)) + r(dir,kp,e)
                    end do
                end do

                ! Calculate the average value of x and y in the intersection face
                xbar = sum(Rvecs(1,1:3))/3
                ybar = sum(Rvecs(2,1:3))/3

                ! Calculate the area of the intersection face using Heron's formula
                Heron(1) = norm2(Rvecs(:,1) - Rvecs(:,2))
                Heron(2) = norm2(Rvecs(:,1) - Rvecs(:,3))
                Heron(3) = norm2(Rvecs(:,2) - Rvecs(:,3))
                HS = 0.5*sum(Heron)

                area = sqrt(HS*(HS-Heron(1))*(HS-Heron(2))*(HS-Heron(3)))
            else if (count(numline) .eq. 2) then
                ! Intersection face is a trapezoid.
                ! The four edges that connect one true and one false node are the intersected
                ! edges.

                tfvals(1:2) = pack([(k,k=1,NKe)], mask = numline) ! indices of true
                tfvals(3:4) = pack([(k,k=1,NKe)], mask = .not. numline) ! indices of false

                do s = 1, 4
                    k = tfvals(((-1)**s + 3)/2)
                    kp = tfvals(ceiling(0.5*s) + 2)

                    ! Determine the intersection points by parametrizing the edge
                    t = (z0-r(3,kp,e))/(r(3,k,e)-r(3,kp,e))

                    ! Determine the intersection points themselves
                    do dir = 1, 2
                        Rvecs(dir,s) = t*(r(dir,k,e) - r(dir,kp,e)) + r(dir,kp,e)
                    end do
                end do

                ! Calculate the average value of x and y in the intersection face
                xbar = sum(Rvecs(1,1:4))/4
                ybar = sum(Rvecs(2,1:4))/4

                ! Calculate the area of the intersection face by breaking trapezoid into two triangles.
                ! First have to organize the nodes clockwise or counterclockwise.
                ! A bit convoluted but this works.

                do s = 1, 4
                    angle(s) = atan2(Rvecs(2,s)-ybar, Rvecs(1,s)-xbar)
                end do

                x(1) = Rvecs(1,minloc(angle,1))
                y(1) = Rvecs(2,minloc(angle,1))
                x(4) = Rvecs(1,maxloc(angle,1))
                y(4) = Rvecs(2,maxloc(angle,1))

                svals = 0
                do s = 1, 4
                    if (s .eq. minloc(angle,1) .or. s .eq. maxloc(angle,1)) cycle
                    svals(s) = s
                end do

                ind1 = 5
                ind2 = 0
                do s = 1, 4
                    if (svals(s) .eq. 0) cycle
                    ind1 = min(ind1,svals(s))
                    ind2 = max(ind2,svals(s))
                end do

                ind1 = merge(ind1, ind2, angle(ind1) .lt. angle(ind2))
                ind2 = 10 - minloc(angle,1) - maxloc(angle,1) - ind1

                x(2) = Rvecs(1,ind1)
                y(2) = Rvecs(2,ind1)
                x(3) = Rvecs(1,ind2)
                y(3) = Rvecs(2,ind2)

                ! Now do Heron's formula on the two triangles formed by connecting vertex 2 and 4 (so, 1-2-4 and 2-3-4)
                Rvecs(:,1) = [x(1), y(1)]
                Rvecs(:,2) = [x(2), y(2)]
                Rvecs(:,3) = [x(4), y(4)]

                Heron(1) = norm2(Rvecs(:,1) - Rvecs(:,2))
                Heron(2) = norm2(Rvecs(:,1) - Rvecs(:,3))
                Heron(3) = norm2(Rvecs(:,2) - Rvecs(:,3))
                HS = 0.5*sum(Heron)

                area = sqrt(HS*(HS-Heron(1))*(HS-Heron(2))*(HS-Heron(3)))

                Rvecs(:,1) = [x(2), y(2)]
                Rvecs(:,2) = [x(3), y(3)]
                Rvecs(:,3) = [x(4), y(4)]

                Heron(1) = norm2(Rvecs(:,1) - Rvecs(:,2))
                Heron(2) = norm2(Rvecs(:,1) - Rvecs(:,3))
                Heron(3) = norm2(Rvecs(:,2) - Rvecs(:,3))
                HS = 0.5*sum(Heron)

                area = area + sqrt(HS*(HS-Heron(1))*(HS-Heron(2))*(HS-Heron(3)))
            else
                ! count is 4 or 0.
                ! Element shouldn't be here. Could be here due to floating point error,
                ! but in any case, it will contribute very small area if any.
                cycle
            end if

            ! Add up contributions from each node interpolated to this face
            contrib = 0
            do k = 1, NKe
                contrib = contrib + &
                    F(e,k)*(a(k,e) + bvector(1,k,e)*xbar + bvector(2,k,e)*ybar + bvector(3,k,e)*z0)
            end do
            contrib = contrib*area/(6*vol(e))

            Fbar(i) = Fbar(i) + contrib
        end do
        deallocate(elslist)
    end do

    print *, "-----------------------------------------------------------------"
    print *, "### POST-PROCESSING: Planar integrals along depth, depth values."
    do i = 1, N
        print *, d0 + (i-1)*(dN-d0)/(N-1)
    end do
    print *, "### POST-PROCESSING: Planar integrals along depth, integrals."
    do i = 1, N
        print *, Fbar(i)
    end do
    print *, "-----------------------------------------------------------------"
end subroutine planar_integrals

end module post_processing