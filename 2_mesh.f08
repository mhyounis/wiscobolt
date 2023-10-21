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

module mesh
    use OMP_LIB
    use math
    use physics
    use user_input
implicit none

! WIP LIST:
! Fix ray tracing for spherical open beams (what's wrong??)
! Add capability for faces parallel to ray

contains
subroutine old_import_nodes(rglobal)
    implicit none
    real, dimension(:,:), allocatable, intent(inout) :: rglobal

    allocate(rglobal(3,NK))

    open(1, file = "Mesh files/"//selected_mesh//"/use_nodes.dat")
    read(1,*) rglobal
    close(1)
end subroutine old_import_nodes

subroutine old_import_elements(Cekk)
    implicit none
    integer, dimension(:,:), allocatable, intent(inout) :: Cekk

    allocate(Cekk(NKe,NE))

    open(1, file = "Mesh files/"//selected_mesh//"/use_elements.dat")
    read(1,*) Cekk
    close(1)

    Cekk = transpose(Cekk)
end subroutine old_import_elements

subroutine old_mesh_analysis(Cekk, Cfkf, Ceff, eprime, fprime, intfc, bdyfc, bdyel)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk

    integer, dimension(:,:), allocatable, intent(inout) :: Cfkf
    integer, dimension(:,:), allocatable, intent(inout) :: Ceff
    integer, dimension(:,:), allocatable, intent(inout) :: eprime
    integer, dimension(:,:), allocatable, intent(inout) :: fprime
    integer, dimension(:,:,:), allocatable, intent(inout) :: intfc
    integer, dimension(:,:), allocatable, intent(inout) :: bdyfc
    integer, dimension(:), allocatable, intent(inout) :: bdyel

    integer :: e, f, fg
    integer :: NF
    integer :: NBF
    integer, dimension(:,:), allocatable :: B_Cfkf
    integer, dimension(:,:,:), allocatable :: B_intfc
    integer, dimension(:,:,:), allocatable :: C_intfc
    integer, dimension(:,:), allocatable :: B_bdyfc
    integer, dimension(:), allocatable :: B_bdyel
    integer, dimension(:,:), allocatable :: cyclematrix
    integer :: maxind
    integer :: minind
    integer :: fgval
    integer :: nextfg
    integer, dimension(3) :: temprow
    integer(kind = 8), dimension(:), allocatable :: tagarray
    integer(kind = 8) :: tag

    allocate(Ceff(NE,NFe))
    allocate(eprime(NFe,NE))
    allocate(fprime(NFe,NE))
    allocate(B_Cfkf(NKf,NE*NFe))
    allocate(B_intfc(2,NE*NFe,2))
    allocate(cyclematrix(NKf,NFe))
    allocate(tagarray(NE*NFe)) ! TRY not allocating it to this large (~2/14/23)

    cyclematrix(:,1) = [2,3,4]
    cyclematrix(:,2) = [1,3,4]
    cyclematrix(:,3) = [1,2,4]
    cyclematrix(:,4) = [1,2,3]

    B_Cfkf = 0
    eprime = 0
    fprime = 0
    B_intfc = 0
    tagarray = 0

    nextfg = 1
    do e = 1, NE
        do f = 1, NFe
            temprow = Cekk(e,cyclematrix(1:3,f))
            maxind = maxloc(temprow, dim = 1)
            minind = minloc(temprow, dim = 1)
            temprow = [temprow(maxind), temprow(6-maxind-minind), temprow(minind)]
            !call qsort_int(temprow)
            tag = temprow(1) + int(NK,8)*(temprow(2)-1) + int(NK,8)**2*(temprow(3)-1)
            if (any(tagarray(1:nextfg) .eq. tag)) then
                fgval = findloc(tagarray, tag, dim = 1)
                B_intfc(:,fgval,2) = [e,f]
                eprime(f,e) = B_intfc(1,fgval,1)
                fprime(f,e) = B_intfc(2,fgval,1)
                eprime(fprime(f,e),eprime(f,e)) = e
                fprime(fprime(f,e),eprime(f,e)) = f
                Ceff(e,f) = fgval
            else
                B_Cfkf(:,nextfg) = temprow
                B_intfc(:,nextfg,1) = [e,f]
                Ceff(e,f) = nextfg
                tagarray(nextfg) = tag
                nextfg = nextfg + 1
            end if
        end do
    end do
    deallocate(tagarray)

    NF = nextfg - 1

    allocate(Cfkf, source = transpose(B_Cfkf(:,1:NF)))
    deallocate(B_Cfkf)

    allocate(C_intfc, source = B_intfc(:,1:NF,:))
    deallocate(B_intfc)

    allocate(intfc, source = reshape(C_intfc, [size(C_intfc,2),size(C_intfc,1),size(C_intfc,3)], order = [2,1,3])) ! (fg, 1 for el or 2 for face, 1 for primary or 2 for secondary face)
    deallocate(C_intfc)

    allocate(B_bdyfc(NF,3))
    nextfg = 1
    do fg = 1, NF
        if (all(intfc(fg,:,2) .eq. 0)) then
            B_bdyfc(nextfg,1:2) = intfc(fg,:,1)
            B_bdyfc(nextfg,3) = fg
            nextfg = nextfg + 1
        end if
    end do

    NBF = nextfg - 1

    allocate(bdyfc, source = B_bdyfc(1:NBF,:))
    deallocate(B_bdyfc)

    allocate(bdyel(0))
    do e = 1, NE
        if (any(eprime(:,e) .eq. 0)) then
            allocate(B_bdyel(size(bdyel)+1))
            B_bdyel(1:size(bdyel)) = bdyel
            B_bdyel(size(B_bdyel)) = e
            call move_alloc(B_bdyel, bdyel)
        end if
    end do
end subroutine old_mesh_analysis

subroutine GMSH_mesh_import_and_analysis &
(rglobal, Cekk, Cfkf, Ceff, eprime, fprime, intfc, bdyfc, bdyel)
    implicit none
    real, dimension(:,:), allocatable, intent(inout) :: rglobal
    integer, dimension(:,:), allocatable, intent(inout) :: Cekk
    integer, dimension(:,:), allocatable, intent(inout) :: Cfkf
    integer, dimension(:,:), allocatable, intent(inout) :: Ceff
    integer, dimension(:,:), allocatable, intent(inout) :: eprime
    integer, dimension(:,:), allocatable, intent(inout) :: fprime
    integer, dimension(:,:,:), allocatable, intent(inout) :: intfc
    integer, dimension(:,:), allocatable, intent(inout) :: bdyfc
    integer, dimension(:), allocatable, intent(inout) :: bdyel

    integer :: e, f, fg, i, j, k, mat
    integer :: fgval
    integer :: nextfg
    real :: trash
    integer :: trashint
    integer :: b_trashint
    integer :: c_trashint
    integer, dimension(4) :: trashint4
    integer, dimension(4) :: b_trashint4
    integer, dimension(3) :: trashint3
    integer, dimension(100) :: long_trash
    integer, dimension(:), allocatable :: kept_trash
    character(60) :: trashchar
    integer, dimension(:), allocatable :: voltomat
    integer, dimension(:), allocatable :: voltobeam
    integer :: numvols
    character(10) :: teststr
    integer(kind = 8) :: tag
    integer(kind = 8), dimension(:), allocatable :: tagarray
    integer, dimension(:,:), allocatable :: B_Cfkf
    integer, dimension(:,:,:), allocatable :: B_intfc
    integer, dimension(:,:,:), allocatable :: C_intfc
    integer, dimension(:,:), allocatable :: B_bdyfc
    integer, dimension(:), allocatable :: B_bdyel
    integer, dimension(:,:), allocatable :: cyclemat
    logical, dimension(:), allocatable :: logic_nodesinbeam
    logical, dimension(:), allocatable :: logic_elsinbeam
    logical, dimension(:,:), allocatable :: logic_nodesinmat
    integer :: NF
    integer :: NBF
    real :: ompstart, ompfinish

    ! MUST CHECK WHAT CHANGES MUST BE MADE WHEN USING FOR EXAMPLE SPHERES, CYLINDERS

    allocate(kept_trash(10))
    kept_trash = 0

    ! NOTE: Can I extract boundary info directly? Like number of bounding surfaces, etc.????
    ! Thoroughly look for where I can simplify construction of eprime and all of that

    ompstart = omp_get_wtime()
    open(1, file = "Mesh files/"//trim(adjustl(selected_mesh))//"/mesh.msh", status = "old", action = "read")
    ! $MeshFormat
    read(1,'(A)') trashchar
    read(1,*) trash, trashint, trashint
    read(1,'(A)') trashchar
    ! $EndMeshFormat

    ! $Entities or PhysicalNames
    read(1,'(A)') trashchar
    if (trashchar .eq. "$PhysicalNames") then
        ! $PhysicalNames
        read(1,*) trashint
        numvols = trashint ! Just a placeholder
        do i = 1, numvols
            read(1,*) trashint, trashint, trashchar
        end do
        read(1,'(A)') trashchar
        ! $EndPhysicalNames

        read(1,'(A)') trashchar
        go to 1
    end if

    ! $Entities
1   continue
    read(1,*) trashint4 ! Number of [points, curves, surfaces, volumes]
    numvols = trashint4(4)
    allocate(voltomat(numvols))
    allocate(voltobeam(numvols))
    voltomat = 0
    voltobeam = 0
    do i = 1, 3 ! Reading points/curves/surfaces
        do j = 1, trashint4(i) ! Reading points (i = 1), curves (i = 2), or surfaces (i = 3)
            read(1,*) trashint ! Point/Curve/Surface tag
        end do
    end do

    do j = 1, numvols ! Reading volumes
        read(1,*) trashint, trash, trash, trash, trash, trash, trash, &
            b_trashint, kept_trash(1:b_trashint), c_trashint, long_trash(1:c_trashint)

        if (b_trashint > 1+Nmats) then
            print *, "STOP!"
            print *, "message wip"
            stop
        end if

        do k = 1, 1+Nmats ! MAX PHYSICAL GROUPS __MUST__ BE Nbeams+Nmats OR LESS...
            if (kept_trash(k) .eq. 0) cycle
            if (kept_trash(k) > 2000) then
                voltobeam(j) = kept_trash(k)
            else
                voltomat(j) = kept_trash(k)
            end if
        end do
        kept_trash = 0
    end do

    read(1,'(A)') trashchar
    ! $EndEntities

    if (all(voltomat .eq. 1) .or. all(voltomat .eq. 0)) then
        if (Nmats > 1) then
            print *, "STOP!"
            print *, "MODULE 2, SUBROUTINE GMSH_mesh_import_and_analysis:"
            print *, "It was found that your GMSH mesh volumes have no physical group, or"
            print *, "all have the same physical group, yet you are using more than one material in"
            print *, "your transport problem. You need to use physical groups to specify the materials"
            print *, "in your mesh by the volume to which they correspond. These physical groups"
            print *, "need to be labeled 1, ..., until the number of materials you specify in this solver,"
            print *, "so that the solver can assign the correct material to the set of elements found"
            print *, "in a volume."
            print *, "PROGRAM ENDING."
            stop
        end if
    end if

    ! $Nodes
    read(1,'(A)') trashchar
    if (trashchar .eq. "$PartitionedEntities") then
        print *, "STOP!"
        print *, "MODULE 2, SUBROUTINE GMSH_mesh_import_and_analysis:"
        print *, "$PartitionedEntities block was found in the GMSH mesh file."
        print *, "At the moment, PartitionedEntities are not allowed (I am not sure what they do"
        print *, "nor how they are created). Either provide a mesh without PartitionedEntities, or"
        print *, "perhaps just delete the PartitionedEntities block with a text editor."
        print *, "If you need to use PartitionedEntities, contact the author of this program (myounis@terpmail.umd.edu)."
        print *, "PROGRAM ENDING."
        stop
    end if
    k = 0
    read(1,*) trashint4
    NK = trashint4(2)

    allocate(rglobal(3,NK))
    do
        read(1,*) b_trashint4
        do j = 1, b_trashint4(4)
            read(1,*) trashint
        end do
        do j = 1, b_trashint4(4)
            k = k + 1
            read(1,*) rglobal(1:3,k)
        end do
        if (k .eq. NK) exit
    end do

    do
        read(1,'(A)') trashchar
        if (trashchar .eq. "$EndNodes") exit ! Why do I need to do this? Why do some mesh files have an extra line after the nodes and others don't???
    end do
    ! $EndNodes

    ompfinish = omp_get_wtime()
    !print *, "PART 1:", ompfinish - ompstart, "seconds"

    ompstart = omp_get_wtime()

    ! $Elements
    read(1,'(A)') trashchar
    read(1,*) trashint4
    NE = trashint4(4)
    allocate(Cekk(NE,NKe))

    allocate(tagarray(NE*NFe))
    allocate(cyclemat(NKf,NFe))
    cyclemat(:,1) = [2,3,4]
    cyclemat(:,2) = [1,3,4]
    cyclemat(:,3) = [1,2,4]
    cyclemat(:,4) = [1,2,3]
    allocate(Ceff(NE,NFe))
    allocate(eprime(NFe,NE))
    allocate(fprime(NFe,NE))
    allocate(B_Cfkf(NKf,NE*NFe))
    allocate(B_intfc(2,NE*NFe,2))
    allocate(eltomat(NE))
    if (Nmats .eq. 1) then
        eltomat = 1
    else
        eltomat = 0
    end if
    B_Cfkf = 0
    eprime = 0
    fprime = 0
    B_intfc = 0
    tagarray = 0
    if (beams_defined_in_mesh) then
        allocate(logic_nodesinbeam(NK))
        allocate(logic_elsinbeam(NE))
        logic_nodesinbeam = .false.
        logic_elsinbeam = .false.
    end if
    allocate(logic_nodesinmat(NK,Nmats))
    logic_nodesinmat = .false.
    elementloop: do
        read(1,*) trashint4
        if (trashint4(1) .eq. 3) then ! Elements can be read
            e = 0
            nextfg = 1
            do i = 1, numvols
                if (Nmats .gt. 1) eltomat(e+1:e+trashint4(4)) = voltomat(i) ! CHECK
                do j = 1, trashint4(4)
                    read(1,*) trashint, b_trashint4
                    e = e + 1
                    !print *, trashint, b_trashint4
                    Cekk(e,1:NKe) = b_trashint4

                    ! Mesh analysis
                    !call qsort_int(b_trashint4) ! Instead, read 1 by 1 and put in specific location?
                    do f = 1, NFe
                        trashint3 = b_trashint4(cyclemat(:,f))
                        call qsort_int(trashint3)
                        tag = trashint3(1) + int(NK,8)*(trashint3(2)-1) + int(NK,8)**2*(trashint3(3)-1)
                        !tag = b_trashint4(cyclemat(1,f)) + &
                        !int(NK,8)*(b_trashint4(cyclemat(2,f))-1) + &
                        !int(NK,8)**2*(b_trashint4(cyclemat(3,f))-1)
                        if (any(tagarray(1:nextfg) .eq. tag)) then
                            fgval = findloc(tagarray, tag, dim = 1)
                            B_intfc(:,fgval,2) = [e,f]
                            eprime(f,e) = B_intfc(1,fgval,1)
                            fprime(f,e) = B_intfc(2,fgval,1)
                            eprime(fprime(f,e),eprime(f,e)) = e
                            fprime(fprime(f,e),eprime(f,e)) = f
                            Ceff(e,f) = fgval
                        else
                            B_Cfkf(1:NKf,nextfg) = trashint3
                            B_intfc(:,nextfg,1) = [e,f]
                            Ceff(e,f) = nextfg
                            tagarray(nextfg) = tag
                            nextfg = nextfg + 1
                        end if

                        !tag(f,e) =
                        ! Store tags then work with them at the very end??
                    end do

                    !Cekk(e,1:NKe) = b_trashint4
                    if (voltobeam(i) > 0 .and. beams_defined_in_mesh) then ! at the moment I don't care WHICH beam has the nodes.
                        logic_nodesinbeam(b_trashint4) = .true.
                        logic_elsinbeam(e) = .true.
                    end if
                    if (voltomat(i) > 0) then
                        logic_nodesinmat(b_trashint4,voltomat(i)) = .true.
                    end if
                end do

                if (i < numvols) read(1,*) trashint4
            end do

            exit elementloop
        else
            do i = 1, trashint4(4) ! The user either specified curves, lines, or points in physical groups, or they didn't specify physical groups
                read(1,*) trashint
            end do
        end if

        if (trashint4(1) > 3) then
            print *, "STOP!"
            print *, "MODULE 2, SUBROUTINE GMSH_mesh_import_and_analysis:"
            print *, "I think your GMSH mesh file has 4D or higher elements. Not sure how you"
            print *, "managed to do this, or if it's even possible. If you are not using physical"
            print *, "groups, go ahead and create one single physical group in your mesh and put"
            print *, "in all of the volumes you have. This should ensure that your GMSH mesh file"
            print *, "will only contain 3D elements."
            print *, "PROGRAM ENDING."
            stop
        end if
    end do elementloop
    read(1,'(A)') trashchar
    ! $EndElements
    close(1)

    ompfinish = omp_get_wtime()
    !print *, "PART 2:", ompfinish - ompstart, "seconds"

    ompstart = omp_get_wtime()

    deallocate(tagarray)

    ! Do work on global face arrays
    NF = nextfg - 1

    allocate(Cfkf, source = transpose(B_Cfkf(:,1:NF)))
    deallocate(B_Cfkf)

    allocate(C_intfc, source = B_intfc(:,1:NF,:))
    deallocate(B_intfc)

    allocate(intfc, source = reshape(C_intfc, [size(C_intfc,2),size(C_intfc,1),size(C_intfc,3)], order = [2,1,3])) ! (fg, 1 for el or 2 for face, 1 for primary or 2 for secondary face)
    deallocate(C_intfc)

    allocate(B_bdyfc(NF,3))
    nextfg = 1
    do fg = 1, NF
        if (all(intfc(fg,:,2) .eq. 0)) then
            B_bdyfc(nextfg,1:2) = intfc(fg,:,1)
            B_bdyfc(nextfg,3) = fg
            nextfg = nextfg + 1
        end if
    end do

    NBF = nextfg - 1

    allocate(bdyfc, source = B_bdyfc(1:NBF,:))
    deallocate(B_bdyfc)

    ompfinish = omp_get_wtime()
    !print *, "PART 3:", ompfinish - ompstart, "seconds"

    ompstart = omp_get_wtime()

    ! Do work on nodesinbeam
    if (beams_defined_in_mesh) then
        allocate(nodesinbeam, source = pack([(k,k=1,NK)], mask = logic_nodesinbeam))
        allocate(elsinbeam, source = pack([(e,e=1,NE)], mask = logic_elsinbeam))
        if (size(nodesinbeam) .eq. 0) then
            print *, "STOP!"
            print *, "MODULE 2, SUBROUTINE GMSH_mesh_import_and_analysis:"
            print *, "You have specified that the beam has been defined in the mesh,"
            print *, "but there appear to be no nodes in the beam. This could mean"
            print *, "you shouldn't have chosen this option, or there is something"
            print *, "wrong with your mesh. Check your input file and your mesh."
            print *, "PROGRAM ENDING"
            stop
        end if
        if (size(elsinbeam) .eq. 0) then
            print *, "STOP!"
            print *, "MODULE 2, SUBROUTINE GMSH_mesh_import_and_analysis:"
            print *, "You have specified that the beam has been defined in the mesh,"
            print *, "but there appear to be no elements in the beam. This could mean"
            print *, "you shouldn't have chosen this option, or there is something"
            print *, "wrong with your mesh. Check your input file and your mesh."
            print *, "PROGRAM ENDING"
            stop
        end if
    end if

    ! Do work on nodesinmat
    if (Nmats .eq. 1) then
        deallocate(logic_nodesinmat)
    else
        allocate(Nelnodes(NK,Nmats))
        Nelnodes = 0
        do k = 1, NKe
            do e = 1, NE
                Nelnodes(Cekk(e,k),eltomat(e)) = Nelnodes(Cekk(e,k),eltomat(e)) + 1
            end do
        end do


        ! SHOULD I MAKE NELNODES OF SIZE (NKmat,Nmats)????

        allocate(NKmat(Nmats))
        do mat = 1, Nmats
            NKmat(mat) = count(logic_nodesinmat(1:NK,mat))
        end do

        allocate(nodesinmat(maxval(NKmat),Nmats))
        nodesinmat = 0
        do mat = 1, Nmats
            nodesinmat(1:NKmat(mat),mat) = pack([(k,k=1,NK)], mask = logic_nodesinmat(:,mat))
        end do
    end if

    ! Do work on elsinmat
    ! REORDER ELEMENTS BASED ON elsinmat????
    ! This way, when it comes time to iterate, I can just do like start:finish in elements.
    if (Nmats .gt. 1) then
        allocate(NEmat(Nmats))
        do mat = 1, Nmats
            NEmat(mat) = count(eltomat .eq. mat)
        end do

        allocate(elsinmat(maxval(NEmat),Nmats))
        elsinmat = 0
        do mat = 1, Nmats
            elsinmat(1:NEmat(mat),mat) = pack([(e,e=1,NE)], mask = eltomat .eq. mat)
        end do
    end if

    ompfinish = omp_get_wtime()
    !print *, "PART 4:", ompfinish - ompstart, "seconds"

    ! If I want to revisit storing the tags while reading so as not to slow it down !! What did I mean by this...
    !do e = 1, NE
    !    do f = 1, NFe
    !        if (tag(f,e) .eq. 0) cycle
    !        sharedface = findloc(tag(f,e), tag(1:NFe,e:NE))
    !        if (sharedface .eq. 0) then
    !        tag(f,e) = 0
    !    end do
    !end do

    allocate(bdyel(0))
    do e = 1, NE
        if (any(eprime(:,e) .eq. 0)) then
            allocate(B_bdyel(size(bdyel)+1))
            B_bdyel(1:size(bdyel)) = bdyel
            B_bdyel(size(B_bdyel)) = e
            call move_alloc(B_bdyel, bdyel)
        end if
    end do
end subroutine GMSH_mesh_import_and_analysis

subroutine construct_elemental_mesh(rglobal, Cekk, r)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    integer, dimension(:,:), intent(in) :: Cekk

    real, dimension(:,:,:), allocatable, intent(inout) :: r

    integer :: e, k

    allocate(r(3,NKe,NE))

    do e = 1, NE
        do k = 1, NKe
            r(:,k,e) = rglobal(:,Cekk(e,k))
        end do
    end do
end subroutine construct_elemental_mesh

subroutine determine_overlapping_nodes(Cekk)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk

    integer :: k
    integer :: mink

    allocate(Nglobal(NK))

    do k = 1, NK
        Nglobal(k) = count(Cekk .eq. k)
    end do

    if (any(Nglobal .eq. 0)) then
        print *, "STOP!"
        print *, "MODULE 2, SUBROUTINE determine_overlapping_nodes:"
        print *, "There are some nodes in your mesh which belong to no elements, or do not exist in"
        print *, "the mesh at all yet are still being indexed. I'm not sure why this happens sometimes."
        print *, "I believe this is a problem with GMSH. I recommend re-generating this mesh with"
        print *, "slightly more or less mesh density and re-running. The bad nodes are:"
        do k = 1, NK
            if (Nglobal(k) .eq. 0) print *, k
        end do
        print *, "PROGRAM ENDING"
        stop
    end if
end subroutine determine_overlapping_nodes

subroutine find_identical_nodes(Cekk, eprime, kprime)
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:,:), intent(in) :: eprime

    integer, dimension(:,:,:), allocatable, intent(inout) :: kprime

    integer :: e, f, k
    integer, dimension(NKe) :: tempA

    allocate(kprime(NKe,NFe,NE))

    do e = 1, NE ! Can be optimized.
        do f = 1, NFe
            tempA = Cekk(eprime(f,e),1:NKe)
            do k = 1, NKe
                kprime(k,f,e) = findloc(tempA, Cekk(e,k), dim = 1)
            end do
        end do
    end do

    !do e = 1, NE
    !    do f = 1, NFe
    !        do k = 1, NKe
    !            if (kprime(k,f,e) /= 0) then
    !                if (any(r(e,k,:) /= r(eprime(e,f),kprime(e,f,k),:))) print *, "PROBLEM"
    !            end if
    !        end do
    !    end do
    !end do
end subroutine find_identical_nodes

!subroutine construct_normal_vectors(r, normal)
!    implicit none
!    real, dimension(:,:,:), intent(in) :: r
!
!    real, dimension(:,:,:), allocatable, intent(inout) :: normal
!
!    integer :: e, f
!    integer, dimension(:,:), allocatable :: fmatrix
!
!    allocate(normal(3,NFe,NE))
!    allocate(fmatrix(NFe,3))
!
!    fmatrix(1,:) = [2,3,4]
!    fmatrix(2,:) = [1,3,4]
!    fmatrix(3,:) = [1,2,4]
!    fmatrix(4,:) = [1,2,3]
!
!    do f = 1, NFe
!        do e = 1, NE
!            normal(:,f,e) = cross_product(r(:,fmatrix(f,3),e)-r(:,fmatrix(f,1),e),r(:,fmatrix(f,2),e)-r(:,fmatrix(f,1),e))
!            normal(:,f,e) = normal(:,f,e)/length(normal(:,f,e))
!            if (dot_product(normal(:,f,e),r(:,fmatrix(f,1),e)-r(:,f,e)) < 0.) normal(:,f,e) = -normal(:,f,e)
!        end do
!    end do
!
!    !do e = 1, NE
!    !    do f = 1, NFe
!    !        if (eprime(e,f) /= 0) then
!    !            if (any(abs(normal(:,f,e) + normal(:,fprime(e,f),eprime(e,f))) > 1E-6)) then
!    !                print *, abs(normal(:,f,e) + normal(:,fprime(e,f),eprime(e,f)))
!    !            end if
!    !        end if
!    !    end do
!    !end do
!end subroutine construct_normal_vectors

subroutine construct_normal_vectors(rglobal, r, Cekk, Ceff, eprime, globalnormal, normal, sgn)
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:,:,:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: Cekk
    integer, dimension(:,:), intent(in) :: Ceff
    integer, dimension(:,:), intent(in) :: eprime

    real, dimension(:,:), allocatable, intent(inout) :: globalnormal ! (dir,fg)
    real, dimension(:,:,:), allocatable, intent(inout) :: normal ! (dir,f,e)
    logical, dimension(:,:), allocatable, intent(inout) :: sgn

    integer :: e, f
    integer, dimension(:,:), allocatable :: fmatrix

    allocate(globalnormal(3,maxval(Ceff)))
    allocate(normal(3,NFe,NE))
    allocate(sgn(NFe,NE))

    sgn = .false.

    !do fg = 1, size(Cfkf,1)
    !    globalnormal(1:3,fg) = cross_product(rglobal(1:3,Cfkf(fg,2))-rglobal(1:3,Cfkf(fg,1)),&
    !        rglobal(1:3,Cfkf(fg,3))-rglobal(1:3,Cfkf(fg,1)))
    !
    !    globalnormal(1:3,fg) = sign(1.,dot_product(globalnormal(1:3,fg),rglobal(1:3,Cekk(intfc(fg,1,1),intfc(fg,2,1)))))*& ! NORMALIZATION CORRESPONDS TO "PRIMARY FACES," AS PER intfc ARRAY
    !        globalnormal(1:3,fg)/length(globalnormal(1:3,fg))
    !    sgn(intfc(fg,2,1),intfc(fg,1,1)) = .true.
    !end do

    allocate(fmatrix(NFe,3))

    fmatrix(1,:) = [2,3,4]
    fmatrix(2,:) = [1,3,4]
    fmatrix(3,:) = [1,2,4]
    fmatrix(4,:) = [1,2,3]

    do f = 1, NFe
        do e = 1, NE
            normal(:,f,e) = cross_product(r(:,fmatrix(f,3),e)-r(:,fmatrix(f,1),e),r(:,fmatrix(f,2),e)-r(:,fmatrix(f,1),e))
            normal(:,f,e) = normal(:,f,e)/norm2(normal(:,f,e))
            if (dot_product(normal(:,f,e),r(:,fmatrix(f,1),e)-r(:,f,e)) .lt. 0) normal(:,f,e) = -normal(:,f,e)

            if (e > eprime(f,e)) then
                globalnormal(:,Ceff(e,f)) = normal(:,f,e)
                sgn(f,e) = .true.
            end if
        end do
    end do
end subroutine construct_normal_vectors

subroutine modify_Cfkf(rglobal, r, bdyfc, Cfkf)
    ! Probably inefficient. Figure out how to not need to do this (~5/28/23)
    ! Reorganizes Cfkf for boundary elements so that their normal vectors constructed with:
    ! [r(3)-r(1)]x[r(2)-r(1)]
    ! Are oriented correctly
    ! This means that, while ray tracing, normal vectors are pointing outwards
    ! so det < 0 means a face is kept.
    implicit none
    real, dimension(:,:), intent(in) :: rglobal
    real, dimension(:,:,:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: bdyfc

    integer, dimension(:,:), intent(inout) :: Cfkf

    integer :: f, k
    real, dimension(3) :: testnormal
    real :: testval
    integer :: tempstore

    do f = 1, size(bdyfc,1)
        testnormal = cross_product( &
            rglobal(:,Cfkf(bdyfc(f,3),3))-rglobal(:,Cfkf(bdyfc(f,3),1)),&
            rglobal(:,Cfkf(bdyfc(f,3),2))-rglobal(:,Cfkf(bdyfc(f,3),1)))
        testval = dot_product(testnormal,rglobal(:,Cfkf(bdyfc(f,3),1))-r(:,bdyfc(f,2),bdyfc(f,1)))
        if (testval .lt. 0) then
            tempstore = Cfkf(bdyfc(f,3),3)
            Cfkf(bdyfc(f,3),3) = Cfkf(bdyfc(f,3),2)
            Cfkf(bdyfc(f,3),2) = tempstore
        end if
    end do
end subroutine modify_Cfkf

subroutine scale_and_translate_mesh(rglobal, r) ! SWAP ORDER OF SCALE AND TRANSLATE ONCE MESH IS DEFINED AS CENTERED
    implicit none
    real, dimension(:,:), allocatable, intent(inout) :: rglobal
    real, dimension(:,:,:), allocatable, intent(inout) :: r

    rglobal(1,:) = scale_factor(1)*rglobal(1,:) + translate_vector(1)
    rglobal(2,:) = scale_factor(2)*rglobal(2,:) + translate_vector(2)
    rglobal(3,:) = scale_factor(3)*rglobal(3,:) + translate_vector(3)

    r(1,:,:) = scale_factor(1)*r(1,:,:) + translate_vector(1)
    r(2,:,:) = scale_factor(2)*r(2,:,:) + translate_vector(2)
    r(3,:,:) = scale_factor(3)*r(3,:,:) + translate_vector(3)
end subroutine scale_and_translate_mesh

subroutine center_of_face(r, rfcm)
    implicit none
    real, dimension(:,:,:), intent(in) :: r

    real, dimension(:,:,:), allocatable, intent(inout) :: rfcm

    integer :: e, f
    integer, dimension(:,:), allocatable :: cyclematrix

    allocate(rfcm(3,NFe,NE))
    allocate(cyclematrix(NFe,3))

    cyclematrix(1,:) = [2,3,4]
    cyclematrix(2,:) = [1,3,4]
    cyclematrix(3,:) = [1,2,4]
    cyclematrix(4,:) = [1,2,3]

    do f = 1, NFe
        do e = 1, NE
            rfcm(:,f,e) = (r(:,cyclematrix(f,1),e) + r(:,cyclematrix(f,2),e) + r(:,cyclematrix(f,3),e))/3.
        end do
    end do

    deallocate(cyclematrix)
end subroutine center_of_face

subroutine determine_sparsity_pattern(Cekk, Ia, Ja, NonZ)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk

    integer, dimension(:), allocatable, intent(inout) :: Ia
    integer, dimension(:), allocatable, intent(inout) :: Ja
    integer, intent(inout) :: NonZ

    integer :: e, i, k, kp
    integer, dimension(:,:), allocatable :: Nodecounter

    allocate(Ia(NK+1))
    allocate(Nodecounter(NK,NK))

    ! TO TURN AN ARRAY INTO ITS SPARSE COUNTER PART, USE PACK. (WITH TRANSPOSE? I THINK I FIGURED A WAY TO DO IT WITHOUT BUT I FORGOT)

    Nodecounter = 0
    do kp = 1, NKe
        do k = 1, NKe
            do e = 1, NE
                Nodecounter(Cekk(e,k),Cekk(e,kp)) = Nodecounter(Cekk(e,k),Cekk(e,kp)) + 1
            end do
        end do
    end do
    ! If a global node k is within an element with another global node kp, then A(k,kp) /= 0.
    ! USE WHERE?

    ! Ia(row) counts the number of nonzero els up to BUT NOT including the row "row," + 1.
    Ia(1) = 1
    do k = 2, NK+1
        Ia(k) = Ia(k-1) + count(Nodecounter(:,k-1) /= 0) ! Should typically be Nodecounter(k-1,:), but Nodecounter is symmetric
    end do

    NonZ = Ia(NK+1) - 1

    allocate(Ja(NonZ))

    ! Ja()
    do k = 1, NK
        Ja(Ia(k):Ia(k+1)-1) = pack([(i,i=1,NK)], mask = Nodecounter(:,k) /= 0)
    end do
end subroutine determine_sparsity_pattern

subroutine determine_nodes_in_beam(Cekk, rglobal)
    implicit none
    integer, dimension(:,:), intent(in) :: Cekk
    real, dimension(:,:), intent(in) :: rglobal

    integer :: e, k
    real, dimension(3) :: k0
    real, dimension(3) :: r0
    real :: alpha0
    real :: beta0
    real :: Lx
    real :: Ly
    real :: S
    logical, dimension(:), allocatable :: logicarray
    real, dimension(:), allocatable :: dotprods
    real, dimension(:), allocatable :: projs
    real :: rho0
    real, dimension(:), allocatable :: rho0s
    real :: l_norm
    real :: l_dotprod
    real :: phi
    real, dimension(:), allocatable :: deltas
    real :: apexcos
    real :: tol = 0
    integer, dimension(:), allocatable :: B_elsinbeam
    logical, dimension(:), allocatable :: nodecond

    k0 = beam_coord_sys(:,3)
    r0 = beam_origin
    allocate(logicarray(NK))
    logicarray = .false.

    rho0 = beam_cutout_params(2)
    Lx = beam_cutout_params(2)
    Ly = beam_cutout_params(3)
    S = beam_cutout_params(1)

    if (beam_angular_dist .eq. "spherical") then
        allocate(dotprods(NK))
        do k = 1, NK
            dotprods(k) = dot_product(k0,rglobal(:,k)-r0)/norm2(rglobal(:,k)-r0)
        end do

        if (beam_cutout .eq. "circle") then
            apexcos = 1/sqrt(1+(rho0/S)**2)
            ! Should test this out

            logicarray = dotprods .ge. apexcos
        else if (beam_cutout .eq. "rectangle") then
            allocate(deltas(NK))
            do k = 1, NK
                l_norm = norm2((rglobal(:,k)-r0) - dot_product(k0,(rglobal(:,k)-r0))*k0)
                if (l_norm .eq. 0) then
                    deltas(k) = 0
                else
                    l_dotprod = dot_product(beam_coord_sys(:,1),(rglobal(:,k)-r0))
                    if (l_dotprod .ge. l_norm) then
                        phi = 0 ! Why did I do this again?
                    else
                        phi = acos(l_dotprod/l_norm)
                    end if
                    deltas(k) = cos(atan(1.0/&
                    (2*S*max(cos(phi)/Lx,sin(phi)/Ly,cos(phi-pi)/Lx,sin(phi-pi)/Ly))))
                end if
            end do

            logicarray = dotprods .ge. deltas
        end if
    else if (beam_angular_dist .eq. "planar") then
        allocate(projs(NK))
        do k = 1, NK
            projs(k) = norm2((rglobal(:,k)-r0) - dot_product(k0,(rglobal(:,k)-r0))*k0)
        end do

        if (beam_cutout .eq. "circle") then
            logicarray = projs .le. rho0
        else if (beam_cutout .eq. "rectangle") then
            Lx = beam_cutout_params(2)
            Ly = beam_cutout_params(3)
            allocate(rho0s(NK))
            do k = 1, NK
                if (projs(k) .eq. 0) then
                    rho0s(k) = 0
                else
                    l_dotprod = dot_product(beam_coord_sys(:,1),(rglobal(:,k)-r0))
                    if (l_dotprod .ge. l_norm) then
                        phi = 0
                    else
                        phi = acos(l_dotprod/projs(k))
                    end if
                    rho0s(k) = 1/(2*max(cos(phi)/Lx,sin(phi)/Ly,cos(phi-pi)/Lx,sin(phi-pi)/Ly))
                end if
            end do

            logicarray = projs .le. rho0s
        end if
    else if (beam_angular_dist .eq. "pencil") then
        allocate(dotprods(NK))
        allocate(rho0s(NK))

        do k = 1, NK
            dotprods(k) = dot_product(rglobal(:,k)-r0,k0)
            rho0s(k) = norm2(rglobal(:,k) - r0)
        end do

        logicarray = abs(dotprods - rho0s) .le. tol
        ! MUST APPEND WITH ERROR TOLERANCE
    end if

    allocate(nodesinbeam, source = pack([(k,k=1,NK)], mask = logicarray))

    allocate(nodecond(NKe))
    allocate(elsinbeam(0))
    do e = 1, NE
        do k = 1, NKe
            nodecond(k) = any(nodesinbeam .eq. Cekk(e,k))
        end do
        if (all(nodecond)) then
            allocate(B_elsinbeam(size(elsinbeam)+1))
            B_elsinbeam(1:size(elsinbeam)) = elsinbeam
            B_elsinbeam(size(elsinbeam)+1) = e
            call move_alloc(B_elsinbeam, elsinbeam)
        end if
    end do

    if (size(nodesinbeam) .eq. 0) then
        print *, "STOP!"
        print *, "MODULE 2, SUBROUTINE determine_nodes_in_beam:"
        print *, "Number of nodes in beam is found to be zero. Double check that the mesh"
        print *, "and input parameters are appropriate."
        print *, "PROGRAM ENDING."
        stop
    end if
end subroutine determine_nodes_in_beam

subroutine find_material_interfaces(intfc, bdyfc, matfc)
    implicit none
    integer, dimension(:,:,:), intent(in) :: intfc
    integer, dimension(:,:), intent(in) :: bdyfc

    integer, dimension(:,:), allocatable, intent(inout) :: matfc

    integer :: f, mat
    integer :: NF
    integer :: NBF
    integer, dimension(:,:), allocatable :: B_matfc
    integer :: e1
    integer :: e2
    integer :: mat1
    integer :: mat2
    integer, dimension(:), allocatable :: lengths
    integer :: newlength

    NF = size(intfc,1)
    NBF = size(bdyfc,1)
    allocate(matfc(0,Nmats))
    allocate(lengths(Nmats))

    lengths = 0

    do f = 1, NF
        e1 = intfc(f,1,1)
        e2 = intfc(f,1,2)

        if (e1 .eq. 0 .or. e2 .eq. 0) cycle ! Does intfc include bdy faces??? I think it does. If not, this is not needed

        if (eltomat(e1) .ne. eltomat(e2)) then ! SHOULD I MODIFY CFKF TO MAKE NORMALS RIGHT? NAH, NO NEED. JUST DON'T DO CULLING. ALSO ROOT OUT DET ~ 0.
            mat1 = eltomat(e1)
            mat2 = eltomat(e2)
            lengths(mat1) = lengths(mat1) + 1
            lengths(mat2) = lengths(mat2) + 1

            ! Size of matfc is given by the max number of faces owned by any material.
            ! If the materials to which f is being appended have fewer than this, the newlength is the same as the size of matfc
            ! If either of the materials to which f is being appended now has more than this (it has size(matfc,1)+1), then the newlength
            ! is the now siz(matfc,1)+1, which is going to be lengths(whatever material)
            newlength = max(size(matfc,1),lengths(mat1),lengths(mat2))

            allocate(B_matfc(newlength,Nmats))
            B_matfc(1:size(matfc,1),:) = matfc
            B_matfc(lengths(mat1),:) = f
            B_matfc(lengths(mat2),:) = f
            do mat = 1, Nmats
                ! The newly created row must all be zero for materials which did not gain a new face.
                if (mat .eq. mat1 .or. mat .eq. mat2) cycle
                B_matfc(newlength,mat) = 0
            end do
            call move_alloc(B_matfc, matfc)
        end if
    end do

    ! Now append bdyfc in matfc
    ! If intfc includes bdy faces, I should just include this in the previous loop...
    do f = 1, NBF
        mat1 = eltomat(bdyfc(f,1))

        lengths(mat1) = lengths(mat1) + 1
        newlength = max(size(matfc,1),lengths(mat1))

        allocate(B_matfc(newlength,Nmats))
        B_matfc(1:size(matfc,1),:) = matfc
        B_matfc(lengths(mat1),:) = f
        do mat = 1, Nmats
            ! The newly created row must all be zero for materials which did not gain a new face.
            if (mat .eq. mat1) cycle
            B_matfc(newlength,mat) = 0
        end do
        call move_alloc(B_matfc, matfc)
    end do
end subroutine find_material_interfaces

subroutine wiscoslab_all_transport_preproc(rglobal, Cekk, r, dz)
    ! Temporary? Could keep it like this...
    implicit none
    real, dimension(:,:), allocatable, intent(inout) :: rglobal
    integer, dimension(:,:), allocatable, intent(inout) :: Cekk
    real, dimension(:,:,:), allocatable, intent(inout) :: r
    real, dimension(:), allocatable, intent(inout) :: dz ! (actually delta z)

    integer :: e, k
    real :: L
    real :: tmat1
    real :: tmat2
    integer :: NE1
    integer :: NE2

    beam_angular_dist = "planar"
    beam_cutout = "none"

    if (k0(3) .ne. -1.0) then
        print *, "you're using the wrong beam axis"
        stop
    end if

    NE = 200
    NK = NE + 1
    NKe = 2
    NFe = 2

    allocate(rglobal(1,NK))
    allocate(Cekk(NE,NKe))
    allocate(r(1,NKe,NE))
    allocate(dz(NE))

    do k = 1, NKe
        do e = 1, NE
            Cekk(e,k) = e + k - 1
        end do
    end do

    ! Linearly spaced mesh. Homogeneous.
    !L = 12.7E-2
    L = 0.2
    dz = L/NE
    allocate(eltomat(NE))
    eltomat = 1

    !! Exponentially spaced mesh from z = L to z = 3*L/4 with 3*NE/4 elements there. Homogeneous.
    !L = 10.0
    !dz = (3*L/4)/(NE/4)
    !do e = 1, 3*NE/4
    !    dz(e) = log((exp(3*L/4)-exp(L))*(e-1)/(3*NE/4) + exp(L)) - &
    !        log((exp(3*L/4)-exp(L))*e/(3*NE/4) + exp(L))
    !end do
    !allocate(eltomat(NE))
    !eltomat = 1

    !! Linearly spaced within the slabs. Heterogeneous.
    !L = 0.18
    !tmat1 = 0.015917778
    !tmat2 = 0.010611852
    !NE1 = 50
    !NE2 = 100
    !
    !dz(1:NE1) = tmat1/NE1
    !dz(NE1+1:NE1+NE2) = tmat2/NE2
    !dz(NE1+NE2+1:NE) = (L-tmat1-tmat2)/(NE-NE1-NE2)
    !
    !allocate(eltomat(NE))
    !eltomat(1:NE1) = 1
    !eltomat(NE1+1:NE1+NE2) = 2
    !eltomat(NE1+NE2+1:NE) = 1

    rglobal(1,1) = L
    do k = 2, NK
        rglobal(1,k) = rglobal(1,k-1) - dz(k-1)
    end do
    do k = 1, NKe
        do e = 1, NE
            r(1,k,e) = rglobal(1,Cekk(e,k))
        end do
    end do
end subroutine wiscoslab_all_transport_preproc

end module mesh