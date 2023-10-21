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

module sweep_order
    use math
    use physics
    use user_input
implicit none

contains
!subroutine construct_sweeplist(bdyfc, eprime, normal, khat, sweeplist, sweepbounds, nosweeps)
!    implicit none
!    integer, dimension(:,:), intent(in) :: bdyfc
!    integer, dimension(:,:), intent(in) :: eprime
!    real, dimension(:,:,:), intent(in) :: normal
!    real, dimension(:,:,:), intent(in) :: khat
!
!    integer, dimension(:,:), allocatable, intent(inout) :: sweeplist
!    integer, dimension(:), allocatable, intent(inout) :: sweepbounds
!    integer, intent(inout) :: nosweeps
!
!    integer, dimension(:,:,:,:), allocatable :: streamcond
!    integer, dimension(:,:,:), allocatable :: need
!    integer, dimension(:), allocatable :: linneed
!    integer, dimension(:), allocatable :: linsweeplist
!    integer, dimension(:,:), allocatable :: B_sweeplist
!    integer, dimension(:), allocatable :: B_sweepbounds
!
!    !integer, dimension(:), allocatable :: dsweepsize
!    !logical, dimension(:), allocatable :: checkmatrix
!    !logical, dimension(:), allocatable :: B_checkmatrix
!    !integer :: falsecounter
!
!    allocate(sweeplist(0,3))
!    allocate(sweepbounds(1))
!    allocate(streamcond(Nmu,Nnu,NFe,NE))
!    allocate(need(Nmu,Nnu,NE))
!
!    !allocate(checkmatrix(0))
!
!    do e = 1, NE
!        do f = 1, NFe
!            do j = 1, Nnu
!                do i = 1, Nmu
!                    streamcond(i,j,f,e) = int(sign(1.,dot_product(khat(:,i,j),normal(:,f,e))))
!                    ! If wave is OUTGOING, streamcond = 1
!                    ! If wave is INCOMING, streamcond = -1
!                    ! If wave is PARALLEL, streamcond = 0
!                end do
!            end do
!        end do
!    end do
!    ! Make streamcond logical?
!
!    streamcond = reshape(streamcond, [Nmu,Nnu,NE,NFe], order = [1,2,4,3])
!
!    need(:,:,:) = -sum(streamcond(:,:,:,:)-1,4)/2
!
!    do n = 1, size(bdyfc,1)
!        need(:,:,bdyfc(n,1)) = need(:,:,bdyfc(n,1)) + (streamcond(:,:,bdyfc(n,1),bdyfc(n,2))-1)/2
!    end do
!    deallocate(streamcond)
!
!    need = reshape(need, [NE,Nmu,Nnu], order = [2,3,1])
!
!    allocate(linneed, source = pack(need,.true.))
!
!    ! IN PRINCIPLE:
!    ! need(e,i,j) = linneed(e+(i-1)*NE+(j-1)*NE*Nmu)
!
!    sweep = 1
!    sweepbounds(1) = 0
!    do while (size(sweeplist,1) < NE*Nmu*Nnu)
!        allocate(linsweeplist, source = pack([(ip,ip=1,NE*Nmu*Nnu)], linneed == 0))
!        ! Think of way to cut out linsweeplist
!
!        allocate(B_sweeplist(size(sweeplist,1)+size(linsweeplist),3))
!        B_sweeplist(1:size(sweeplist,1),:) = sweeplist
!
!        ! j(s)
!        B_sweeplist(size(sweeplist,1)+1:size(B_sweeplist,1),3) = 1+&
!        floor(real(linsweeplist-1)/(NE*Nmu))
!
!        ! i(s)
!        B_sweeplist(size(sweeplist,1)+1:size(B_sweeplist,1),2) = 1+floor(real(linsweeplist-1)/NE+&
!        (1-B_sweeplist(size(sweeplist,1)+1:size(B_sweeplist,1),3))*Nmu)
!
!        ! e(s)
!        B_sweeplist(size(sweeplist,1)+1:size(B_sweeplist,1),1) = linsweeplist+&
!        (1-B_sweeplist(size(sweeplist,1)+1:size(B_sweeplist,1),2))*NE+&
!        (1-B_sweeplist(size(sweeplist,1)+1:size(B_sweeplist,1),3))*NE*Nmu
!
!        call move_alloc(B_sweeplist,sweeplist)
!
!        allocate(B_sweepbounds(size(sweepbounds)+1))
!        B_sweepbounds(1:size(sweepbounds)) = sweepbounds
!        B_sweepbounds(size(sweepbounds)+1) = size(sweeplist,1)
!        call move_alloc(B_sweepbounds,sweepbounds)
!
!        !falsecounter = 0
!        !! Check that need(sweep) are all zero
!        !do n = sweepbounds(sweep)+1, sweepbounds(sweep+1)
!        !    if (need(sweeplist(n,1),sweeplist(n,2),sweeplist(n,3)) /= 0) then
!        !        falsecounter = falsecounter + 1
!        !    end if
!        !end do
!        !
!        !if (falsecounter == 0) then
!        !    allocate(B_checkmatrix(size(checkmatrix)+1))
!        !    B_checkmatrix(1:size(checkmatrix)) = checkmatrix(:)
!        !    B_checkmatrix(size(checkmatrix)+1) = .true.
!        !    call move_alloc(B_checkmatrix,checkmatrix)
!        !else if (falsecounter /= 0) then
!        !    allocate(B_checkmatrix(size(checkmatrix)+1))
!        !    B_checkmatrix(1:size(checkmatrix)) = checkmatrix(:)
!        !    B_checkmatrix(size(checkmatrix)+1) = .false.
!        !    call move_alloc(B_checkmatrix,checkmatrix)
!        !end if
!
!        do n = sweepbounds(sweep)+1, sweepbounds(sweep+1)
!            if (need(sweeplist(n,1),sweeplist(n,2),sweeplist(n,3)) == 0) &
!            need(sweeplist(n,1),sweeplist(n,2),sweeplist(n,3)) = -1
!            do f = 1, NFe
!                if (eprime(f,sweeplist(n,1)) /= 0 .and. &
!                need(eprime(f,sweeplist(n,1)),sweeplist(n,2),sweeplist(n,3)) > 0) then
!                    need(eprime(f,sweeplist(n,1)),sweeplist(n,2),sweeplist(n,3)) = &
!                    need(eprime(f,sweeplist(n,1)),sweeplist(n,2),sweeplist(n,3)) - 1
!                end if
!            end do
!        end do
!
!        !do f = 1, NFe
!        !    do n = 1, size(sweeplist,1)
!        !        if (need(sweeplist(n,1),sweeplist(n,2),sweeplist(n,3)) == 0) then
!        !            need(sweeplist(n,1),sweeplist(n,2),sweeplist(n,3)) = -1
!        !        end if
!        !        if (eprime(f,sweeplist(n,1)) /= 0 .and.&
!        !            need(eprime(f,sweeplist(n,1)),sweeplist(n,2),sweeplist(n,3)) > 0) then
!        !            need(eprime(f,sweeplist(n,1)),sweeplist(n,2),sweeplist(n,3)) = &
!        !            need(eprime(f,sweeplist(n,1)),sweeplist(n,2),sweeplist(n,3)) - 1
!        !        end if
!        !    end do
!        !end do
!
!        linneed = pack(need,.true.)
!
!        deallocate(linsweeplist)
!        sweep = sweep + 1
!        !if (size(sweeplist,1) == NE*Nmu*Nnu) exit
!        if (size(sweeplist,1) > NE*Nmu*Nnu) then
!            print *, "sweeplist exceeds max size. Program ended. Look for problem."
!            stop
!        end if
!    end do
!    nosweeps = size(sweepbounds)
!
!    !print *, nosweeps
!
!    allocate(B_sweepbounds(nosweeps+1))
!    B_sweepbounds(1:nosweeps) = sweepbounds(:)
!    B_sweepbounds(nosweeps+1) = NE*Nmu*Nnu
!    call move_alloc(B_sweepbounds,sweepbounds)
!
!    open(1, file = "Storage/Storage_7_"//selected_mesh//"sweeplist.dat", form = "unformatted")
!    write(1) sweeplist
!    close(1)
!
!    open(1, file = "Storage/Storage_7_"//selected_mesh//"sweepbounds.dat", form = "unformatted")
!    write(1) sweepbounds
!    close(1)
!
!    open(1, file = "Storage/Storage_7_"//selected_mesh//"nosweeps.dat", form = "unformatted")
!    write(1) nosweeps
!    close(1)
!
!    !allocate(dsweepsize(nosweeps))
!    !forall(sweep = 1:nosweeps)
!    !    dsweepsize(sweep) = sweepbounds(sweep+1)-sweepbounds(sweep)
!    !end forall
!
!    !print *, all(checkmatrix) -- > TRUE
!    !print *, (all(sweeplist > 0)) --> TRUE
!    !print *, (all(sweeplist <= NE)) --> TRUE
!    !print *, (all(need == -1)) --> TRUE
!
!    !allocate(streamcond(NE,Nmu,Nnu,NFe))
!    !allocate(need(NE,Nmu,Nnu))
!    !allocate(oldneed(NE,Nmu,Nnu))
!    !allocate(visitedcounter(NE,Nmu,Nnu))
!    !
!    !forall(e = 1:NE, i = 1:Nmu, j = 1:Nnu, f = 1:NFe)
!    !    streamcond(e,i,j,f) = int(sign(1.0,dot_product(khat(i,j,:),normal(e,f,:))))
!    !    ! If wave is OUTGOING, streamcond = 1
!    !    ! If wave is INCOMING, streamcond = -1
!    !    ! If wave is PARALLEL, streamcond = 0
!    !end forall
!    !! Make streamcond logical?
!    !
!    !need(:,:,:) = -sum((streamcond(:,:,:,:)-1)/2,4)
!    !
!    !do n = 1, size(bdyfc,1)
!    !   need(bdyfc(n,1),:,:) = need(bdyfc(n,1),:,:) + (streamcond(bdyfc(n,1),:,:,bdyfc(n,2))-1)/2
!    !end do
!
!    ! DIAGNOSES
!
!    ! need - PASSED
!    ! Nothing that has need = 0 has upstream faces other than boundary face - PASSED
!    !do e = 1, NE
!    !    do i = 1, Nmu
!    !        do j = 1, Nnu
!    !            if (need(e,i,j) == 0) then
!    !                do f = 1, NFe
!    !                    if (eprime(f,e) /= 0) then
!    !                        if (dot_product(normal(:,f,e),khat(:,i,j)) < 0) print *, "PROBLEM"
!    !                    end if
!    !                end do
!    !            end if
!    !        end do
!    !    end do
!    !end do
!
!    ! sweep
!    ! 1: Check that every element in the given sweep has needs = 0.
!    ! if .not. needs(e(s),i(s),j(s)) == 0
!    ! 2: Check that every element and face in the sweep is neighboring the elements in the previous sweep.
!    ! if any e(s) is not in prevsweepelist
!    ! 3: Reduce needs of next sweep.
!    ! need(eprime(e(s+1),f),i(s+1),j(s+1)) =-1
!    ! 4: Check that needs of next sweep was reduced by as many faces as are shared with this sweep.
!    ! in previous needs loop, do visitedcounter = 0, then visitedcounter(e,i,j) = visitedcounter(e,i,j) + 1.
!    ! This gives the number of times the SAME element was visited. Ensure that new needs is needs - visitedcounter.
!    ! Check that visitedcounter is in fact the number of shared faces between an element and the els of prev sweep???
!
!    !allocate(B_needlogic(nosweeps))
!    !
!    !!allocate(needlogic(sweepbounds(sweep+1)-sweepbounds(sweep)))
!    !!forall(n = sweepbounds(sweep)+1:sweepbounds(sweep+1))
!    !!    needlogic(n-sweepbounds(sweep)) = need(sweeplist(n,1),&
!    !!    sweeplist(n,2),&
!    !!    sweeplist(n,3)) == 0
!    !!end forall
!    !
!    !do sweep = 1, 1
!    !    allocate(needlogic(sweepbounds(sweep+1)-sweepbounds(sweep)))
!    !    forall(n = sweepbounds(sweep)+1:sweepbounds(sweep+1))
!    !        needlogic(n-sweepbounds(sweep)) = need(sweeplist(n,1),&
!    !        sweeplist(n,2),&
!    !        sweeplist(n,3)) == 0
!    !    end forall
!    !
!    !    B_needlogic(sweep) = all(needlogic(:))
!    !
!    !    deallocate(needlogic)
!    !
!    !    !do f = 1, NFe
!    !    !    do n = sweepbounds(sweep)+1, sweepbounds(sweep+1)
!    !    !        if (eprime(sweeplist(n,1),f) == 0) cycle
!    !    !        need(eprime(sweeplist(n,1),f),&
!    !    !                    sweeplist(n,2),&
!    !    !                    sweeplist(n,3)) = &
!    !    !        need(eprime(sweeplist(n,1),f),&
!    !    !                    sweeplist(n,2),&
!    !    !                    sweeplist(n,3)) - 1
!    !    !        visitedcounter(eprime(sweeplist(n,1),f),&
!    !    !        sweeplist(n,2),&
!    !    !        sweeplist(n,3)) = &
!    !    !        visitedcounter(eprime(sweeplist(n,1),f),&
!    !    !        sweeplist(n,2),&
!    !    !        sweeplist(n,3)) + 1
!    !    !    end do
!    !    !end do
!    !    !if (.not. all(need(:,:,:) == oldneed(:,:,:) - visitedcounter(:,:,:))) print *, "PROBLEM 3"
!    !end do
!    !
!    !print *, all(B_needlogic)
!end subroutine construct_sweeplist

!subroutine read_sweeplist(sweeplist, sweepbounds, nosweeps)
!    implicit none
!    integer, dimension(:,:), allocatable, intent(inout) :: sweeplist
!    integer, dimension(:), allocatable, intent(inout) :: sweepbounds
!    integer :: nosweeps
!
!    allocate(sweeplist(NE*Nmu*Nnu,3))
!
!    open(1, file = "Storage/Storage_7_"//selected_mesh//"sweeplist.dat", form = "unformatted", action = "read")
!    read(1) sweeplist
!    close(1)
!
!    open(1, file = "Storage/Storage_7_"//selected_mesh//"nosweeps.dat", form = "unformatted", action = "read")
!    read(1) nosweeps
!    close(1)
!
!    allocate(sweepbounds(nosweeps+1))
!
!    open(1, file = "Storage/Storage_7_"//selected_mesh//"sweepbounds.dat", form = "unformatted", action = "read")
!    read(1) sweepbounds
!    close(1)
!end subroutine read_sweeplist

!subroutine construct_esweeplist(bdyfc, eprime, normal, khat, nosweeps, esweeplist, esweepbounds, enosweeps)
!    implicit none
!    integer, dimension(:,:), intent(in) :: bdyfc
!    integer, dimension(:,:), intent(in) :: eprime
!    real, dimension(:,:,:), intent(in) :: normal
!    real, dimension(:,:,:), intent(in) :: khat
!
!    integer, intent(inout) ::  nosweeps
!    integer, dimension(:,:,:), allocatable, intent(inout) :: esweeplist
!    integer, dimension(:,:,:), allocatable, intent(inout) :: esweepbounds
!    integer, dimension(:,:), allocatable, intent(inout) :: enosweeps
!
!    integer, dimension(:,:), allocatable :: sweeplist
!    integer, dimension(:), allocatable :: sweepbounds
!    integer, dimension(:), allocatable :: sweeplistindices
!    integer :: starting
!    integer :: ending
!
!    !integer, dimension(:), allocatable :: iplist
!
!    call construct_sweeplist(bdyfc, eprime, normal, khat, sweeplist, sweepbounds, nosweeps)
!
!    !allocate(iplist(NE*Nmu*Nnu))
!    !
!    !iplist = 0
!    !
!    !do ip = 1, NE*Nmu*Nnu
!    !    iplist(sweeplist(ip,1)+(sweeplist(ip,2)-1)*NE+(sweeplist(ip,3)-1)*NE*Nmu) = &
!    !    sweeplist(ip,1)+(sweeplist(ip,2)-1)*NE+(sweeplist(ip,3)-1)*NE*Nmu
!    !end do
!    !
!    !if (any(iplist(:) == 0)) print *, "ZERO FOUND"
!
!    allocate(esweeplist(NE,Nmu,Nnu))
!    allocate(esweepbounds(Nmu,Nnu,nosweeps+1))
!    allocate(enosweeps(Nmu,Nnu))
!
!    esweepbounds(:,:,:) = 0
!
!    do j = 1, Nnu
!        thisloop: do i = 1, Nmu
!            starting = 0
!            ending = 0
!            do sweep = 1, nosweeps
!                allocate(sweeplistindices, source = pack([(k,k=sweepbounds(sweep)+1,sweepbounds(sweep+1))],&
!                ((sweeplist(sweepbounds(sweep)+1:sweepbounds(sweep+1),2) == i) .and. &
!                (sweeplist(sweepbounds(sweep)+1:sweepbounds(sweep+1),3) == j))))
!
!                starting = starting + ending
!                ending = size(sweeplistindices)
!
!                if (ending == 0) then
!                    enosweeps(i,j) = sweep - 1
!                    deallocate(sweeplistindices)
!                    cycle thisloop
!                end if
!
!                esweeplist(starting+1:starting+ending,i,j) = sweeplist(sweeplistindices(:),1)
!
!                esweepbounds(i,j,sweep+1) = starting+ending
!
!                deallocate(sweeplistindices)
!            end do
!        end do thisloop
!    end do
!
!    do k = 1, NE
!        if (.not. any(esweeplist(:,14,12) == k)) then
!            print *, .false., k
!        end if
!    end do
!
!    open(1, file = "Storage/Storage_7_nosweeps.dat", form = "unformatted")
!    write(1) nosweeps
!    close(1)
!
!    open(1, file = "Storage/Storage_7_esweeplist.dat", form = "unformatted")
!    write(1) esweeplist
!    close(1)
!
!    open(1, file = "Storage/Storage_7_esweepbounds.dat", form = "unformatted")
!    write(1) esweepbounds
!    close(1)
!
!    open(1, file = "Storage/Storage_7_enosweeps.dat", form = "unformatted")
!    write(1) enosweeps
!    close(1)
!end subroutine construct_esweeplist
!
!subroutine read_esweeplist(nosweeps, esweeplist, esweepbounds, enosweeps)
!    implicit none
!    integer, intent(inout) :: nosweeps
!    integer, dimension(:,:,:), allocatable, intent(inout) :: esweeplist
!    integer, dimension(:,:,:), allocatable, intent(inout) :: esweepbounds
!    integer, dimension(:,:), allocatable, intent(inout) :: enosweeps
!
!    open(1, file = "Storage/Storage_7_nosweeps.dat", form = "unformatted", action = "read")
!    read(1) nosweeps
!    close(1)
!
!    allocate(esweeplist(Nmu,Nnu,NE))
!    allocate(esweepbounds(Nmu,Nnu,nosweeps))
!    allocate(enosweeps(Nmu,Nnu))
!
!    open(1, file = "Storage/Storage_7_esweeplist.dat", form = "unformatted", action = "read")
!    read(1) esweeplist
!    close(1)
!
!    open(1, file = "Storage/Storage_7_esweepbounds.dat", form = "unformatted", action = "read")
!    read(1) esweepbounds
!    close(1)
!
!    open(1, file = "Storage/Storage_7_enosweeps.dat", form = "unformatted", action = "read")
!    read(1) enosweeps
!    close(1)
!end subroutine read_esweeplist
!
!subroutine sweep_mesh_restructure_map(Cekk, sweeplist, eprime, Cipkk, eijip, ipprime)
!    implicit none
!    integer, dimension(:,:), intent(in) :: Cekk
!    integer, dimension(:,:), intent(in) :: sweeplist
!    integer, dimension(:,:), intent(in) :: eprime
!
!    integer, dimension(:,:), allocatable, intent(inout) :: Cipkk
!    integer, dimension(:,:,:), allocatable, intent(inout) :: eijip
!    integer, dimension(:,:), allocatable, intent(inout) :: ipprime
!
!    allocate(Cipkk(NKe,NE*Nmu*Nnu))
!    allocate(eijip(NE,Nmu,Nnu))
!    allocate(ipprime(NFe,NE*Nmu*Nnu))
!
!    do k = 1, NKe
!        do ip = 1, NE*Nmu*Nnu
!            Cipkk(k,ip) = Cekk(sweeplist(ip,1),k)
!        end do
!    end do
!
!    do ip = 1, NE*Nmu*Nnu
!        eijip(sweeplist(ip,1),sweeplist(ip,2),sweeplist(ip,3)) = ip
!    end do
!
!    do ip = 1, NE*Nmu*Nnu
!        do f = 1, NFe
!            if (eprime(f,sweeplist(ip,1)) == 0) then
!                ipprime(f,ip) = 0
!                cycle
!            end if
!            ipprime(f,ip) = eijip(eprime(f,sweeplist(ip,1)),sweeplist(ip,2),sweeplist(ip,3))
!        end do
!    end do
!end subroutine sweep_mesh_restructure_map
!
!subroutine esweep_mesh_restructure_map(eprime, esweeplist, enew, eprimeij)
!    implicit none
!    integer, dimension(:,:), intent(in) :: eprime
!    integer, dimension(:,:,:), intent(in) :: esweeplist
!
!    integer, dimension(:,:,:), allocatable, intent(inout) :: enew
!    integer, dimension(:,:,:,:), allocatable, intent(inout) :: eprimeij
!
!    allocate(enew(NE,Nmu,Nnu))
!    allocate(eprimeij(NFe,NE,Nmu,Nnu))
!
!    do j = 1, Nnu
!        do i = 1, Nmu
!            do ip = 1, NE
!                enew(esweeplist(ip,i,j),i,j) = ip
!            end do
!        end do
!    end do
!
!    do j = 1, Nnu
!        do i = 1, Nmu
!            do e = 1, NE
!                do f = 1, NFe
!                    if (eprime(f,e) == 0) then
!                        eprimeij(f,enew(e,i,j),i,j) = 0
!                    else
!                        eprimeij(f,enew(e,i,j),i,j) = enew(eprime(f,e),i,j)
!                    end if
!                end do
!            end do
!        end do
!    end do
!end subroutine esweep_mesh_restructure_map
!
!subroutine varsigmas(normal, khat, eijip, varsigmaup, restruc_varsigmadown)
!    implicit none
!    real, dimension(:,:,:), intent(in) :: normal
!    real, dimension(:,:,:), intent(in) :: khat
!    integer, dimension(:,:,:), intent(in) :: eijip
!
!    real, dimension(:,:,:,:), allocatable, intent(inout) :: varsigmaup
!    real, dimension(:,:), allocatable, intent(inout) :: restruc_varsigmadown
!
!    real, dimension(:,:,:,:), allocatable :: varsigmadown
!
!    allocate(varsigmaup(NFe,NE,Nmu,Nnu))
!    allocate(restruc_varsigmadown(NFe,NE*Nmu*Nnu))
!    allocate(varsigmadown(NFe,NE,Nmu,Nnu))
!
!    do j = 1, Nnu
!        do i = 1, Nmu
!            do e = 1, NE
!                do f = 1, NFe
!                    varsigmaup(f,e,i,j) = (1./2)*(dot_product(khat(:,i,j),normal(:,f,e)) + &
!                        abs(dot_product(khat(:,i,j),normal(:,f,e))))
!                    varsigmadown(f,e,i,j) = (1./2)*(dot_product(khat(:,i,j),normal(:,f,e)) - &
!                        abs(dot_product(khat(:,i,j),normal(:,f,e))))
!                end do
!            end do
!        end do
!    end do
!
!    do j = 1, Nnu
!        do i = 1, Nmu
!            restruc_varsigmadown(:,eijip(:,i,j)) = varsigmadown(:,:,i,j)
!        end do
!    end do
!
!    !do e = 1, NE
!    !    do i = 1, Nmu
!    !        do j = 1, Nnu
!    !            do f = 1, NFe
!    !                if (varsigmaup(e,i,j,f) /= 0.0 .and. varsigmadown(e,i,j,f) /= 0.0) print *, "PROBLEM 1"
!    !                if (varsigmaup(e,i,j,f) == 0.0 .and. varsigmadown(e,i,j,f) == 0.0) then
!    !                    if (dot_product(khat(i,j,:),normal(e,f,:)) /= 0) print *, "PROBLEM 2"
!    !                end if
!    !                if (varsigmaup(e,i,j,f) < 0) print *, "PROBLEM 3"
!    !                if (varsigmadown(e,i,j,f) > 0) print *, "PROBLEM 4"
!    !            end do
!    !        end do
!    !    end do
!    !end do
!end subroutine varsigmas

subroutine construct_esweeplist(eprime, fprime, normal, khat, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa)
    implicit none
    integer, dimension(:,:), intent(in) :: eprime
    integer, dimension(:,:), intent(in) :: fprime
    real, dimension(:,:,:), intent(in) :: normal
    real, dimension(:,:,:), intent(in) :: khat

    integer, dimension(:,:,:), allocatable, intent(inout) :: esweeplist
    integer, dimension(:,:,:), allocatable, intent(inout) :: esweepbounds
    !integer, dimension(:,:), allocatable, intent(inout) :: enosweeps
    integer, dimension(:,:,:), allocatable, intent(inout) :: enew
    !integer, dimension(:,:,:,:), allocatable, intent(inout) :: epprime
    integer, dimension(:,:,:,:), allocatable, intent(inout) :: fold
    integer, dimension(:,:,:), allocatable, intent(inout) :: Ia
    integer, dimension(:,:,:), allocatable, intent(inout) :: Ja
    integer, intent(inout) :: noJa

    integer :: e, f, i, ip, j, sweep
    integer, dimension(:,:), allocatable :: enosweeps
    integer, dimension(:,:), allocatable :: esparents
    integer, dimension(:,:), allocatable :: dupeesparents
    integer, dimension(:), allocatable :: zerorows
    integer, dimension(:), allocatable :: localesweeplist
    integer, dimension(:), allocatable :: localesweepbounds
    integer, dimension(:), allocatable :: tempenew
    integer, dimension(:,:), allocatable :: tempepprime
    integer, dimension(:,:), allocatable :: tempfold
    integer, dimension(:), allocatable :: tempIa
    integer, dimension(:), allocatable :: tempJa
    integer :: l_count
    integer :: zeroloc
    integer :: l_minloc
    integer :: midloc
    integer :: l_maxloc
    logical, dimension(:), allocatable :: logarray
    integer :: runningtotal
    integer :: SS0
    integer :: SS1
    integer, dimension(:,:,:), allocatable :: B_esweepbounds
    integer, dimension(:,:,:), allocatable :: B_Ja
    character(100) :: fname
    character(50) :: Nmustr
    character(50) :: Nphistr
    character(100) :: folder
    character(100) :: command
    integer :: status

    allocate(esweeplist(NE,Nmu,Nphi))
    allocate(esweepbounds(NE,Nmu,Nphi))
    allocate(enosweeps(Nmu,Nphi))
    allocate(enew(NE,Nmu,Nphi))
    !allocate(epprime(NFe,NE,Nmu,Nphi))
    allocate(fold(NFe,NE,Nmu,Nphi))
    allocate(Ia(NE+1,Nmu,Nphi))
    allocate(Ja(NFe*NE,Nmu,Nphi))
    allocate(esparents(NFe,NE))
    allocate(dupeesparents(NFe,NE))
    allocate(localesweeplist(NE))
    allocate(localesweepbounds(NE))
    allocate(tempenew(0:NE))
    allocate(tempepprime(NFe,NE))
    allocate(tempfold(NFe,NE))
    allocate(tempIa(NE+1))
    allocate(tempJa(NFe*NE))
    allocate(logarray(NFe))

    esweepbounds = 0
    Ja = 0

    ! PARALLELIZATION DOESN'T GIVE CORRECT SWEEP... FIX??
    !!$OMP parallel do collapse(2) shared(esweeplist, esweepbounds, enosweeps, enew, fold, Ia, Ja, noJa) &
    !!$OMP private(f, e, i, j, ip, sweep) &
    !!$OMP private(esparents, dupeesparents, SS0, SS1, zerorows, localesweepbounds, localesweeplist) &
    !!$OMP private(tempepprime, tempenew, tempfold, tempIa, tempJa, runningtotal, logarray, l_count, l_maxloc, l_minloc, midloc)
    do j = 1, Nphi
        do i = 1, Nmu
            dupeesparents = 0
            do e = 1, NE
                do f = 1, NFe
                    if (dot_product(khat(1:3,i,j),normal(1:3,f,e)) < 0) then
                        dupeesparents(f,e) = eprime(f,e)
                    end if
                    ! DO A LOGARRAY = normal(:,f,e)dot tempkhat(:) < 0 and use that somehow?
                end do
            end do
            esparents = dupeesparents

            sweep = 1
            localesweepbounds = 0
            localesweeplist = 0
            do while (localesweeplist(NE) .eq. 0)
                ! FIND ALL ROWS WHICH ARE WHOLLY ZERO:
                allocate(zerorows, source = pack([(e,e=1,NE)], mask = sum(dupeesparents,1) .eq. 0))
                localesweepbounds(sweep+1) = localesweepbounds(sweep)+size(zerorows)
                if (size(zerorows) .eq. 0) then
                    print *, "PROBLEM: MESSAGE WIP"
                    print *, "sweep is occuring with no elements. WHY???"
                    stop
                end if
                localesweeplist(localesweepbounds(sweep)+1:localesweepbounds(sweep+1)) = zerorows
                deallocate(zerorows) ! SLOW?

                SS0 = localesweepbounds(sweep)+1
                SS1 = localesweepbounds(sweep+1)
                dupeesparents(1:NFe,localesweeplist(SS0:SS1)) = -1
                do ip = SS0, SS1
                    do f = 1, NFe
                        if (eprime(f,localesweeplist(ip)) .eq. 0) cycle ! SLOW? Have a garbage index in dupeesparents and just allow it to be overwritten over and over? would need to do like sum(dupeesparents(1:NFe,1:NE),1) if there is a 0
                        if (dupeesparents(fprime(f,localesweeplist(ip)),eprime(f,localesweeplist(ip))) .eq. -1) cycle ! Should I do this?
                        dupeesparents(fprime(f,localesweeplist(ip)),eprime(f,localesweeplist(ip))) = 0
                        ! Note that this may interfere with setting dupeesparents = -1. Is that a problem for me later? See if () cycle above for solution
                    end do
                end do

                sweep = sweep + 1

                !print *, SS0, SS1

                if (sweep .gt. NE) then ! I think condition is SS1 or SS0 .ge. NE, not sweep.
                    print *, "STOP!"
                    print *, "MODULE 5, SUBROUTINE construct_esweeplist:"
                    print *, "Error while constructing sweep order. The number of sweep steps"
                    print *, "exceeds the number of elements in the mesh. This means that the"
                    print *, "sweep constructed is nonsense, and not the correct sweep order."
                    print *, "Double check the provided mesh and mesh file."
                    print *, "This error is occuring at the discrete ordinate indices:"
                    print *, [i,j]
                    print *, "As a reminder, as of this version wiscobolt does not support"
                    print *, "concave geometries nor elements other than tetrahedral elements."
                    print *, "PROGRAM ENDING."
                    stop
                end if
            end do

            esweeplist(1:NE,i,j) = localesweeplist(1:NE)
            esweepbounds(1:sweep,i,j) = localesweepbounds(1:sweep)
            enosweeps(i,j) = sweep - 1

            tempenew = 0
            do e = 1, NE
                tempenew(localesweeplist(e)) = e
            end do

            !do f = 1, NFe
            !    tempepprime(f,tempenew(1:NE)) = tempenew(esparents(f,1:NE))
            !end do
            do e = 1, NE
                do f = 1, NFe
                    tempepprime(f,tempenew(e)) = tempenew(esparents(f,e))
                end do
            end do

            SS0 = localesweepbounds(2)+1
            tempfold(1:NFe,1:SS0-1) = 0
            tempIa(1:SS0-1) = 0
            tempIa(SS0) = 1
            tempJa = 0
            runningtotal = 1
            do e = SS0, NE
                l_count = count(tempepprime(1:NFe,e) .eq. 0)
                if (l_count .eq. 1) then
                    zeroloc = findloc(tempepprime(1:NFe,e), 0, dim = 1)
                    logarray = .true.
                    logarray(zeroloc) = .false.
                    l_maxloc = maxloc(tempepprime(1:NFe,e), dim = 1)
                    l_minloc = minloc(tempepprime(1:NFe,e), dim = 1, mask = logarray)
                    midloc = 10-l_maxloc-l_minloc-zeroloc
                    tempfold(1,e) = l_minloc
                    tempfold(2,e) = midloc
                    tempfold(3,e) = l_maxloc
                    tempfold(4,e) = zeroloc
                    tempepprime(1:NFe,e) = [tempepprime(l_minloc,e),tempepprime(midloc,e),tempepprime(l_maxloc,e),0]
                    tempJa(runningtotal:runningtotal+2) = tempepprime(1:3,e)
                    runningtotal = runningtotal + 3
                    tempIa(e+1) = runningtotal
                else if (l_count .eq. 2) then
                    logarray = tempepprime(1:NFe,e) /= 0
                    l_maxloc = maxloc(tempepprime(1:NFe,e), dim = 1)
                    l_minloc = minloc(tempepprime(1:NFe,e), dim = 1, mask = logarray)
                    tempfold(1,e) = l_minloc
                    tempfold(2,e) = l_maxloc
                    tempfold(3,e) = 0
                    tempfold(4,e) = 0
                    tempepprime(1:NFe,e) = [tempepprime(l_minloc,e),tempepprime(l_maxloc,e),0,0]
                    tempJa(runningtotal:runningtotal+1) = tempepprime(1:2,e)
                    runningtotal = runningtotal + 2
                    tempIa(e+1) = runningtotal
                else
                    l_maxloc = maxloc(tempepprime(1:NFe,e), dim = 1)
                    tempfold(1,e) = l_maxloc
                    tempfold(l_maxloc,e) = 1
                    tempepprime(1:NFe,e) = [tempepprime(l_maxloc,e),0,0,0]
                    tempJa(runningtotal) = tempepprime(1,e)
                    runningtotal = runningtotal + 1
                    tempIa(e+1) = runningtotal
                end if
            end do

            enew(1:NE,i,j) = tempenew(1:NE)
            !epprime(1:NFe,1:NE,i,j) = tempepprime(1:NFe,1:NE)
            fold(1:NFe,1:NE,i,j) = tempfold
            Ia(1:NE+1,i,j) = tempIa
            noJa = runningtotal-1
            Ja(1:noJa,i,j) = tempJa(1:noJa)
        end do
    end do
    !!$OMP end parallel do

    l_count = maxval(enosweeps)+1
    allocate(B_esweepbounds(l_count,Nmu,Nphi))
    B_esweepbounds = esweepbounds(1:l_count,1:Nmu,1:Nphi)
    call move_alloc(B_esweepbounds, esweepbounds)

    allocate(B_Ja(noJa,Nmu,Nphi))
    B_Ja = Ja(1:noJa,1:Nmu,1:Nphi)
    call move_alloc(B_Ja, Ja)

    write(Nmustr,*) Nmu
    write(Nphistr,*) Nphi

    if (save_sweep) then
        ! Check if folder exists for this Nmu Nphi pair. If not, create one
        folder = "Mesh files/"//trim(adjustl(selected_mesh))//"/sweep sets/"//trim(adjustl(Nmustr))//"x"//trim(adjustl(Nphistr))
        command = 'dir "'//trim(folder)//'" > NUL 2>&1'
        status = system(command)
        if (status .ne. 0) then
            command = 'mkdir "'//trim(folder)//'" > NUL 2>&1'
            status = system(command)
        end if

        fname = trim(adjustl(selected_mesh))//"/sweep sets/"//trim(adjustl(Nmustr))//"x"//trim(adjustl(Nphistr))

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/esweeplist.dat", form = "unformatted")
        write(1) esweeplist
        close(1)

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/esweepbounds.dat", form = "unformatted")
        write(1) esweepbounds
        close(1)

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/enew.dat", form = "unformatted")
        write(1) enew
        close(1)

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/fold.dat", form = "unformatted")
        write(1) fold
        close(1)

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/Ia.dat", form = "unformatted")
        write(1) Ia
        close(1)

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/Ja.dat", form = "unformatted")
        write(1) Ja
        close(1)

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/noJa.dat", form = "unformatted")
        write(1) noJa
        close(1)

        open(1, file = "Mesh files/"//trim(adjustl(fname))//"/maxenosweeps.dat", form = "unformatted")
        write(1) l_count
        close(1)
    end if
end subroutine construct_esweeplist

subroutine read_esweeplist(esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa)
    implicit none
    integer, dimension(:,:,:), allocatable, intent(inout) :: esweeplist
    integer, dimension(:,:,:), allocatable, intent(inout) :: esweepbounds
    integer, dimension(:,:,:), allocatable, intent(inout) :: enew
    integer, dimension(:,:,:,:), allocatable, intent(inout) :: fold
    integer, dimension(:,:,:), allocatable, intent(inout) :: Ia
    integer, dimension(:,:,:), allocatable, intent(inout) :: Ja
    integer, intent(inout) :: noJa

    integer :: maxenosweeps
    character(100) :: fname
    character(50) :: Nmustr
    character(50) :: Nphistr

    write(Nmustr,*) Nmu
    write(Nphistr,*) Nphi
    fname = trim(adjustl(selected_mesh))//"/sweep sets/"//trim(adjustl(Nmustr))//"x"//trim(adjustl(Nphistr))

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/maxenosweeps.dat", form = "unformatted", action = "read")
    read(1) maxenosweeps
    close(1)

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/noJa.dat", form = "unformatted", action = "read")
    read(1) noJa
    close(1)

    allocate(esweeplist(NE,Nmu,Nphi))
    allocate(esweepbounds(maxenosweeps,Nmu,Nphi))
    allocate(enew(NE,Nmu,Nphi))
    allocate(fold(NFe,NE,Nmu,Nphi))
    allocate(Ia(NE+1,Nmu,Nphi))
    allocate(Ja(noJa,Nmu,Nphi))

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/esweeplist.dat", form = "unformatted", action = "read")
    read(1) esweeplist
    close(1)

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/esweepbounds.dat", form = "unformatted", action = "read")
    read(1) esweepbounds
    close(1)

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/enew.dat", form = "unformatted", action = "read")
    read(1) enew
    close(1)

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/fold.dat", form = "unformatted", action = "read")
    read(1) fold
    close(1)

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/Ia.dat", form = "unformatted", action = "read")
    read(1) Ia
    close(1)

    open(1, file = "Mesh files/"//trim(adjustl(fname))//"/Ja.dat", form = "unformatted", action = "read")
    read(1) Ja
    close(1)
end subroutine read_esweeplist

subroutine construct_varsigmas(normal, khat, esweeplist, fold, Ia, noJa, varsigmaup, varsigmadown)
    implicit none
    real, dimension(:,:,:), intent(in) :: normal
    real, dimension(:,:,:), intent(in) :: khat
    integer, dimension(:,:,:), intent(in) :: esweeplist
    integer, dimension(:,:,:,:), intent(in) :: fold
    integer, dimension(:,:,:), intent(in) :: Ia
    integer, intent(in) :: noJa

    real, dimension(:,:,:,:), allocatable, intent(inout) :: varsigmaup
    real, dimension(:,:,:), allocatable, intent(inout) :: varsigmadown

    integer :: e, f, i, j
    integer, dimension(:), allocatable :: l_esweeplist
    integer, dimension(:,:), allocatable :: l_fold
    integer, dimension(:), allocatable :: l_Ia
    integer :: start
    integer :: SS0
    real :: l_dotprod

    allocate(varsigmaup(NFe,NE,Nmu,Nphi))
    allocate(varsigmadown(noJa,Nmu,Nphi))
    allocate(l_esweeplist(NE))
    allocate(l_fold(NFe,NE))
    allocate(l_Ia(NE+1))

    varsigmadown = 0

    !$OMP parallel do collapse(2) shared(varsigmadown, varsigmaup) private(j, i, e, f) &
    !$OMP private(l_esweeplist, l_fold, l_Ia, start, SS0, l_dotprod)
    do j = 1, Nphi
        do i = 1, Nmu
            l_esweeplist = esweeplist(1:NE,i,j)
            l_fold = fold(1:NFe,1:NE,i,j)
            l_Ia = Ia(1:NE+1,i,j)
            start = findloc(l_Ia > 0, .true., dim = 1) ! USE SWEEPBOUNDS INSTEAD
            do e = start, NE
                SS0 = l_Ia(e+1)-l_Ia(e)
                do f = 1, SS0
                    varsigmadown(l_Ia(e)-1+f,i,j) = dot_product(normal(1:3,l_fold(f,e),l_esweeplist(e)),khat(1:3,i,j))
                end do
            end do
            do e = 1, NE
                do f = 1, NFe
                    l_dotprod = dot_product(khat(1:3,i,j),normal(1:3,f,e))
                    varsigmaup(f,e,i,j) = (l_dotprod + abs(l_dotprod))/2
                end do
            end do
        end do
    end do
    !$OMP end parallel do
end subroutine construct_varsigmas

subroutine esweeplist_verification(r, eprime, normal, khat, esweeplist)
    implicit none
    real, dimension(:,:,:), intent(in) :: r
    integer, dimension(:,:), intent(in) :: eprime
    real, dimension(:,:,:), intent(in) :: normal
    real, dimension(:,:,:), intent(in) :: khat
    integer, dimension(:,:,:), intent(in) :: esweeplist

    integer :: f, i, ip, j, k
    integer, dimension(:), allocatable :: l_esweeplist
    integer :: l_e
    integer :: l_eprime
    integer :: localindex
    real, dimension(:,:,:,:), allocatable :: dotprods
    integer, dimension(:), allocatable :: checked
    integer, dimension(:), allocatable :: B_checked
    real :: l_dotprod

    allocate(l_esweeplist(NE))
    !allocate(dotprods(NFe,NE,Nmu,Nphi))

    !do j = 1, Nphi
    !    do i = 1, Nmu
    !        do e = 1, NE
    !            do f = 1, NFe
    !                dotprods(f,e,i,j) = dot_product(khat(1:3,i,j),normal(1:3,f,e))
    !            end do
    !        end do
    !    end do
    !end do

    !do j = 1, Nphi
    !    do i = 1, Nmu
    !        l_esweeplist = esweeplist(1:NE,i,j)
    !        do ip = 1, NE
    !            l_e = l_esweeplist(ip)
    !            face_iter: do f = 1, NFe
    !                if (dot_product(khat(1:3,i,j),normal(1:3,f,l_e)) < 0) then
    !                    l_eprime = eprime(f,l_e)
    !                    if (l_eprime == 0) cycle face_iter
    !                    localindex = findloc(l_esweeplist(1:ip), l_eprime, dim = 1)
    !                    if (localindex == 0) then
    !                        print *, "Error found in sweeplist."
    !                        print *, "Direction:", [i,j]
    !                        print *, "Allegedly downstream element and face:", [l_e,f]
    !                        print *, "Location in sweep:", ip
    !                        print *, "Allegedly upstream element:", l_eprime
    !                        print *, "Angle cosine between direction and face normal:", dot_product(khat(1:3,i,j),normal(1:3,f,l_e))
    !                        print *, "Coordinates of allegedly downstream element:"
    !                        do k = 1, NKe
    !                            print *, r(:,k,l_e)
    !                        end do
    !                        print *, "Coordinates of allegedly upstream element:"
    !                        do k = 1, NKe
    !                            print *, r(:,k,l_eprime)
    !                        end do
    !                        print *, "Coordinates of sweep direction:"
    !                        print *, khat(:,i,j)
    !                        print *, "Sweeplist up to errant element:"
    !                        print *, l_esweeplist(1:ip)
    !                        print *, "Check whether the allegedly upstream element is."
    !                        stop
    !                    end if
    !                end if
    !            end do face_iter
    !        end do
    !    end do
    !end do

    do j = 1, Nphi
        do i = 1, Nmu
            l_esweeplist = esweeplist(1:NE,i,j)
            allocate(checked(0))
            do ip = 1, NE
                l_e = l_esweeplist(ip)
                face_iter2: do f = 1, NFe
                    l_dotprod = dot_product(khat(1:3,i,j),normal(1:3,f,l_e))
                    if (l_dotprod < 0) then
                        l_eprime = eprime(f,l_e)
                        if (l_eprime == 0) cycle face_iter2
                        localindex = findloc(checked, l_eprime, dim = 1)
                        if (localindex == 0) then
                            print *, "ERROR", [j,i,l_e]
                            stop
                        end if
                    end if
                    allocate(B_checked(size(checked)+1))
                    B_checked(1:size(checked)) = checked
                    B_checked(size(checked)+1) = l_e
                    call move_alloc(B_checked, checked)
                end do face_iter2
            end do
            deallocate(checked)
        end do
    end do

    ! OR: visit each element, find all of its upstream elements, check if they exist in a particular stored array,
    ! and if their upstream elements exist in this array, write this element down, proceed.

    print *, "---", "MODULE 5: esweeplist_verification: No issues found."
end subroutine esweeplist_verification

!subroutine restructure_mesh(esweeplist, eprime, varsigmadown, Tprimepinv, Fsweep, enew, eprimeij)
!    implicit none
!    integer, dimension(:,:,:), intent(in) :: esweeplist
!    integer, dimension(:,:), intent(in) :: eprime
!
!    real, dimension(:,:,:,:), allocatable, intent(inout) :: varsigmadown
!    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Tprimepinv
!    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Fsweep
!    integer, dimension(:,:,:), allocatable, intent(inout) :: enew
!    integer, dimension(:,:,:,:), allocatable, intent(inout) :: eprimeij
!
!    real, dimension(:,:,:,:), allocatable :: B_varsigmadown
!    real, dimension(:,:,:,:,:,:), allocatable :: B_Tprimepinv
!    real, dimension(:,:,:,:,:,:), allocatable :: B_Fsweep
!
!    allocate(enew(NE,Nmu,Nnu))
!    allocate(eprimeij(NE,NFe,Nmu,Nnu))
!    allocate(B_varsigmadown(NE,Nmu,Nnu,NFe))
!    allocate(B_Tprimepinv(NE,Gp,Nmu,Nnu,NKe,NKe))
!    allocate(B_Fsweep(NE,Nmu,Nnu,NKe,NKe,NFe))
!
!    do j = 1, Nnu
!        do i = 1, Nmu
!            do ip = 1, NE
!                enew(esweeplist(i,j,ip),i,j) = ip
!            end do
!        end do
!    end do
!
!    do f = 1, NFe
!        do j = 1, Nnu
!            do i = 1, Nmu
!                do e = 1, NE
!                    if (eprime(e,f) == 0) then
!                        eprimeij(enew(e,i,j),f,i,j) = 0
!                    else
!                        eprimeij(enew(e,i,j),f,i,j) = enew(eprime(e,f),i,j)
!                    end if
!                end do
!            end do
!        end do
!    end do
!
!    do j = 1, Nnu
!        do i = 1, Nmu
!            do e = 1, NE
!                B_Tprimepinv(enew(e,i,j),:,i,j,:,:) = Tprimepinv(e,:,i,j,:,:)
!                B_Fsweep(enew(e,i,j),i,j,:,:,:) = Fsweep(e,i,j,:,:,:)
!                B_varsigmadown(enew(e,i,j),i,j,:) = varsigmadown(e,i,j,:)
!            end do
!        end do
!    end do
!
!    call move_alloc(B_Tprimepinv, Tprimepinv)
!    call move_alloc(B_Fsweep, Fsweep)
!    call move_alloc(B_varsigmadown, varsigmadown)
!end subroutine restructure_mesh
!
!subroutine lin_restructure_mesh&
!(varsigmadown, eprime, Tprimepinv, sweeplist, Fsweep, eijip,&
!lin_varsigmadown, lin_eprime, lin_Tprimepinv, lin_Fsweep)
!    implicit none
!    real, dimension(:,:,:,:), allocatable, intent(inout) :: varsigmadown
!    integer, dimension(:,:), intent(in) :: eprime
!    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Tprimepinv
!    integer, dimension(:,:), intent(in) :: sweeplist
!    real, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: Fsweep
!
!    real, dimension(:,:), allocatable, intent(inout) :: lin_varsigmadown
!    integer, dimension(:,:,:), allocatable, intent(inout) :: eijip
!    integer, dimension(:,:), allocatable, intent(inout) :: lin_eprime
!    real, dimension(:,:,:,:), allocatable, intent(inout) :: lin_Tprimepinv
!    real, dimension(:,:,:,:), allocatable, intent(inout) :: lin_Fsweep
!
!    allocate(eijip(NE,Nmu,Nnu))
!    allocate(lin_varsigmadown(NE*Nmu*Nnu,NFe))
!    allocate(lin_eprime(NE*Nmu*Nnu,NFe))
!    allocate(lin_Tprimepinv(NE*Nmu*Nnu,Gp,NKe,NKe))
!    allocate(lin_Fsweep(NE*Nmu*Nnu,NKe,NKe,NFe))
!
!    do ip = 1, NE*Nmu*Nnu
!        eijip(sweeplist(ip,1),sweeplist(ip,2),sweeplist(ip,3)) = ip
!    end do
!
!    do f = 1, NFe
!        do ip = 1, NE*Nmu*Nnu
!            if (eprime(sweeplist(ip,1),f) == 0) then
!                lin_eprime(ip,f) = 0
!            else
!                lin_eprime(ip,f) = eijip(eprime(sweeplist(ip,1),f),sweeplist(ip,2),sweeplist(ip,3))
!            end if
!        end do
!    end do
!
!    do ip = 1, NE*Nmu*Nnu
!        lin_varsigmadown(ip,:) = varsigmadown(sweeplist(ip,1),sweeplist(ip,2),sweeplist(ip,3),:)
!        lin_Tprimepinv(ip,:,:,:) = Tprimepinv(sweeplist(ip,1),:,sweeplist(ip,2),sweeplist(ip,3),:,:)
!        lin_Fsweep(ip,:,:,:) = Fsweep(sweeplist(ip,1),sweeplist(ip,2),sweeplist(ip,3),:,:,:)
!    end do
!
!    deallocate(varsigmadown)
!    deallocate(Tprimepinv)
!    deallocate(Fsweep)
!end subroutine lin_restructure_mesh

end module sweep_order