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

module energy_discretization
    use math
    use physics
    use user_input
implicit none

contains
subroutine construct_exceptions(particle, exc)
    implicit none
    character(*), intent(in) :: particle

    real, dimension(:), allocatable, intent(inout) :: exc

    integer :: ele, mat
    integer :: shell
    integer :: Nexc
    real, dimension(:), allocatable :: B_exc
    real, dimension(:), allocatable :: peWvals
    real, dimension(:), allocatable :: kept_peWvals
    real, dimension(:), allocatable :: trash_Espec
    real, dimension(:,:), allocatable :: trash_pexs
    real, dimension(:), allocatable :: trash_degen
    real, dimension(:), allocatable :: Wvals
    real, dimension(:), allocatable :: kept_Wvals
    real, dimension(:), allocatable :: trash_kinetic
    integer :: l_n

    ! E0 pm DeltaE CAN'T BE CONSIDERED SIMPLY AN EXCEPTION, BECAUSE I NEED FOR ONE GROUP TO BE IN THE RANGE E0 - DeltaE TO E0 + DeltaE.

    if (particle .eq. "photon") then
        if (all(user_photon_exc .eq. 0.0) .and. Epmax .ge. 1.022) then ! User did not provide exceptions
            Nexc = 1
            allocate(exc(1))
            exc(1) = 1.022
        else if (Epmax .ge. 1.022) then ! User did provide exceptions
            Nexc = N_user_photon_exc + 1
            allocate(exc, source = user_photon_exc)
            allocate(B_exc(Nexc))
            B_exc(1:Nexc-1) = exc
            B_exc(Nexc) = 1.022
            call move_alloc(B_exc, exc)
        else
            Nexc = N_user_photon_exc
            allocate(exc, source = user_photon_exc)
        end if
        if (photon_exc_shells) then ! User wants photoionization shell energies
            do mat = 1, Nmats
                do ele = 1, Nmaxeles
                    if (eles_list(mat,ele) .eq. 0) cycle
                    call read_photoelectric_effect(eles_list(mat,ele), Epmin, Epmax, peWvals, trash_Espec, trash_pexs)
                    allocate(kept_peWvals, &
                        source = pack(peWvals, mask = peWvals .ge. Epmin .and. peWvals .le. Epmax))

                    Nexc = Nexc + size(kept_peWvals)
                    allocate(B_exc(Nexc))
                    if (allocated(exc)) then !! DOUBLE CHECK THIS, IT'S RECURSIVE!!
                        B_exc(1:size(exc)) = exc !! DOUBLE CHECK THIS, IT'S RECURSIVE!!
                        B_exc(size(exc)+1:Nexc) = kept_peWvals !! DOUBLE CHECK THIS, IT'S RECURSIVE!!
                    else !! DOUBLE CHECK THIS, IT'S RECURSIVE!!
                        B_exc = kept_peWvals !! DOUBLE CHECK THIS, IT'S RECURSIVE!!
                    end if !! DOUBLE CHECK THIS, IT'S RECURSIVE!!
                    call move_alloc(B_exc, exc)

                    deallocate(peWvals)
                    deallocate(kept_peWvals)
                    deallocate(trash_Espec)
                    deallocate(trash_pexs)
                end do
            end do
        end if
    else if (particle .eq. "electron") then
        if (all(user_electron_exc .eq. 0.0)) then ! User did not provide exceptions
            if (electron_exc_shells) then ! User wants impact ionization shell energies, so empty array is created
                Nexc = 0
                allocate(exc(0))
            else ! User does not want impact ionization shell energies, so no exceptions exist
                allocate(exc(1))
                exc = 0.0
            end if
        else
            Nexc = N_user_electron_exc
            allocate(exc, source = user_electron_exc)
        end if
        if (electron_exc_shells) then
            do mat = 1, Nmats
                do ele = 1, Nmaxeles
                    if (eles_list(mat,ele) .eq. 0) cycle
                    call read_atomic_parameters(eles_list(mat,ele), trash_degen, Wvals, trash_kinetic)
                    allocate(kept_Wvals, &
                        source = pack(Wvals, mask = Wvals .ge. Eemin .and. Wvals .le. Eemax))

                    Nexc = Nexc + size(kept_Wvals)
                    allocate(B_exc(Nexc))
                    B_exc(1:size(exc)) = exc
                    B_exc(size(exc)+1:Nexc) = kept_Wvals
                    call move_alloc(B_exc, exc)

                    deallocate(Wvals)
                    deallocate(kept_Wvals)
                    deallocate(trash_degen)
                    deallocate(trash_kinetic)
                end do
            end do
        end if
    else
        if (beam_energy_dist .ne. "boxcar") then
            allocate(exc(1))
            exc = 0
        end if
    end if

    if (beam_energy_dist .eq. "boxcar" .and. &
        transport_mode .eq. "external photon beam" .and. &
        particle .eq. "photon") then ! User provided boxcar beam
        if (size(exc) .gt. 0) then
            allocate(B_exc(size(exc)+1))
            B_exc(1:size(exc)) = exc
            B_exc(size(exc)+1) = Epmax - deltaE
            call move_alloc(B_exc, exc)
        else
            allocate(B_exc(1))
            B_exc(1) = Epmax - deltaE
            call move_alloc(B_exc, exc)
        end if
    end if

    if (beam_energy_dist .eq. "boxcar" .and. &
        transport_mode .eq. "external electron beam" .and. &
        particle .eq. "electron") then ! User provided boxcar beam
        if (size(exc) .gt. 0) then
            allocate(B_exc(size(exc)+1))
            B_exc(1:size(exc)) = exc
            B_exc(size(exc)+1) = Eemax - deltaE
            call move_alloc(B_exc, exc)
        else
            allocate(B_exc(1))
            B_exc(1) = Eemax - deltaE
            call move_alloc(B_exc, exc)
        end if
    end if

    if (size(exc) .gt. 1 .and. exc(1) .eq. 0) then
        allocate(B_exc(size(exc)-1))
        B_exc(1:size(exc)-1) = exc(2:size(exc))
        call move_alloc(B_exc, exc)
    end if
end subroutine construct_exceptions

subroutine linear_energy_groups(l_G, Emin, Emax, exceptions, E, dE, Emid)
    implicit none
    integer, intent(in) :: l_G
    real, intent(in) :: Emin
    real, intent(in) :: Emax
    real, dimension(:), intent(in) :: exceptions

    real, dimension(:), allocatable, intent(inout) :: E
    real, dimension(:), allocatable, intent(inout) :: dE
    real, dimension(:), allocatable, intent(inout) :: Emid

    integer :: g, i
    real, dimension(:), allocatable :: l_exceptions
    real, dimension(:), allocatable :: B_l_exceptions
    integer :: hypothetical
    real, dimension(:), allocatable :: B_E
    integer :: counter

    allocate(E(l_G+1))
    allocate(dE(l_G))
    allocate(Emid(l_G))
    allocate(l_exceptions, source = exceptions)

    if (size(exceptions) .gt. 0) call qsort(l_exceptions)

    if (l_exceptions(1) .eq. 0) then
        do g = 1, l_G+1
            E(g) = -((Emax-Emin)/(l_G))*g + &
            ((1+1.0/(l_G))*Emax-(1.0/(l_G))*Emin)
        end do
    else
        allocate(B_E(l_G+1))

        hypothetical = l_G+1-size(l_exceptions)
        do g = 1, hypothetical
            B_E(g) = -((Emax-Emin)/(hypothetical-1))*g + &
            ((1+1.0/(hypothetical-1))*Emax-(1.0/(hypothetical-1))*Emin)
        end do

        do i = 1, size(l_exceptions)
            B_E(hypothetical+i) = l_exceptions(i)
        end do

        call qsort(B_E)

        E = B_E(l_G+1:1:-1)
    end if

    do g = 1, l_G
        dE(g) = E(g) - E(g+1)
        Emid(g) = 0.5*(E(g) + E(g+1))
    end do

    if (beam_energy_dist .eq. "boxcar") then
        if (l_G .eq. Gp .and. &
        transport_mode .eq. "external photon beam") then
            boxg = findloc(E, Epmax-deltaE, dim = 1) - 1
        else if (l_G .eq. Ge .and. &
        transport_mode .eq. "external electron beam") then
            boxg = findloc(E, Eemax-deltaE, dim = 1) - 1
        end if
    end if
end subroutine linear_energy_groups

subroutine logarithmic_energy_groups(l_G, Emin, Emax, exceptions, E, dE, Emid)
    implicit none
    integer, intent(in) :: l_G
    real, intent(in) :: Emin
    real, intent(in) :: Emax
    real, dimension(:), intent(in) :: exceptions

    real, dimension(:), allocatable, intent(inout) :: E
    real, dimension(:), allocatable, intent(inout) :: dE
    real, dimension(:), allocatable, intent(inout) :: Emid

    integer :: g, i
    real, dimension(:), allocatable :: l_exceptions
    real, dimension(:), allocatable :: B_l_exceptions
    integer :: hypothetical
    real, dimension(:), allocatable :: B_E

    allocate(E(l_G+1))
    allocate(dE(l_G))
    allocate(Emid(l_G))
    allocate(l_exceptions, source = exceptions)

    if (Emin .eq. 0.0) then
        print *, "STOP!"
        print *, "MODULE 7, SUBROUTINE logarithmic_energy_groups:"
        print *, "Minimum energy is provided as 0. This can not be done with logarithmic energy groups."
        print *, "Use linear energy groups."
        print *, "PROGRAM ENDING."
        stop
    end if

    if (size(exceptions) > 0) call qsort(l_exceptions)

    if (l_exceptions(1) .eq. 0) then
        do g = 1, l_G+1
            E(g) = exp(((log(Emin)-log(Emax))/l_G)*g + (1+1.0/l_G)*log(Emax)-(1.0/l_G)*log(Emin))
        end do
    else
        if (any(l_exceptions .lt. Emin) .or. any(l_exceptions .gt. Emax)) then
            print *, "WARNING!"
            print *, "MODULE 7, SUBROUTINE logarithmic_energy_groups:"
            print *, "Given exception is outside of energy range. The exception is:"
            print *, l_exceptions(&
                findloc(l_exceptions .lt. Emin .or. l_exceptions .gt. Emax, .true., dim = 1))
            print *, "This exception will be removed. Next time you run the program, you should remove it."
            print *, "Calculation will proceed."
            allocate(B_l_exceptions, source = pack(l_exceptions, mask = l_exceptions .ne. &
                l_exceptions(&
                findloc(l_exceptions .lt. Emin .or. l_exceptions .gt. Emax, .true., dim = 1))))
            call move_alloc(B_l_exceptions, l_exceptions)
        end if

        allocate(B_E(l_G+1))

        hypothetical = l_G+1-size(l_exceptions)
        do g = 1, hypothetical
            B_E(g) = exp(((log(Emin)-log(Emax))/(hypothetical-1))*g + &
            ((1+1.0/(hypothetical-1))*log(Emax)-(1.0/(hypothetical-1))*log(Emin)))
        end do

        do i = 1, size(l_exceptions)
            B_E(hypothetical+i) = l_exceptions(i)
        end do

        call qsort(B_E)

        E = B_E(l_G+1:1:-1)

    end if

    do g = 1, l_G
        dE(g) = E(g) - E(g+1)
        Emid(g) = 0.5*(E(g) + E(g+1))
    end do

    if (beam_energy_dist .eq. "boxcar") then
        if (l_G .eq. Gp .and. &
        transport_mode .eq. "external photon beam") then
            boxg = findloc(E, Epmax-deltaE, dim = 1) - 1
        else if (l_G .eq. Ge .and. &
        transport_mode .eq. "external electron beam") then
            boxg = findloc(E, Eemax-deltaE, dim = 1) - 1
        end if
    end if
end subroutine logarithmic_energy_groups

subroutine exponential_energy_groups(l_G, Emin, Emax, exceptions, E, dE, Emid)
    implicit none
    integer, intent(in) :: l_G
    real, intent(in) :: Emin
    real, intent(in) :: Emax
    real, dimension(:), intent(in) :: exceptions

    real, dimension(:), allocatable, intent(inout) :: E
    real, dimension(:), allocatable, intent(inout) :: dE
    real, dimension(:), allocatable, intent(inout) :: Emid

    integer :: g, i
    real, dimension(:), allocatable :: l_exceptions
    real, dimension(:), allocatable :: B_l_exceptions
    integer :: hypothetical
    real, dimension(:), allocatable :: B_E

    allocate(E(l_G+1))
    allocate(dE(l_G))
    allocate(Emid(l_G))
    allocate(l_exceptions, source = exceptions)

    if (size(exceptions) > 0) call qsort(l_exceptions)

    if (l_exceptions(1) .eq. 0) then
        do g = 1, l_G+1
            E(g) = log(((exp(Emin)-exp(Emax))/l_G)*g + (1+1.0/l_G)*exp(Emax)-(1.0/l_G)*exp(Emin))
        end do
    else
        if (any(l_exceptions .lt. Emin) .or. any(l_exceptions .gt. Emax)) then
            print *, "WARNING!"
            print *, "MODULE 7, SUBROUTINE exponential_energy_groups:"
            print *, "Given exception is outside of energy range. The exception is:"
            print *, l_exceptions(&
                findloc(l_exceptions .lt. Emin .or. l_exceptions .gt. Emax, .true., dim = 1))
            print *, "This exception will be removed. Next time you run the program, you should remove it."
            print *, "Calculation will proceed."
            allocate(B_l_exceptions, source = pack(l_exceptions, mask = l_exceptions .ne. &
                l_exceptions(&
                findloc(l_exceptions .lt. Emin .or. l_exceptions .gt. Emax, .true., dim = 1))))
            call move_alloc(B_l_exceptions, l_exceptions)
        end if

        allocate(B_E(l_G+1))

        hypothetical = l_G+1-size(l_exceptions)
        do g = 1, hypothetical
            B_E(g) = log(((exp(Emin)-exp(Emax))/(hypothetical-1))*g + &
            ((1+1.0/(hypothetical-1))*exp(Emax)-(1.0/(hypothetical-1))*exp(Emin)))
        end do

        do i = 1, size(l_exceptions)
            B_E(hypothetical+i) = l_exceptions(i)
        end do

        call qsort(B_E)

        E = B_E(l_G+1:1:-1)
    end if

    do g = 1, l_G
        dE(g) = E(g) - E(g+1)
        Emid(g) = 0.5*(E(g) + E(g+1))
    end do

    if (beam_energy_dist .eq. "boxcar") then
        if (l_G .eq. Gp .and. &
        transport_mode .eq. "external photon beam") then
            boxg = findloc(E, Epmax-deltaE, dim = 1) - 1
        else if (l_G .eq. Ge .and. &
        transport_mode .eq. "external electron beam") then
            boxg = findloc(E, Eemax-deltaE, dim = 1) - 1
        end if
    end if
end subroutine exponential_energy_groups

subroutine uncollided_beam_energy_distribution(fEspec, fE)
    implicit none
    real, dimension(:), allocatable, intent(inout) :: fEspec ! x values (energies)
    real, dimension(:), allocatable, intent(inout) :: fE ! y values (distribution at energy)

    integer :: Nspec
    real, dimension(2) :: reader2
    integer :: eof
    real, dimension(:), allocatable :: B_fEspec
    real, dimension(:), allocatable :: B_fE
    character(70) :: file_name
    real :: normalize

    allocate(fEspec(0))
    allocate(fE(0))

    file_name = "Linacs/"//trim(adjustl(linac))//".txt"


    open(1, file = trim(file_name))
    do
        read(1,*, iostat = eof) reader2

        if (eof < 0) exit

        allocate(B_fEspec(size(fEspec)+1))
        B_fEspec(1:size(fEspec)) = fEspec
        B_fEspec(size(fEspec)+1) = reader2(1)

        allocate(B_fE(size(fE)+1))
        B_fE(1:size(fE)) = fE
        B_fE(size(fE)+1) = reader2(2)

        call move_alloc(B_fEspec, fEspec)
        call move_alloc(B_fE, fE)
    end do
    close(1)

    call trapezoidal_integration_1D &
        (0, fE, fEspec, minval(fEspec), maxval(fEspec), normalize)

    fE = fE/normalize
end subroutine uncollided_beam_energy_distribution

subroutine verify_beam_spectrum(E, fEspec, fE)
    implicit none
    real, dimension(:), intent(in) :: E
    real, dimension(:), intent(in) :: fEspec
    real, dimension(:), intent(in) :: fE

    integer :: g, n
    integer :: l_G
    real, dimension(:), allocatable :: fEg


    print *, "---", "VERIFICATION OF BEAM: USER-DEFINED SPECTRUM ENERGY VALUES PRINTOUT STARTING:"
    do g = 1, size(fEspec)
        print *, fEspec(g)
    end do
    print *, "---", "VERIFICATION OF BEAM: USER-DEFINED SPECTRUM ENERGY VALUES PRINTOUT FINISHED."

    print *, "---", "VERIFICATION OF BEAM: USER-DEFINED SPECTRUM PRINTOUT STARTING:"
    do g = 1, size(fEspec)
        print *, fE(g)
    end do
    print *, "---", "VERIFICATION OF BEAM: USER-DEFINED SPECTRUM PRINTOUT FINISHED."

    l_G = size(E) - 1
    if (energy_discretization_method .eq. "FEXS") then
        allocate(fEg(l_G+1))
        do g = 1, l_G+1
            fEg(g) = interp1d(fEspec, fE, E(g))
        end do

        print *, "---", "VERIFICATION OF BEAM: ENERGY NODE PRINTOUT STARTING:"
        do g = 1, l_G+1
            print *, E(g)
        end do
        print *, "---", "VERIFICATION OF BEAM: ENERGY NODE PRINTOUT FINISHED."
        print *, "---", "VERIFICATION OF BEAM: ENERGY NODE PRINTOUT STARTING:"
        do g = 1, l_G+1
            print *, fEg(g)
        end do
        print *, "---", "VERIFICATION OF BEAM: ENERGY NODE PRINTOUT FINISHED."
    else if (energy_discretization_method .eq. "MGXS") then
        allocate(fEg(l_G))
        do g = 1, l_G
            call trapezoidal_integration_1D &
                (0, fE, fEspec, E(g+1), E(g), fEg(g))
        end do

        print *, "---", "VERIFICATION OF BEAM: ENERGY BIN ENDPOINT PRINTOUT STARTING:"
        do g = 1, l_G
            do n = 1, 2
                print *, E(g+n-1)
            end do
        end do
        print *, "---", "VERIFICATION OF BEAM: ENERGY BIN ENDPOINT PRINTOUT FINISHED."
        print *, "---", "VERIFICATION OF BEAM: ENERGY BIN DISTRIBUTION PRINTOUT STARTING:"
        do g = 1, l_G
            do n = 1, 2
                print *, fEg(g)/(E(g)-E(g+1))
            end do
        end do
        print *, "---", "VERIFICATION OF BEAM: ENERGY BIN DISTRIBUTION PRINTOUT FINISHED."
    end if

    print *, "PROGRAM ENDING."
    stop
end subroutine verify_beam_spectrum

end module energy_discretization