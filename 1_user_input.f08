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

module user_input
    use math
    use physics
implicit none

! UNDER DEVELOPMENT ! Do not change.
character(*), parameter :: mesh_reader = "GMSH"
logical, parameter :: verify_sweep = .false.
logical, parameter :: calculate_dose = .true.
logical, parameter :: estimate_storage = .true.
real, parameter :: max_storage = 16.0 ! GB
logical, parameter :: stolim = .false.
logical, parameter :: uncollided_fluence_correction = .false.
character(*), parameter :: FEM_type = "default"
logical, parameter :: results_as_dat = .true.
integer, parameter :: N_user_photon_exc = 1
real, dimension(N_user_photon_exc) :: user_photon_exc = [0.0]
!integer, parameter :: N_user_electron_exc = 3
!real, dimension(N_user_electron_exc) :: user_electron_exc = [5.7, 6.1, 6.5]
integer, parameter :: N_user_electron_exc = 1
real, dimension(N_user_electron_exc) :: user_electron_exc = [0.0]
integer :: NKe = 4
integer :: NKf = 3
integer :: NFe = 4
logical, parameter :: wiscoslab = .false.
logical :: MMS_bdy = .false.

! ### Solution parameters ###
character(50) :: transport_mode = "NULL"
logical :: FCS = .true.
real :: convergence = -1.0
character(50) :: photon_solution_method = "NULL"
integer :: p_pmax = -1
integer :: p_mmax = -1
character(50) :: electron_solution_method = "NULL"
integer :: e_pmax = -1
integer :: e_mmax = -1
character(50) :: MMS_solution_method = "NULL"
integer :: m_pmax = -1
integer :: m_mmax = -1
logical :: beams_defined_in_mesh = .true.
logical :: solve_photons = .false.
logical :: solve_electrons = .false.

! ### Mesh ###
character(100) :: selected_mesh = "NULL"
integer :: NK = -1
integer :: NE = -1

! ### Angular discretization ###
character(10) :: photon_angular = "NULL"
character(10) :: electron_angular = "NULL"
character(10) :: MMS_angular = "NULL"
integer :: Nmu = -1
integer :: Nphi = -1
integer :: NL = -1

! ### Energy discretization ###
character(10) :: energy_discretization_method = "NULL"
real :: Epmin = -1.0
real :: Epmax = -1.0
integer :: Gp = -1
character(20) :: photon_energy_structure = "NULL"
logical :: photon_exc_shells = .true.
real :: Eemin = -1.0
real :: Eemax = -1.0
integer :: Ge = -1
character(20) :: electron_energy_structure = "NULL"
logical :: electron_exc_shells = .false.
real :: Emmin = -1.0
real :: Emmax = -1.0

! ### Physics ###
character(20) :: electron_inelastic = "NULL"
character(20) :: electron_RCSDA = "NULL"
character(20) :: electron_elastic = "NULL"
real :: ccut = -1.0
logical :: elastic_etc = .true.
logical :: inelastic_etc = .true.
logical :: exact_RCSDA_angular = .true.
real :: MMS_attn = 1.0
logical :: customXS = .false.

! ### Materials ###
! # Material atomic numbers
integer :: Nmats
integer :: Nmaxeles
integer, dimension(:,:), allocatable :: eles_list
integer, dimension(:,:), allocatable :: Natoms
real, dimension(:), allocatable :: mat_rho

! ### Beam definition ###
character(20) :: beam_energy_dist = "NULL"
character(100) :: linac = "NULL"
real :: deltaE = -1.0
real, dimension(3) :: beam_origin
real, dimension(3) :: beam_axis
logical :: wide_open_collimator = .false.
character(50) :: beam_angular_dist = "NULL"
character(50) :: beam_cutout = "NULL"
logical :: beam_isotropic = .true.
real, dimension(4) :: beam_cutout_params = [0.0,0.0,0.0,0.0]

! ### Extra options ###
logical :: gen_post_proc = .false.
logical :: doslices = .false.
integer :: mshplane = -1
real :: mshslice = -1.0
logical :: doEDEP = .false.
logical :: doCDEP = .false.
logical :: doDDEP = .false.
logical :: doline = .false.
logical :: doplanarints = .false.
integer :: planarNP = -1
real :: planarx0 = -1.0
real :: planarxN = -1.0
logical :: pp_PN_to_SN = .false.
logical :: fresh_photon_run = .true.
logical :: fresh_electron_run = .true.
logical :: fresh_sweep_run = .true.
logical :: save_sweep = .false.
logical :: check_residuals = .false.
real, dimension(3) :: translate_vector = [0.0, 0.0, 0.0]
real, dimension(3) :: scale_factor = [1.0, 1.0, 1.0]
logical :: verify_beam_spec = .false.

! Public variables
logical :: isMGXS = .false.
logical :: isFEXS = .false.
logical :: isSN = .false.
logical :: isPN = .false.
logical :: t_isSN = .false.
logical :: t_isPN = .false.
logical :: isSI = .false.
logical :: isGMRESm = .false.
character(50) :: output_fname = "NULL"
integer, dimension(:), allocatable :: eltomat
integer, dimension(:,:), allocatable :: nodesinmat
integer, dimension(:), allocatable :: NKmat
integer, dimension(:,:), allocatable :: elsinmat
integer, dimension(:), allocatable :: NEmat
integer, dimension(:), allocatable :: nodesinbeam
integer, dimension(:), allocatable :: elsinbeam
integer, dimension(:,:), allocatable :: Nelnodes
integer, dimension(:), allocatable :: Nglobal
real, dimension(:), allocatable :: total_GAM
real, dimension(:,:), allocatable :: mass_frac
real, dimension(:), allocatable :: mat_rho_e
real, dimension(:,:), allocatable :: mat_rho_a
real, dimension(99) :: GAM
real, dimension(3,3) :: beam_coord_sys
real, dimension(3) :: k0
real :: aperture
real :: bt0
real :: bp0
real :: bd0
integer :: boxg

contains
subroutine read_input_file
    implicit none
    integer :: mat, i
    integer, dimension(99) :: anum
    integer, dimension(99) :: Nelespermat
    integer :: io_status
    integer, dimension(:,:), allocatable :: B_eles_list
    character(50) :: trashchar

    open(1, file = "input.in")

    ! ### Solution parameters ###
    read(1,'(A)') trashchar

    ! # Problem
    read(1,'(A)') trashchar
    read(1,'(A)') transport_mode
    ! CHECK: # Problem
    if (transport_mode .ne. "external photon beam" .and. &
        transport_mode .ne. "external photon beam coupled" .and. &
        transport_mode .ne. "external electron beam" .and. &
        transport_mode .ne. "T MMS" .and. &
        transport_mode .ne. "full MMS") then
        print *, "STOP!"
        print *, "Input file: You have entered an invalid input for # Problem:"
        print *, transport_mode
        print *, "Please choose either: external photon beam, external photon beam coupled"
        print *, "or external electron beam."
        print *, "PROGRAM ENDING."
        stop
    end if
    if (transport_mode .eq. "external photon beam") then
        solve_photons = .true.
    else if (transport_mode .eq. "external photon beam coupled") then
        solve_photons = .true.
        solve_electrons = .true.
    else if (transport_mode .eq. "external electron beam") then
        solve_electrons = .true.
    end if
    read(1,'(A)') trashchar ! BLANK

    ! # Use FCS
    read(1,'(A)') trashchar
    read(1,'(A)') trashchar
    if (trashchar .eq. "false") FCS = .false.
    read(1,'(A)') trashchar ! BLANK

    ! #* [particle] solution method
    read(1,'(A)') trashchar
    if (trashchar .eq. "### Mesh ###") then
        print *, "STOP!"
        print *, "Input file: Your input file does not specify particle types nor their solution methods."
        print *, "Please provide this information in the proper format."
        print *, "PROGRAM ENDING."
        stop
    end if
    1 continue
    if (trashchar .eq. "#* Photon solution method") then
        ! #* Photon solution method
        read(1,'(A)') photon_solution_method
        ! #--> Max number of iterates
        if (photon_solution_method .eq. "SI") then
            read(1,'(A)') trashchar
            read(1,*) p_pmax
        else if (photon_solution_method .eq. "GMRESm") then
            read(1,'(A)') trashchar
            read(1,*) p_pmax, p_mmax
        else if (solve_photons) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid photon solution method:"
            print *, photon_solution_method
            print *, "Please choose either SI or GMRESm."
            print *, "PROGRAM ENDING."
            stop
        else
            read(1,'(A)') trashchar
            read(1,'(A)') trashchar
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 1
    else if (trashchar .eq. "#* Electron solution method") then
        ! #* Electron solution method
        read(1,'(A)') electron_solution_method
        ! #--> Max number of iterates
        if (electron_solution_method .eq. "SI") then
            read(1,'(A)') trashchar
            read(1,*) e_pmax
        else if (electron_solution_method .eq. "GMRESm") then
            read(1,'(A)') trashchar
            read(1,*) e_pmax, e_mmax
        else if (solve_electrons) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid electron solution method:"
            print *, Electron_solution_method
            print *, "Please choose either SI or GMRESm."
            print *, "PROGRAM ENDING."
            stop
        else
            read(1,'(A)') trashchar
            read(1,'(A)') trashchar
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 1
    else if (trashchar .eq. "#* MMS solution method") then
        ! #* MMS solution method
        read(1,'(A)') MMS_solution_method
        ! #--> Max number of iterates
        if (MMS_solution_method .eq. "SI") then
            read(1,'(A)') trashchar
            read(1,*) m_pmax
        else if (MMS_solution_method .eq. "GMRESm") then
            read(1,'(A)') trashchar
            read(1,*) m_pmax, m_mmax
        else if (transport_mode .ne. "T MMS" .or. transport_mode .ne. "full MMS") then ! Must update when more MMS modes are implemented
            print *, "STOP!"
            print *, "Input file: You have entered an invalid MMS solution method:"
            print *, MMS_solution_method
            print *, "Please choose either SI or GMRESm."
            print *, "PROGRAM ENDING."
            stop
        else
            read(1,'(A)') trashchar
            read(1,'(A)') trashchar
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 1
    end if

    ! CHECK: Verify that problem has all of the necessary solution methods
    if (transport_mode .eq. "external photon beam") then
        if (photon_solution_method .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified external photon beam as your problem,"
            print *, "but you have not specified the photon solution method."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    else if (transport_mode .eq. "external photon beam coupled") then
        if (photon_solution_method .eq. "NULL" .or. &
            electron_solution_method .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified external photon beam coupled as your problem,"
            print *, "but you have not specified either the photon solution method,"
            print *, "the electron solution method, or both."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    else if (transport_mode .eq. "external electron beam") then
        if (electron_solution_method .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified external electron beam as your problem,"
            print *, "but you have not specified the electron solution method."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    else if (transport_mode .eq. "MMS") then
        if (MMS_solution_method .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified MMS as your problem,"
            print *, "but you have not specified the MMS solution method."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    end if

    ! CHECK: Verify that the max number of iterates are properly entered
    if (solve_photons .and. photon_solution_method .eq. "SI" .and. p_pmax .le. 0) then
        print *, "STOP!"
        print *, "Input file: You have entered an invalid number of iterates for"
        print *, "photon source iteration, or you have not entered a number at all."
        print *, "At the moment, you would be using:"
        print *, p_pmax
        print *, "Please enter a value greater than zero."
        print *, "PROGRAM ENDING."
        stop
    else if (solve_photons .and. photon_solution_method .eq. "GMRESm" .and. &
        (p_pmax .le. -1 .or. p_mmax .le. -1)) then
        print *, "STOP!"
        print *, "Input file: You have entered an invalid number of iterates for"
        print *, "photon GMRESm, or you have not entered a number at all."
        print *, "At the moment, you would be using:"
        print *, p_pmax, p_mmax
        print *, "Please enter two values greater than zero."
        print *, "PROGRAM ENDING."
        stop
    end if
    if (solve_electrons .and. electron_solution_method .eq. "SI" .and. e_pmax .eq. -1) then
        print *, "STOP!"
        print *, "Input file: You have entered an invalid number of iterates for"
        print *, "electron source iteration, or you have not entered a number at all."
        print *, "At the moment, you would be using:"
        print *, e_pmax
        print *, "Please enter a value greater than zero."
        print *, "PROGRAM ENDING."
        stop
    else if (solve_electrons .and. electron_solution_method .eq. "GMRESm" .and. &
        (e_pmax .le. -1 .or. e_mmax .le. -1)) then
        print *, "STOP!"
        print *, "Input file: You have entered an invalid number of iterates for"
        print *, "electron GMRESm, or you have not entered a number at all."
        print *, "At the moment, you would be using:"
        print *, e_pmax, e_mmax
        print *, "Please enter two values greater than zero."
        print *, "PROGRAM ENDING."
        stop
    end if
    if (MMS_solution_method .eq. "SI" .and. m_pmax .eq. -1) then
        print *, "STOP!"
        print *, "Input file: You have entered an invalid number of iterates for"
        print *, "MMS source iteration, or you have not entered a number at all."
        print *, "At the moment, you would be using:"
        print *, m_pmax
        print *, "Please enter a value greater than zero."
        print *, "PROGRAM ENDING."
        stop
    else if (MMS_solution_method .eq. "GMRESm" .and. &
        (m_pmax .le. -1 .or. m_mmax .le. -1)) then
        print *, "STOP!"
        print *, "Input file: You have entered an invalid number of iterates for"
        print *, "MMS GMRESm, or you have not entered a number at all."
        print *, "At the moment, you would be using:"
        print *, m_pmax, m_mmax
        print *, "Please enter two values greater than zero."
        print *, "PROGRAM ENDING."
        stop
    end if

    ! # Convergence criterion (already read)
    read(1,*) convergence

    ! CHECK: # Convergence criterion
    if (convergence .lt. 0) then
        print *, "STOP!"
        print *, "Input file: You have entered a negative value for # Convergence criterion:"
        print *, convergence
        print *, "Please enter a positive value."
        print *, "PROGRAM ENDING."
    end if
    read(1,'(A)') trashchar ! BLANK

    ! ### Mesh ###
    read(1,'(A)') trashchar
    ! # Selected mesh
    read(1,'(A)') trashchar
    read(1,'(A)') selected_mesh
    !selected_mesh = trim(adjustl(selected_mesh))
    read(1,'(A)') trashchar ! BLANK

    ! # Beam in mesh
    read(1,'(A)') trashchar
    read(1,'(A)') trashchar
    if (trashchar .eq. "true") then
        beams_defined_in_mesh = .true.
    else
        beams_defined_in_mesh = .false.
    end if
    read(1,'(A)') trashchar ! BLANK

    ! ### Angular discretization ###
    read(1,'(A)') trashchar

    ! #* [particle] angular treatment
    read(1,'(A)') trashchar
    if (trashchar .eq. "### Energy discretization ###") then
        print *, "STOP!"
        print *, "Input file: Your input file does not specify angular discretization methods"
        print *, "for the relevant particle types."
        print *, "Please provide this information in the proper format."
        print *, "PROGRAM ENDING."
        stop
    end if
    2 continue
    if (trashchar .eq. "#* Photon angular treatment") then
        read(1,'(A)') photon_angular
        if (solve_photons .and. photon_angular .ne. "SN" .and. &
            photon_angular .ne. "PN") then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid photon angular treatment:"
            print *, photon_angular
            print *, "Please choose either SN or PN."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 2
    else if (trashchar .eq. "#* Electron angular treatment") then
        read(1,'(A)') electron_angular
        if (solve_electrons .and. electron_angular .ne. "SN" .and. &
            electron_angular .ne. "PN") then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid electron angular treatment:"
            print *, electron_angular
            print *, "Please choose either SN or PN."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 2
    else if (trashchar .eq. "#* MMS angular treatment") then
        read(1,'(A)') MMS_angular
        if ((transport_mode .eq. "T MMS" .or. transport_mode .eq. "full MMS") .and. &
            MMS_angular .ne. "SN" .and. &
            MMS_angular .ne. "PN") then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid MMS angular treatment:"
            print *, MMS_angular
            print *, "Please choose either SN or PN."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 2
    end if

    ! CHECK: Verify that the user has specified all of the necessary angular treatments
    if (transport_mode .eq. "external photon beam") then
        if (photon_angular .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified external photon beam as your problem,"
            print *, "but you have not specified the photon angular treatment."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    else if (transport_mode .eq. "external photon beam coupled") then
        if (photon_angular .eq. "NULL" .or. &
            electron_angular .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified external photon beam as your problem,"
            print *, "but you have not specified either the photon angular treatment,"
            print *, "the electron angular treatment, or both."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    else if (transport_mode .eq. "external electron beam") then
        if (electron_angular .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified external electron beam as your problem,"
            print *, "but you have not specified the electron angular treatment."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    else if (transport_mode .eq. "MMS") then
        if (MMS_angular .eq. "NULL") then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified MMS as your problem,"
            print *, "but you have not specified the MMS angular treatment."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    end if

    ! # Angular discretization parameter (already read)
    read(1,*) Nmu
    Nphi = 2*Nmu
    NL = Nmu-1
    read(1,'(A)') trashchar ! BLANK

    ! ### Energy discretization ###
    read(1,'(A)') trashchar
    ! # Energy discretization method
    read(1,'(A)') trashchar
    read(1,'(A)') energy_discretization_method
    isFEXS = energy_discretization_method .eq. "FEXS"
    isMGXS = .not. isFEXS
    if (isFEXS .and. FCS) then
        print *, "WARNING!"
        print *, "Input file: You have turned on the FEXS energy discretization method,"
        print *, "but also indicated that you wish to use the FCS method. At the moment,"
        print *, "these are incompatible. The solution will proceed with FCS off."
        FCS = .false.
    end if

    read(1,'(A)') trashchar ! BLANK

    ! #* [particle] energy min and max
    read(1,'(A)') trashchar
    if (trashchar .eq. "### Physics ###") then
        print *, "STOP!"
        print *, "Input file: Your input file does not specify energy discretization methods"
        print *, "for the relevant particle types."
        print *, "Please provide this information in the proper format."
        print *, "PROGRAM ENDING."
        stop
    end if
    3 continue
    if (trashchar .eq. "#* Photon energy min and max") then
        if (solve_photons) then
            read(1,*) Epmin, Epmax
        else
            read(1,'(A)') trashchar
        end if
        if (solve_photons .and. Epmin .lt. 1.0E-5) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid photon minimum energy:"
            print *, Epmin
            print *, "Please enter no less than 1.0E-5."
            print *, "PROGRAM ENDING."
            stop
        end if
        if (solve_photons .and. Epmax .gt. 100.0) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid photon maximum energy:"
            print *, Epmax
            print *, "Please enter no larger than 100.0."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        ! #* Number of photon energy groups
        read(1,'(A)') trashchar
        if (trashchar .ne. "#* Number of photon energy groups") then
            print *, "STOP!"
            print *, "Input file: Your input file is improperly formatted."
            print *, "It's missing the #* Number of photon energy groups subheader."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
        end if
        if (solve_photons) then
            read(1,*) Gp
        else
            read(1,'(A)') trashchar
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! #* Photon energy structure
        if (trashchar .ne. "#* Photon energy structure") then
            print *, "STOP!"
            print *, "Input file: Your input file is improperly formatted."
            print *, "It's missing the #* Photon energy structure subheader."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
        end if
        read(1,'(A)') photon_energy_structure
        if (solve_photons .and. photon_energy_structure .ne. "logarithmic" .and. &
            photon_energy_structure .ne. "linear" .and. &
            photon_energy_structure .ne. "exponential") then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid photon energy structure:"
            print *, photon_energy_structure
            print *, "At the moment, your three choices are:"
            print *, "logarithmic"
            print *, "linear"
            print *, "exponential"
            print *, "Please enter one of these choices."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        ! #* Use photoelectric shell energies
        read(1,'(A)') trashchar
        if (trashchar .ne. "#* Use photoelectric shell energies") then
            print *, "STOP!"
            print *, "Input file: Your input file is improperly formatted."
            print *, "It's missing the #* Use photoelectric shell energies subheader."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
        end if
        read(1,'(A)') trashchar
        if (trashchar .eq. "true") then
            photon_exc_shells = .true.
        else if (trashchar .eq. "false") then
            photon_exc_shells = .false.
        else if (solve_photons) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid choice for"
            print *, "#* Use photoelectric shell energies. You have entered:"
            print *, trashchar
            print *, "Please enter either true or false."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 3
    else if (trashchar .eq. "#* Electron energy min and max") then
        if (solve_electrons) then
            read(1,*) Eemin, Eemax
        else
            read(1,'(A)') trashchar
        end if
        if (solve_electrons .and. Eemin .lt. 1.0E-5) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid electron minimum energy:"
            print *, Eemin
            print *, "Please enter no less than 1.0E-5."
            print *, "PROGRAM ENDING."
            stop
        end if
        if (solve_electrons .and. Eemax .gt. 100.0) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid electron maximum energy:"
            print *, Eemax
            print *, "Please enter no larger than 100.0."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        ! #* Number of electron energy groups
        read(1,'(A)') trashchar
        if (trashchar .ne. "#* Number of electron energy groups") then
            print *, "STOP!"
            print *, "Input file: Your input file is improperly formatted."
            print *, "It's missing the #* Number of electron energy groups subheader."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
        end if
        if (solve_electrons) then
            read(1,*) Ge
        else
            read(1,'(A)') trashchar
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! #* Electron energy structure
        if (trashchar .ne. "#* Electron energy structure") then
            print *, "STOP!"
            print *, "Input file: Your input file is improperly formatted."
            print *, "It's missing the #* Electron energy structure subheader."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
        end if
        read(1,'(A)') electron_energy_structure
        if (solve_electrons .and. electron_energy_structure .ne. "logarithmic" .and. &
            electron_energy_structure .ne. "linear" .and. &
            electron_energy_structure .ne. "exponential") then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid electron energy structure:"
            print *, electron_energy_structure
            print *, "At the moment, your three choices are:"
            print *, "logarithmic"
            print *, "linear"
            print *, "exponential"
            print *, "Please enter one of these choices."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        ! #* Use impact ionization shell energies
        read(1,'(A)') trashchar
        if (trashchar .ne. "#* Use impact ionization shell energies") then
            print *, "STOP!"
            print *, "Input file: Your input file is improperly formatted."
            print *, "It's missing the #* Use impact ionization shell energies subheader."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
        end if
        read(1,'(A)') trashchar
        if (trashchar .eq. "true") then
            electron_exc_shells = .true.
        else if (trashchar .eq. "false") then
            electron_exc_shells = .false.
        else if (solve_electrons) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid choice for"
            print *, "#* Use impact ionization shell energies. You have entered:"
            print *, trashchar
            print *, "Please enter either true or false."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 3
    end if

    ! ### Physics ### (already read)
    read(1,'(A)') trashchar
    if (trashchar .eq. "") go to 40
    4 continue
    if (trashchar .eq. "#* Electron inelastic scattering cross section") then
        read(1,'(A)') electron_inelastic
        if (solve_electrons .and. &
            electron_inelastic .ne. "Moller" .and. electron_inelastic .ne. "RBED") then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid choice for"
            print *, "#* Electron inelastic scattering cross section. You have entered:"
            print *, electron_inelastic
            print *, "Please enter either Moller or RBED."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    else if (trashchar .eq. "#* Electron RCSDA cross section (MGXS only)") then
        read(1,'(A)') electron_RCSDA
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    else if (trashchar .eq. "#* Electron elastic scattering cross section") then
        read(1,'(A)') electron_elastic
        if (solve_electrons .and. &
            electron_elastic .ne. "Wentzel-Moliere" .and. electron_elastic .ne. "fitted" .and. &
            electron_elastic .ne. "ELSEPA") then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid choice for"
            print *, "#* Electron elastic scattering cross section. You have entered:"
            print *, electron_elastic
            print *, "Please enter either Wentzel-Moliere, fitted, or ELSEPA."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    else if (trashchar .eq. "#* Electron RCSDA cutoff (FEXS only)") then
        if (solve_electrons .and. isFEXS) then
            read(1,*) ccut
        else
            ccut = 1.0
            read(1,'(A)') trashchar
        end if
        if (solve_electrons .and. ccut .le. 0) then
            print *, "STOP!"
            print *, "Input file: You have entered a negative value for"
            print *, "#* Electron RCSDA cutoff (FEXS only). You have entered:"
            print *, ccut
            print *, "Please enter a positive value."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    else if (trashchar .eq. "#* Extended transport correction (elastic)") then
        read(1,'(A)') trashchar
        if (trashchar .eq. "true") then
            elastic_etc = .true.
        else if (trashchar .eq. "false") then
            elastic_etc = .false.
        else if (solve_electrons) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid choice for"
            print *, "#* Extended transport correction (elastic). You have entered:"
            print *, trashchar
            print *, "Please enter either true or false."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    else if (trashchar .eq. "#* Extended transport correction (inelastic)") then
        read(1,'(A)') trashchar
        if (trashchar .eq. "true") then
            inelastic_etc = .true.
        else if (trashchar .eq. "false") then
            inelastic_etc = .false.
        else if (solve_electrons) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid choice for"
            print *, "#* Extended transport correction (inelastic). You have entered:"
            print *, trashchar
            print *, "Please enter either true or false."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    else if (trashchar .eq. "#* Exact RCSDA angular treatment") then
        read(1,'(A)') trashchar
        if (trashchar .eq. "true") then
            exact_RCSDA_angular = .true.
        else if (trashchar .eq. "false") then
            exact_RCSDA_angular = .false.
        else if (solve_electrons) then
            print *, "STOP!"
            print *, "Input file: You have entered an invalid choice for"
            print *, "#* Exact RCSDA angular treatment. You have entered:"
            print *, trashchar
            print *, "Please enter either true or false."
            print *, "PROGRAM ENDING."
            stop
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    else if (trashchar .eq. "#* MMS attenuation coefficient") then
        if (transport_mode .eq. "T MMS" .or. transport_mode .eq. "full MMS") then
            read(1,*) MMS_attn
        else
            read(1,'(A)') trashchar
        end if
        if (MMS_attn .lt. 0) then
            print *, "WARNING!"
            print *, "Input file: You have provided a negative MMS attenuation coefficient."
            print *, "Please note that this solver is not verified for negative attenuation coefficients."
            print *, "Take the results with a grain of salt."
        end if
        read(1,'(A)') trashchar ! BLANK
        read(1,'(A)') trashchar ! Next entry
        go to 4
    end if

    40 continue
    ! CHECK: Verify that the user has specified the required entries for their problem type.
    if (transport_mode .eq. "external electron beam" .or. &
        transport_mode .eq. "external photon beam coupled") then
        if (electron_inelastic .eq. "NULL" .or. &
            ccut .eq. -1.0) then
            print *, "STOP!"
            print *, "Input file: You are missing inputs."
            print *, "You have specified external electron beam or external photon beam coupled"
            print *, "as your problem, but you have not specified either the electron inelastic"
            print *, "mechanism or RCSDA cutoff energy."
            print *, "Please fix the input file."
            print *, "PROGRAM ENDING."
            stop
        end if
    end if

    ! ### Materials ### (already read)
    ! # Material atomic numbers
    read(1,'(A)') trashchar
    allocate(eles_list(0,0))
    mat = 0
    Nmats = 0
    Nmaxeles = 0
    do
        i = 0
        read(1,'(A)') trashchar
        do
            i = i + 1
            read(trashchar,*) anum(1:i)
            if (anum(i) .eq. -1 .and. i .eq. 1) exit
            if (anum(i) .eq. -1) exit

            Nmaxeles = max(i, Nmaxeles)
        end do
        if (anum(i) .eq. -1 .and. i .eq. 1) exit

        Nmats = Nmats + 1
        Nelespermat(Nmats) = i-1
        allocate(B_eles_list(Nmats,Nmaxeles))
        B_eles_list = 0
        do mat = 1, Nmats-1
            B_eles_list(mat,1:Nelespermat(mat)) = eles_list(mat,1:Nelespermat(mat))
        end do
        B_eles_list(Nmats,1:i-1) = anum(1:i-1)

        call move_alloc(B_eles_list, eles_list)
    end do
    read(1,'(A)') trashchar ! BLANK

    ! # Corresponding number of atoms
    read(1,'(A)') trashchar
    if (trashchar .ne. "# Corresponding number of atoms") then
        print *, "STOP!"
        print *, "Input file: Your input file is improperly formatted."
        print *, "The subheader # Material atomic numbers was potentially read or formatted incorrectly."
        print *, "This is what the program has read, by material:"
        do mat = 1, Nmats
            print *, eles_list(mat,1:Nelespermat(mat))
        end do
        print *, "If this is correct, ensure that you have created an extra entry with only -1 at"
        print *, "the end of your materials list. If this is not correct, ensure that you have"
        print *, "included -1 at the end of each entry of your materials list. See the example"
        print *, "input file for a properly formatted materials list."
        print *, "PROGRAM ENDING."
    end if
    allocate(Natoms(Nmats,Nmaxeles))
    Natoms = 0
    do mat = 1, Nmats
        read(1,*) Natoms(mat,1:Nelespermat(mat))
    end do
    read(1,'(A)') trashchar ! BLANK

    ! # Material densities
    read(1,'(A)') trashchar
    if (trashchar .ne. "# Material densities") then
        print *, "STOP!"
        print *, "Input file: Your input file is improperly formatted."
        print *, "The subheader # Corresponding number of atoms was potentially read or"
        print *, "formatted incorrectly. This is what the program has read, by material:"
        do mat = 1, Nmats
            print *, Natoms(mat,1:Nelespermat(mat))
        end do
        print *, "If this is correct, ensure that you have created an extra entry with only -1 at"
        print *, "the end of your number of atoms list. If this is not correct, ensure that"
        print *, "the number of rows you have entered under # Material atomic numbers"
        print *, "matches the number of rows you have entered under # Corresponding number of atoms."
        print *, "See the example input file for a properly formatted materials list."
        print *, "PROGRAM ENDING."
    end if
    allocate(mat_rho(Nmats))
    do mat = 1, Nmats
        read(1,*) mat_rho(mat)
    end do
    read(1,'(A)') trashchar ! BLANK

    ! CHECK: Verify that no elements are greater than 99
    ! NOTE: add custom materials past 99.
    if (any(eles_list .gt. 99)) then
        print *, "STOP!"
        print *, "You have entered a material with an element having atomic number"
        print *, "exceeding 99. Please use only elements up to 99."
        print *, "PROGRAM ENDING."
        stop
    end if

    ! ### Beam definition ###
    read(1,'(A)') trashchar
    ! # Beam energy distribution
    read(1,'(A)') trashchar
    read(1,'(A)') beam_energy_dist
    if (beam_energy_dist .ne. "polychromatic" .and. &
        beam_energy_dist .ne. "boxcar" .and. &
        beam_energy_dist .ne. "monochromatic") then
        print *, "STOP!"
        print *, "Input file: You have entered an beam energy distribution type:"
        print *, beam_energy_dist
        print *, "At the moment, your two choices are:"
        print *, "polychromatic"
        print *, "boxcar"
        print *, "monochromatic"
        print *, "Please enter one of these choices."
        print *, "PROGRAM ENDING."
        stop
    end if
    read(1,'(A)') trashchar
    5 continue
    if (trashchar .eq. "#--> Linac") then
        read(1,'(A)') linac
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 5
    else if (trashchar .eq. "#--> Energy width") then
        if (beam_energy_dist .eq. "boxcar") then
            read(1,*) deltaE
        else
            read(1,'(A)') trashchar
        end if
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 5
    end if

    ! Make changes to user input for monochromatic beam energy distribution
    if (beam_energy_dist .eq. "monochromatic") then
        if (transport_mode .eq. "external photon beam" .or. &
            transport_mode .eq. "external photon beam coupled") then
            Epmax = Epmax + (Epmax - Epmin)/(2*Gp-1)
        else if (transport_mode .eq. "external electron beam") then
            Eemax = Eemax + (Eemax - Eemin)/(2*Ge-1)
        end if
    end if

    ! # Beam origin
    read(1,'(A)') trashchar
    read(1,*) beam_origin
    read(1,'(A)') trashchar ! BLANK

    ! # Beam axis
    read(1,'(A)') trashchar
    read(1,*) beam_axis
    read(1,'(A)') trashchar ! BLANK

    ! # Beam angular distribution
    read(1,'(A)') trashchar
    read(1,'(A)') beam_angular_dist
    if (beam_angular_dist .ne. "spherical" .and. &
        beam_angular_dist .ne. "planar") then
        print *, "STOP!"
        print *, "Input file: You have entered an beam angular distribution type:"
        print *, beam_angular_dist
        print *, "At the moment, your two choices are:"
        print *, "spherical"
        print *, "planar"
        print *, "Please enter one of these choices."
        print *, "PROGRAM ENDING."
        stop
    end if
    read(1,'(A)') trashchar ! BLANK

    ! # Beam cutout
    read(1,'(A)') trashchar
    read(1,'(A)') beam_cutout
    if (beam_cutout .eq. "none") then
        wide_open_collimator = .true.
        read(1,'(A)') trashchar ! Next entry OR BLANK
        if (trashchar .eq. "#--> Beam cutout parameters") then
            read(1,'(A)') trashchar
        end if
    else if (beam_cutout .eq. "rectangle") then
        read(1,'(A)') trashchar
        beam_cutout_params = 0
        read(1,*) beam_cutout_params(1:4)
    else if (beam_cutout .eq. "circle") then
        read(1,'(A)') trashchar
        beam_cutout_params = 0
        read(1,*) beam_cutout_params(1:2)
    else
        print *, "STOP!"
        print *, "Input file: You have entered an invalid beam cutout:"
        print *, beam_cutout
        print *, "At the moment, your three choices are:"
        print *, "rectangle"
        print *, "circle"
        print *, "none"
        print *, "Please enter one of these choices."
        print *, "PROGRAM ENDING."
        stop
    end if
    read(1,'(A)') trashchar ! BLANK

    ! ### Extra options ###
    read(1,'(A)') trashchar
    read(1,'(A)') trashchar
    6 continue
    if (trashchar .eq. "post processing") then
        gen_post_proc = .true.
        read(1,*) mshplane, mshslice
        if (mshplane .gt. 0) doslices = .true.
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "calculate energy deposition") then
        doEDEP = .true.
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "calculate charge deposition") then
        doCDEP = .true.
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "calculate dose deposition") then
        doDDEP = .true.
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "calculate line distributions") then
        doline = .true.
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "planar integrals along depth") then
        doplanarints = .true.
        read(1,*) planarNP, planarx0, planarxN
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "PN to SN") then
        pp_PN_to_SN = .true.
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "use stored photon solution") then
        fresh_photon_run = .false.
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "use stored electron solution") then
        fresh_electron_run = .false.
        read(1,'(A)') trashchar
        go to 6
    else if (trashchar .eq. "use stored sweep") then
        fresh_sweep_run = .false.
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "store sweep") then
        save_sweep = .true.
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "check residuals") then
        check_residuals = .true.
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "scale mesh") then
        read(1,*) scale_factor
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "translate mesh") then
        read(1,*) translate_vector
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "verify beam spectrum") then
        verify_beam_spec = .true.
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "use custom physics") then
        customXS = .true.
        read(1,'(A)') trashchar ! Next entry OR BLANK
        go to 6
    else if (trashchar .eq. "") then
        go to 60
    else
        print *, "STOP!"
        print *, "Input file: You have entered an invalid extra option."
        print *, "You have entered:"
        print *, trashchar
        print *, "Please enter only valid options."
        print *, "PROGRAM ENDING."
        stop
    end if
    60 continue
    close(1)
end subroutine read_input_file

subroutine check_for_WIP
    implicit none

    if (photon_angular .eq. "PN" .or. &
        electron_angular .eq. "PN") then
        print *, "STOP!"
        print *, "You have entered PN for some particle type. At the moment,"
        print *, "this has not been validated in wiscobolt. Please use SN."
        print *, "PROGRAM ENDING."
        stop
    end if
    if (.not. any(eles_list .eq. 1) .and. &
        .not. any(eles_list .eq. 8) .and. &
        .not. any(eles_list .eq. 13) .and. &
        .not. any(eles_list .eq. 26) .and. &
        .not. any(eles_list .eq. 29) .and. &
        .not. any(eles_list .eq. 47) .and. &
        .not. any(eles_list .eq. 79) .and. &
        .not. any(eles_list .eq. 83)) then
        print *, "STOP!"
        print *, "You have entered an element other than (atomic number):"
        print *, 1, 8, 13, 26, 29, 47, 79, 83
        print *, "At the moment, these are the only elements implemented."
        print *, "PROGRAM ENDING."
        stop
    end if
    if (Nmats .gt. 1) then
        print *, "STOP!"
        print *, "You are using more than 1 material in the problem."
        print *, "At the moment, this is not ready for use, because of"
        print *, "significant errors in the solutions."
        print *, "PROGRAM ENDING."
        stop
    end if
end subroutine check_for_WIP

subroutine prepare_globals
    implicit none
    integer :: ele, i, mat
    real :: bx
    real :: by
    real :: bz
    character(50) :: test

    open(1, file = "Physics data/Atomic parameters/atomic_masses.txt")
    do i = 1, 99
        read(1,*) GAM(i)
    end do
    close(1)

    allocate(total_GAM(Nmats))
    allocate(mass_frac(Nmats,Nmaxeles))
    allocate(mat_rho_e(Nmats))
    allocate(mat_rho_a(Nmats,Nmaxeles))
    total_GAM = 0
    mat_rho_e = 0
    mat_rho_a = 0
    mass_frac = 0
    do mat = 1, Nmats
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            total_GAM(mat) = total_GAM(mat) + Natoms(mat,ele)*GAM(eles_list(mat,ele))
        end do
        do ele = 1, Nmaxeles
            if (eles_list(mat,ele) .eq. 0) cycle
            mat_rho_e(mat) = mat_rho_e(mat) + &
                Avogadro*mat_rho(mat)*Natoms(mat,ele)*eles_list(mat,ele)/total_GAM(mat)
            mat_rho_a(mat,ele) = Avogadro*mat_rho(mat)*Natoms(mat,ele)/total_GAM(mat)
            mass_frac(mat,ele) = Natoms(mat,ele)*GAM(eles_list(mat,ele))/total_GAM(mat)
        end do
    end do

    if (beam_angular_dist .eq. "planar") then
        if (beam_cutout .eq. "circle") then ! FINAL UNITS: s0^(-1)
            aperture = twopi*beam_cutout_params(2)**2
        else if (beam_cutout .eq. "rectangle") then ! FINAL UNITS: s0^(-1)
            aperture = beam_cutout_params(2)*beam_cutout_params(3)
        else if (beam_cutout .eq. "none") then ! FINAL UNITS: (s0/beam area)^(-1) (so, areal density of particles)
            aperture = 1
        end if
    else if (beam_angular_dist .eq. "spherical") then
        if (beam_cutout .eq. "cicle") then
            bx = beam_cutout_params(1) ! axial distance
            by = beam_cutout_params(2) ! full base
            bz = bx/sqrt(bx**2+by**2/4) ! cos(theta)

            aperture = twopi*(1-bz)
        else if (beam_cutout .eq. "rectangle") then
            bx = beam_cutout_params(2)
            by = beam_cutout_params(3)
            bz = beam_cutout_params(1)

            aperture = 4*atan(bx*by/(2*bz*sqrt(4*bz**2+bx**2+by**2)))
        else if (beam_cutout .eq. "none") then ! FINAL UNITS: (s0/beam sr)^(-1) (so, angular density of particles)
            aperture = 1
        end if
    end if

    k0 = beam_axis/norm2(beam_axis)
    bt0 = acos(k0(3))
    bp0 = atan2(k0(2),k0(1))
    bd0 = beam_cutout_params(4)
    beam_coord_sys(:,1) = [cos(bp0)*cos(bt0)*cos(bd0) + sin(bp0)*sin(bd0), &
                           (sin(bp0)*cos(bt0)*cos(bd0) - cos(bp0)*sin(bd0)), &
                           -sin(bt0)*cos(bd0)]
    beam_coord_sys(:,2) = [cos(bp0)*cos(bt0)*sin(bd0) - sin(bp0)*cos(bd0), &
                           (sin(bp0)*cos(bt0)*sin(bd0) + cos(bp0)*cos(bd0)), &
                           -sin(bt0)*sin(bd0)]
    beam_coord_sys(:,3) = k0

    if (command_argument_count() .eq. 1) then
        call get_command_argument(1, output_fname)
    else
        output_fname = "Results/Unnamed run"
    end if
end subroutine prepare_globals

subroutine alpha_vers_warnings
    implicit none
    print *, "WARNING:"

    print *, "You are using the version 0.1 of wiscobolt."
    print *, "Prior to using this for any calculations whatsoever, you should know:"

    print *, "** The 'physics', i.e., the attenuation coefficients and how they have been"
    print *, "implemented, has not yet been verified against test cases. For photons,"
    print *, "the attenuation coefficients themselves can be, and have been, verified in isolation."
    print *, "What then remains is to verify solution methods, of which SN is reliable but PN"
    print *, "is not quite so (see below). Still, it is best to compare with Monte Carlo or known"
    print *, "test cases. For electrons, this is a substantial problem, because of the difficulty"
    print *, "of the problem of electron transport. Electron physics requires a number of"
    print *, "approximations and special treatments so that reasonable convergence can be"
    print *, "achieved. If these treatments are not correctly applied, the solution may be"
    print *, "inaccurate. It is necessary to compare to experimental data or Monte Carlo"
    print *, "simulations, using identical problem geometry and beam quality/modelling,"
    print *, "and this has yet to be done for wiscobolt."

    print *, "** PN angular discretization is not currently available."
    print *, "At the moment, it provides solutions which are not accurate on the volume boundary."
    print *, "The cause of this is likely related to construction of certain matrix elements"
    print *, "involved in the PN system matrix. This is under investigation and should be"
    print *, "resolved shortly."

    print *, "** The program is generally unoptimized, except for SN angular, MGXS energy,"
    print *, "and the SI solution method."

    print *, "** SN FEXS GMRESm is not yet implemented."

    print *, "For these reasons, the developer does not recommend the use of wiscobolt for"
    print *, "research until it advances past version 0.1, which will happen when"
    print *, "all of the physics can be verified."
end subroutine alpha_vers_warnings

subroutine write_input_file
    implicit none
    character(50) :: trashchar
    integer :: io_stat

    open(1, file = "input.in")
    open(2, file = trim(adjustl(output_fname))//"/input.in")

    do
        read(1,'(A)', iostat = io_stat) trashchar
        if (io_stat .ne. 0) exit
        write(2,'(A)') trim(adjustl(trashchar))
    end do

    close(1)
    close(2)
end subroutine write_input_file

real function T_MMS_solution(x,y,z,mu,phi)
    implicit none
    real, intent(in) :: x
    real, intent(in) :: y
    real, intent(in) :: z
    real, intent(in) :: mu
    real, intent(in) :: phi

    real :: A = 100
    real :: alpha = 2*log(2.0+sqrt(3.0))
    real :: angular
    real :: Xt
    real :: Yt
    real :: Zt
    real :: Rad = 1
    real :: falloff = 1.0/1000
    real :: C
    real :: r
    real :: rdotk

    angular = exp(-cos(phi))*exp(mu)
    Xt = exp(-A*(x-0.5)**2)
    Yt = exp(-A*(y-0.5)**2)
    Zt = 2-cosh(alpha*z-alpha/2)
    T_MMS_solution = angular*Xt*Yt*Zt
    MMS_bdy = .false.

    !C = -4*log(falloff)
    !r = sqrt(x**2+y**2+z**2)
    !T_MMS_solution = 1.0/sqrt(fourpi)*&!0.5*sqrt(15/pi)*mu*sqrt(1-mu**2)*sin(phi)*&
    !    exp(-C*r**2)
    !MMS_bdy = .true.

    !C = -log(falloff)/(4*Rad**2)
    !r = sqrt(x**2+y**2+z**2)
    !rdotk = x*cos(phi)*sqrt(1-mu**2) + &
    !        y*sin(phi)*sqrt(1-mu**2) + &
    !        z*mu
    !T_MMS_solution = exp(-C*r**2)*exp(-C*Rad**2)*exp(2*C*Rad*rdotk)
    !MMS_bdy = .true.

    !angular = exp(-cos(phi))*exp(mu)
    !Xt = exp(-A*(x-0.5)**2)
    !Yt = exp(-A*(y-0.5)**2)
    !Zt = 2-cos(alpha*z-alpha/2)
    !T_MMS_solution = angular*Xt*Yt*Zt
    !MMS_bdy = .true.
end function T_MMS_solution

real function T_MMS_source(x,y,z,mu,phi)
    implicit none
    real, intent(in) :: x
    real, intent(in) :: y
    real, intent(in) :: z
    real, intent(in) :: mu
    real, intent(in) :: phi

    real :: psi
    real :: angular
    real :: Xt
    real :: Yt
    real :: Zt1
    real :: Zt2
    real :: A = 100
    real :: alpha = 2*log(2.0+sqrt(3.0))
    real :: Rad = 1
    real :: falloff = 1.0/1000
    real :: C
    real :: rdotk

    psi = T_MMS_solution(x,y,z,mu,phi)

    angular = exp(-cos(phi))*exp(mu)
    Xt = exp(-A*(x-0.5)**2)
    Yt = exp(-A*(y-0.5)**2)
    Zt1 = 2-cosh(alpha*z-alpha/2)
    Zt2 = sinh(alpha*z-alpha/2)
    T_MMS_source = &
        - A*sqrt(1-mu**2)*cos(phi)*(2*x-1)*angular*Xt*Yt*Zt1 &
        - A*sqrt(1-mu**2)*sin(phi)*(2*y-1)*angular*Xt*Yt*Zt1 &
        - mu*alpha*angular*Xt*Yt*Zt2 &
        + MMS_attn*psi

    !C = -4*log(falloff)
    !rdotk = x*cos(phi)*sqrt(1-mu**2) + &
    !        y*sin(phi)*sqrt(1-mu**2) + &
    !        z*mu
    !T_MMS_source = psi*(MMS_attn - 2*C*rdotk)

    !C = -log(falloff)/(4*Rad**2)
    !rdotk = x*cos(phi)*sqrt(1-mu**2) + &
    !        y*sin(phi)*sqrt(1-mu**2) + &
    !        z*mu
    !T_MMS_source = psi*(MMS_attn - 2*C*(rdotk-Rad))

    !angular = exp(-cos(phi))*exp(mu)
    !Xt = exp(-A*(x-0.5)**2)
    !Yt = exp(-A*(y-0.5)**2)
    !Zt1 = 2-cos(alpha*z-alpha/2)
    !Zt2 = sin(alpha*z-alpha/2)
    !T_MMS_source = &
    !    - A*sqrt(1-mu**2)*cos(phi)*(2*x-1)*angular*Xt*Yt*Zt1 &
    !    - A*sqrt(1-mu**2)*sin(phi)*(2*y-1)*angular*Xt*Yt*Zt1 &
    !    + mu*alpha*angular*Xt*Yt*Zt2 &
    !    + MMS_attn*psi
end function T_MMS_source

real function full_MMS_solution(x,y,z,mu,phi)
    implicit none
    real, intent(in) :: x
    real, intent(in) :: y
    real, intent(in) :: z
    real, intent(in) :: mu
    real, intent(in) :: phi

    real :: A = 100
    real :: alpha = 2*log(2.0+sqrt(3.0))
    real :: Rad = 1
    real :: falloff = 1.0/1000
    real :: C
    real :: r
    real :: rdotk

    full_MMS_solution = exp(-cos(phi))*exp(mu)*&
        exp(-A*(x-0.5)**2)*exp(-A*(y-0.5)**2)*(2-cosh(alpha*z-alpha/2))
end function full_MMS_solution

real function full_MMS_non_K_source(x,y,z,mu,phi)
    implicit none
    real, intent(in) :: x
    real, intent(in) :: y
    real, intent(in) :: z
    real, intent(in) :: mu
    real, intent(in) :: phi

    real :: psi
    real :: A = 100
    real :: alpha = 2*log(2.0+sqrt(3.0))
    real :: Rad = 1
    real :: falloff = 1.0/1000
    real :: C
    real :: rdotk

    psi = full_MMS_solution(x,y,z,mu,phi)

    full_MMS_non_K_source = &
        psi*(- A*sqrt(1-mu**2)*cos(phi)*(2*x-1) &
             - A*sqrt(1-mu**2)*sin(phi)*(2*y-1) &
             + MMS_attn) &
             - exp(-cos(phi))*exp(mu)*&
             exp(-A*(x-0.5)**2)*exp(-A*(y-0.5)**2)*mu*alpha*sinh(alpha*z-alpha/2)
end function full_MMS_non_K_source

real function full_MMS_scattering_XS(mu)
    implicit none
    real, intent(in) :: mu

    full_MMS_scattering_XS = exp(mu)
end function full_MMS_scattering_XS

subroutine storage_estimation
    implicit none
    real :: permanent = 0
    real :: temporary = 0
    real :: solution = 0

    !! Aup and Adn, elemental forms only (find way to estimate boundary faces for global forms too? Both will be allocated at same time)
    !permanent = permanent + 2*NE*NQ*NQ*NFe*4
    !
    !! Temporary solution arrays, avg'd over SN and PN
    !temporary = temporary + 0.5*(NE*Nmu*Nphi*NKe+NE*NQ*NKe)*4
    !
    !! psip
    !solution = solution + NK*Nmu*Nphi*Gp*4
    !
    !! psie
    !solution = solution + NK*NQ*Ge*4
    !
    !
    !permanent = permanent*10.0**(-9)
    !temporary = temporary*10.0**(-9)
    !solution = solution*10.0**(-9)

    !print *, "PERMANENT ARRAYS:", permanent, "GB"
    !print *, "TEMPORARY ARRAYS:", temporary, "GB"
    !print *, "TOTAL:", permanent+temporary, "GB"
    !
    !print *, "SOLUTION ONLY:", , "GB"
end subroutine storage_estimation

end module user_input