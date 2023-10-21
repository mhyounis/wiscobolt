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

program solver
    use OMP_LIB
    use math
    use physics
    use user_input
    use mesh
    use shape_functions
    use angular_discretization
    use sweep_order
    use spatial_inner_products
    use energy_discretization
    use MGXS
    use FEXS
    use sources
    use ele_sources
    use iteration
    use t1_ele_iteration
    use t2_ele_iteration
    use post_processing
! gfortran math.f08 -g physics.f08 -g 1_user_input.f08 -g 2_mesh.f08 -g 3_shape_functions.f08 -g 4_angular_discretization.f08 -g 5_sweep_order.f08 -g 6_spatial_inner_products.f08 -g 7_energy_discretization.f08 -g MGXS.f08 -g FEXS.f08 -g 8_sources.f08 -g 8_ele_sources.f08 -g 9_iteration.f08 -g 9_t1_ele_iteration.f08 -g 9_t2_ele_iteration.f08 -g post_processing.f08 -g solver.f08 -g -o wiscobolt -fopenmp -Wl,--stack,100000000 -fdefault-real-8
implicit none

! MAIN
real :: start, finish, firststart, lastfinish
real :: ompstart, ompfinish, firstompstart, lastompfinish, outerompstart, outerompfinish
integer :: dir, f, e, ele, g, gpr, i, j, k, kp, l, lp, m, mat, n, np, particle, q, s

! 2
real, dimension(:,:), allocatable :: rglobal !(dir,k)
integer, dimension(:,:), allocatable :: Cekk
real, dimension(:,:,:), allocatable :: r ! (dir,k,e)
integer, dimension(:,:), allocatable :: Cfkf ! (fg,kf) --> kg
integer, dimension(:,:), allocatable :: Ceff ! (e,f) --> fg
integer, dimension(:,:), allocatable :: eprime ! (f,e)
integer, dimension(:,:), allocatable :: fprime ! (f,e)
integer, dimension(:,:,:), allocatable :: intfc !(fg, 1 for el or 2 for face, 1 for primary or 2 for secondary face), "primary" and "secondary" face is not ordered particularly.
integer, dimension(:,:), allocatable :: bdyfc ! (fbdy, 1 for el or 2 for face or 3 for global face)
integer, dimension(:), allocatable :: bdyel
integer, dimension(:,:,:), allocatable :: kprime ! (k,f,e)
real, dimension(:,:,:), allocatable :: rfcm ! (dir,f,e)
real, dimension(:,:,:), allocatable :: normal ! (dir,f,e)
real, dimension(:,:), allocatable :: globalnormal ! (dir,fg) sign of normal determined by the correct sign for the element with the bigger index.
logical, dimension(:,:), allocatable :: sgn ! (f,e) gives true if (f,e) is "primary" face, false if not.
integer, dimension(:,:), allocatable :: intsgn ! (e,f) gives +1 if (e,f) is "primary" face, -1 if not.
integer, dimension(:,:), allocatable :: matfc ! (fmat,mat)

! 3
real, dimension(:,:), allocatable :: area ! (f,e)
real, dimension(:), allocatable :: vol
real, dimension(:,:), allocatable :: a ! (k,e)
real, dimension(:,:,:), allocatable :: bvector ! (dir,k,e)

! 4
integer :: NQ
real, dimension(:), allocatable :: mu
real, dimension(:), allocatable :: wL
real, dimension(:), allocatable :: phiC
real, dimension(:,:,:), allocatable :: khat ! (dir,i,j)
real, dimension(:,:,:), allocatable :: Pa ! (i,l,m)
real, dimension(:,:,:), allocatable :: cosmat ! (j,jp,m)
real, dimension(:,:), allocatable :: factorialmat ! (l,m)
real, dimension(:,:), allocatable :: qPa ! (i,q)
real, dimension(:), allocatable :: qfactorialmat
!real, dimension(:,:), allocatable :: Awl !
real, dimension(:,:,:), allocatable :: Awl !
!real, dimension(:,:,:,:), allocatable :: Aup ! sign of normal determined by the "primary" face definition in intfc receiving priority. Map between primary and secondary is:
!real, dimension(:,:,:,:), allocatable :: Adn !
real, dimension(:,:,:), allocatable :: Aup
real, dimension(:,:,:), allocatable :: Adn
integer, dimension(:,:), allocatable :: qIa
integer, dimension(:,:), allocatable :: qJa
real, dimension(:,:,:), allocatable :: RSH !
real, dimension(:,:), allocatable :: RSH_unc !
real, dimension(:,:,:), allocatable :: Pmlk ! (l,m,k)
real, dimension(:,:,:), allocatable :: cmjk ! (j,m,k)

! 5
integer, dimension(:,:), allocatable :: sweeplist
integer, dimension(:), allocatable :: sweepbounds
integer :: nosweeps
integer, dimension(:,:,:), allocatable :: esweeplist
integer, dimension(:,:,:), allocatable :: esweepbounds
integer, dimension(:,:), allocatable :: enosweeps
integer, dimension(:,:,:), allocatable :: enew
integer, dimension(:,:,:,:), allocatable :: epprime
integer, dimension(:,:,:,:), allocatable :: fold
integer, dimension(:,:,:), allocatable :: Ia
integer, dimension(:,:,:), allocatable :: Ja
integer :: noJa
integer, dimension(:,:), allocatable :: Cipkk ! (k,ip)
integer, dimension(:,:,:), allocatable :: eijip ! (e,i,j)
integer, dimension(:,:), allocatable :: ipprime ! (f,ip)
real, dimension(:,:,:,:), allocatable :: varsigmaup ! (f,e,i,j)
real, dimension(:,:,:), allocatable :: varsigmadown ! (f,e,i,j)

! 6
real, dimension(:,:,:,:), allocatable :: I1invI2 ! (k,kp,f,e)
real, dimension(:,:,:,:), allocatable :: I1invI2f ! (kp,k,f,e)
real, dimension(:,:,:,:), allocatable :: I1invI3vector ! (k,kp,dir,e)
real, dimension(:,:,:,:), allocatable :: PN_I1invI2 ! (e,kp,f,k)
real, dimension(:,:,:,:), allocatable :: PN_I1invI2f ! (e,kp,f,k)
real, dimension(:,:,:,:), allocatable :: PN_I1invI3vector ! (e,kp,dir,k)

! 7
real, dimension(:), allocatable :: photon_exc
real, dimension(:), allocatable :: electron_exc
real, dimension(:), allocatable :: Ep
real, dimension(:), allocatable :: dEp
real, dimension(:), allocatable :: Epmid
real, dimension(:), allocatable :: Ee
real, dimension(:), allocatable :: dEe
real, dimension(:), allocatable :: Eemid
real, dimension(:), allocatable :: Em
real, dimension(:), allocatable :: dEm
real, dimension(:), allocatable :: Emmid
real, dimension(:), allocatable :: fEspec
real, dimension(:), allocatable :: fE

! MGXS
! Electrons
real, dimension(:,:,:,:), allocatable :: MS_ee ! (gpr,g,l+1,mat)
real, dimension(:,:), allocatable :: MSt_e ! (g,mat)
real, dimension(:,:), allocatable :: MSa_e
real, dimension(:,:), allocatable :: MS2_e
real, dimension(:,:), allocatable :: MEDEP_e
real, dimension(:,:), allocatable :: MCDEP_e
!  REWORK/ RENAME EVERYTHING IN MGXS AND FEXS THAT IS BELOW
! Misc + MMS
real, dimension(:,:,:,:), allocatable :: Sigma
real, dimension(:,:), allocatable :: Sigmat
real, dimension(:), allocatable :: Sigmamt ! (g)
! Photons
real, dimension(:,:,:,:), allocatable :: M_Sigmapp ! (gpr,g,l+1,mat)
real, dimension(:,:), allocatable :: p_Sigmapt ! (g,mat)
real, dimension(:,:,:,:), allocatable :: M_Sigmape ! (gpr,g,l+1,mat)
real, dimension(:,:), allocatable :: M_Sigmapt ! (g,mat)
real, dimension(:,:), allocatable :: MS_pt ! (g,mat)
real, dimension(:,:), allocatable :: MS_pa
real, dimension(:,:), allocatable :: MS_p2
real, dimension(:,:), allocatable :: MEDEP_p
real, dimension(:,:), allocatable :: MCDEP_p

! FEXS
real, dimension(:,:), allocatable :: i_Sigmapt ! (g,mat)
real, dimension(:,:), allocatable :: i_Sigmaet ! (g,mat)
real, dimension(:,:), allocatable :: pi_Sigmaet ! (g,mat)
real, dimension(:,:,:,:), allocatable :: Fp_Sigmapt ! (g,mat)
real, dimension(:,:,:,:), allocatable :: Fp_Sigmaet ! (g,mat)
real, dimension(:,:,:,:,:,:), allocatable :: F_Sigmapp ! (np,n,gpr,g,l+1,mat)
real, dimension(:,:,:,:,:,:), allocatable :: F_Sigmape ! (np,n,gpr,g,l+1,mat)
real, dimension(:,:,:,:,:,:), allocatable :: F_Sigmaee ! (np,n,gpr,g,l+1,mat)
real, dimension(:,:,:,:), allocatable :: F_Sigmapt ! (np,n,g,mat)
real, dimension(:,:,:,:), allocatable :: F_Sigmaet ! (np,n,g,mat)
! FEXS temporary cross sections
real, dimension(:,:,:,:,:), allocatable :: tF_Sigma
real, dimension(:,:,:), allocatable :: tF_Sigmat
real, dimension(:,:,:,:), allocatable :: tF_Sigmagg
! FEXS MMS
real, dimension(:,:,:), allocatable :: F_Sigmamt ! (np,n,g)
! Deposition and other cross sections!
real, dimension(:,:,:), allocatable :: F_Scol
real, dimension(:,:,:), allocatable :: F_muen

! 8
real, dimension(:,:), allocatable :: ells ! (k,mat)
real, dimension(:,:,:), allocatable :: phipu ! (e,k,g)
real, dimension(:,:,:), allocatable :: phieu ! (e,k,g)
real, dimension(:,:,:), allocatable :: noquad_phieu ! (e,k,g). Used to store the phieu prior to correction, this is ACTUALLY phi(r(e,k)), whereas phieu is I1inv*int(u(e,k)*phi(r)) following beam quadrature
real, dimension(:,:,:), allocatable :: onlyquad_phieu ! (e,k,g). Used to store the phieu AFTER correction, but before the ETC correction. Then, take phieu - onlyquad_phieu, to get the portion of the fluence due only to the ETC correction.
real, dimension(:), allocatable :: M_fEg
real, dimension(:,:), allocatable :: F_fEg
real, dimension(:,:,:,:), allocatable :: bdysrc
real, dimension(:,:,:), allocatable :: bdy
! SN MGXS sources
real, dimension(:,:,:,:), allocatable :: SNM_sourcep ! (k,i,j,g)
real, dimension(:,:,:,:), allocatable :: SNM_sourcee ! (k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SNM_sourcep ! (e,k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SNM_sourcee ! (e,k,i,j,g)
! SN FEXS sources
real, dimension(:,:,:,:,:), allocatable :: SNF_sourcep ! (k,i,j,n,g)
real, dimension(:,:,:,:,:), allocatable :: SNF_sourcee ! (k,i,j,n,g)
! PN MGXS temporary sources
real, dimension(:,:,:), allocatable :: tPNM_e_source ! (e,k,q)
real, dimension(:,:), allocatable :: tPNM_source ! (k,q)
! PN FEXS temporary sources
real, dimension(:,:,:,:), allocatable :: tPNF_e_source ! (e,k,q,n)
real, dimension(:,:,:), allocatable :: tPNF_source ! (k,q,n)
! PN MGXS sources
real, dimension(:,:,:), allocatable :: PNM_sourcep ! (k,q,g)
real, dimension(:,:,:), allocatable :: PNM_sourcee ! (k,q,g)
! PN FEXS sources
real, dimension(:,:,:,:), allocatable :: PNF_sourcep ! (k,q,n,g)
real, dimension(:,:,:,:), allocatable :: PNF_sourcee ! (k,q,n,g)
! MMS sources
real, dimension(:,:,:,:), allocatable :: SN_sourcem ! (k,i,j,g)
real, dimension(:,:,:,:), allocatable :: SN_sourcemm ! (k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SN_sourcem ! (e,k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SN_sourcemm ! (e,k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: sweepnormcheck ! (e,k,i,j,g)
real, dimension(:,:,:), allocatable :: PN_sourcem ! (k,q,g)
real, dimension(:,:,:), allocatable :: PN_sourcemm ! (k,q,g)

! 9
integer :: total
real, dimension(:,:,:,:,:), allocatable :: M_Tprimeinv ! (kp,k,e,i,j)
real, dimension(:,:,:,:,:,:,:), allocatable :: F_Tprimeinv ! (kp,k,np,n,e,i,j)
real, dimension(:,:,:), allocatable :: phim ! (e,k,g)
real, dimension(:,:,:), allocatable :: phip ! (e,k,g)
real, dimension(:,:,:), allocatable :: phie ! (e,k,g)
! SN MGXS temporary angular fluences
real, dimension(:,:,:), allocatable :: tSNM_psi ! (k,i,j)
! SN FEXS temporary angular fluences
real, dimension(:,:,:,:), allocatable :: tSNF_psi ! (k,i,j,n)
! SN fluences
real, dimension(:,:,:,:), allocatable :: SN_psip ! (k,i,j,g)
real, dimension(:,:,:,:), allocatable :: SN_psie ! (k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SN_psip ! (e,k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SN_psie ! (e,k,i,j,g)
! PN MGXS temporary fluences
real, dimension(:,:,:), allocatable :: tPNM_e_psi ! (e,k,q)
! PN FEXS temporary fluences
real, dimension(:,:,:,:), allocatable :: tPNF_e_psi ! (e,k,q,n)
! PN fluences
real, dimension(:,:,:), allocatable :: PN_psip ! (k,q,g)
real, dimension(:,:,:), allocatable :: PN_psie ! (k,q,g)
! MMS fluences
real, dimension(:,:,:,:), allocatable :: SN_psim ! (k,i,j,g)
real, dimension(:,:,:,:), allocatable :: SN_psimm ! (k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SN_psim ! (e,k,i,j,g)
real, dimension(:,:,:,:,:), allocatable :: e_SN_psimm ! (e,k,i,j,g)
real, dimension(:,:,:), allocatable :: PN_psim ! (k,q,g)
real, dimension(:,:,:), allocatable :: PN_psimm ! (k,q,g)

! POST-PROCESSING
real, dimension(:,:), allocatable :: g_phipu
real, dimension(:,:), allocatable :: g_phip
real, dimension(:,:), allocatable :: g_phieu
real, dimension(:,:), allocatable :: g_phie
real, dimension(:,:), allocatable :: g_phim
real, dimension(:), allocatable :: DOSE_p
real, dimension(:), allocatable :: ENERGY_p
real, dimension(:), allocatable :: CHARGE_p
real, dimension(:), allocatable :: DOSE_e
real, dimension(:), allocatable :: ENERGY_e
real, dimension(:), allocatable :: CHARGE_e
real, dimension(:), allocatable :: DOSE
real, dimension(:), allocatable :: ENERGY
real, dimension(:), allocatable :: CHARGE
real, dimension(:,:), allocatable :: e_ENERGY_p
real, dimension(:,:), allocatable :: e_DOSE_p
real, dimension(:,:), allocatable :: e_CHARGE_p
real, dimension(:,:), allocatable :: e_ENERGY_e
real, dimension(:,:), allocatable :: e_DOSE_e
real, dimension(:,:), allocatable :: e_CHARGE_e
real, dimension(:,:), allocatable :: e_ENERGY
real, dimension(:,:), allocatable :: e_DOSE
real, dimension(:,:), allocatable :: e_CHARGE
integer, dimension(:), allocatable :: slicek
real, dimension(:,:), allocatable :: slicerglobal
integer, dimension(:,:), allocatable :: sliceCfkf

! TEMP
real, dimension(:), allocatable :: degen
real, dimension(:), allocatable :: W
real, dimension(:), allocatable :: U
real, dimension(:), allocatable :: RCSDA_MMS_f
real, dimension(:), allocatable :: RCSDA_MMS_s
real, dimension(:,:), allocatable :: RCSDA_MMS_solved
real, dimension(:), allocatable :: RCSDA_MMS_solved_1D
real, dimension(:), allocatable :: RCSDA_MMS_s_1D
real, dimension(:,:,:,:,:), allocatable :: B_Sigma
integer :: g1, g2
real :: S1, S2, S3, S4, S5
real, dimension(:), allocatable :: rnd
integer :: Gm
real, dimension(:), allocatable :: xq1, wq1
real, dimension(:), allocatable :: xq2, wq2
real, dimension(:,:), allocatable :: coeff
real, dimension(:), allocatable :: P

! wiscoslab
real, dimension(:,:,:,:), allocatable :: GpFterm
real, dimension(:,:,:,:,:), allocatable :: Fsweep

print *, &
"-------------------------------------------------------------------------------------------------"
print *, &
"### wiscobolt --- Developed by Muhsin Hytham Younis (myounis@wisc.edu). First published 10-21-23."
print *, "VERSION: 0.1 (published 10-21-23)"
print *, "*** Licensed under the GNU General Public License 3.0."
print *, "*** See LICENSE/gpl-3.0.txt or visit https://www.gnu.org/licenses/ for more information."
print *, &
"-------------------------------------------------------------------------------------------------"

! 1
firstompstart = omp_get_wtime()
call alpha_vers_warnings
call read_input_file
call check_for_WIP
call prepare_globals

print *, "INPUT FILE SUCCESSFULLY READ."

if (output_fname .ne. "Results/Unnamed run") call write_input_file
!if (estimate_storage) call storage_estimation

if (verify_beam_spec) then
    Gm = merge(Ge, Gp, transport_mode .eq. "external electron beam")
    Emmin = merge(Eemin, Epmin, transport_mode .eq. "external electron beam")
    Emmax = merge(Eemax, Epmax, transport_mode .eq. "external electron beam")
    electron_energy_structure = merge(electron_energy_structure, &
        photon_energy_structure, &
        transport_mode .eq. "external electron beam")

    if (transport_mode .eq. "external electron beam") then
        call construct_exceptions("electron", electron_exc)
    else
        call construct_exceptions("photon", photon_exc)
        allocate(electron_exc, source = photon_exc)
    end if

    if (electron_energy_structure .eq. "logarithmic") then
        call logarithmic_energy_groups(Gm, Emmin, Emmax, electron_exc, Em, dEm, Emmid)
    else if (electron_energy_structure .eq. "linear") then
        call linear_energy_groups(Gm, Emmin, Emmax, electron_exc, Em, dEm, Emmid)
    else if (electron_energy_structure .eq. "exponential") then
        call exponential_energy_groups(Gm, Emmin, Emmax, electron_exc, Em, dEm, Emmid)
    end if

    if (beam_energy_dist .eq. "polychromatic") then
        call uncollided_beam_energy_distribution(fEspec, fE)
    end if

    call MGXS_discretize_energy_spectrum(Em, fEspec, fE, M_fEg)

    print *, "### Verification of beam energy distribution."
    if (allocated(fE)) then
        print *, "User-defined beam spectrum, energy nodes:"
        do i = 1, size(fE)
            print *, fEspec(i)
        end do
        print *, "User-defined beam spectrum, distribution values:"
        do i = 1, size(fE)
            print *, fE(i)
        end do
    else
        print *, "User has defined non-polychromatic beam."
        print *, "There is no meaningful analog of the un-discretized"
        print *, "beam energy distribution."
    end if
    print *, "MGXS-discretized beam spectrum, energy bin endpoints:"
    do g = 1, Gm
        do n = 1, 2
            print *, Em(g+n-1)
        end do
    end do
    print *, "MGXS-discretized beam spectrum, bin distribution values:"
    do g = 1, Gm
        do n = 1, 2
            print *, M_fEg(g)/dEm(g)
        end do
    end do
    print *, "PROGRAM ENDING."
    stop
end if

if (wiscoslab) go to 1100

print *, "PRE-PROCESSING STARTED."

! 2
ompstart = omp_get_wtime()
if (mesh_reader .eq. "old") then
    NE = 9312
    NK = 2188
    call old_import_nodes(rglobal)
    call old_import_elements(Cekk)
    call old_mesh_analysis(Cekk, Cfkf, Ceff, eprime, fprime, intfc, bdyfc, bdyel)
else if (mesh_reader .eq. "GMSH") then
    call GMSH_mesh_import_and_analysis &
    (rglobal, Cekk, Cfkf, Ceff, eprime, fprime, intfc, bdyfc, bdyel)
end if
call construct_elemental_mesh(rglobal, Cekk, r)
call determine_overlapping_nodes(Cekk)
call find_identical_nodes(Cekk, eprime, kprime)
call scale_and_translate_mesh(rglobal, r)
call construct_normal_vectors(rglobal, r, Cekk, Ceff, eprime, globalnormal, normal, sgn)
call modify_Cfkf(rglobal, r, bdyfc, Cfkf)
call center_of_face(r, rfcm)
if (Nmats .gt. 1) call find_material_interfaces(intfc, bdyfc, matfc)
if (transport_mode .ne. "T MMS" .and. transport_mode .ne. "full MMS") then
    if (wide_open_collimator) then
        if (allocated(nodesinbeam)) deallocate(nodesinbeam)
        allocate(nodesinbeam, source = [(k,k=1,NK)])

        if (allocated(elsinbeam)) deallocate(elsinbeam)
        allocate(elsinbeam, source = [(e,e=1,NE)])
    else
        if (.not. beams_defined_in_mesh) then
            call determine_nodes_in_beam(Cekk, rglobal)
        end if
    end if
end if
ompfinish = omp_get_wtime()
print *, "MODULE 2:", real(ompfinish-ompstart), "seconds"

!call single_material_ray_trace_mesh(rglobal, Cfkf, bdyfc, 1.0E-16, ells)
!call TEMP_MANUAL_ray_tracing(rglobal, coeff)
!
!do k = 1, size(nodesinbeam)
!    if (abs(coeff(k,1) - ells(k,1)) .gt. 1.0E-8) then
!        print *, k, coeff(k,1), ells(k,1)
!    end if
!end do
!
!stop

! 3
ompstart = omp_get_wtime()
if (NKe .eq. 4) then
    call tetrahedral_geometry_and_shape_functions(r, area, vol, a, bvector)
end if
ompfinish = omp_get_wtime()
print *, "MODULE 3:", real(ompfinish-ompstart), "seconds"

! 4
ompstart = omp_get_wtime()
call SN_angular_quadrature_abscissae_and_weights(Nmu, Nphi, mu, wL, phiC)
call construct_discrete_ordinates(mu, phiC, khat)
call populate_real_spherical_harmonics(NL, mu, phiC, RSH)
call uncollided_RSH_terms(NL, rglobal, RSH_unc)
if ((photon_angular .eq. "PN" .and. transport_mode .eq. "external photon beam") .or. &
    (electron_angular .eq. "PN" .and. transport_mode .eq. "external photon beam coupled") .or. &
    (electron_angular .eq. "PN" .and. transport_mode .eq. "external electron beam") .or. &
    (MMS_angular .eq. "PN" .and. transport_mode .eq. "T MMS") .or. &
    (MMS_angular .eq. "PN" .and. transport_mode .eq. "full MMS")) then
    call PN_angular_integrals_via_quadrature &
        (NL, Ceff, globalnormal, sgn, 4*Nmu, Awl, qIa, qJa, Aup, Adn)
end if
if (photon_angular .eq. "SN" .or. "electron_angular" .eq. "SN" .or. MMS_angular .eq. "SN") then
    call SN_scattering_angular_quadrature_terms &
        (mu, phiC, cosmat, Pa, factorialmat, qPa, qfactorialmat)
    call SN_uncollided_angular_quadrature_terms(rglobal, phiC, Pmlk, cmjk)
end if
NQ = (NL+1)**2
ompfinish = omp_get_wtime()
print *, "MODULE 4:", real(ompfinish-ompstart), "seconds"

! 5
if (photon_angular .eq. "SN" .or. "electron_angular" .eq. "SN" .or. MMS_angular .eq. "SN") then
    ompstart = omp_get_wtime()
    if (fresh_sweep_run) then
        call construct_esweeplist(eprime, fprime, normal, khat, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa)
    else
        call read_esweeplist(esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa)
    end if
    call construct_varsigmas(normal, khat, esweeplist, fold, Ia, noJa, varsigmaup, varsigmadown)
    ompfinish = omp_get_wtime()
    print *, "MODULE 5:", real(ompfinish-ompstart), "seconds"
else
    print *, "NOTHING TO DO WITH MODULE 5"
end if

if (verify_sweep) then
    call esweeplist_verification(r, eprime, normal, khat, esweeplist)
    stop
end if

! 6
ompstart = omp_get_wtime()
call inner_products_with_I1inv &
(eprime, fprime, kprime, area, vol, bvector, I1invI2, I1invI2f, I1invI3vector)
ompfinish = omp_get_wtime()
print *, "MODULE 6:", real(ompfinish-ompstart), "seconds"

print *, "PRE-PROCESSING COMPLETE."

! DEALLOCATIONS
deallocate(fprime)
deallocate(kprime)

if (transport_mode .eq. "T MMS" .or. transport_mode .eq. "full MMS") then
    if (transport_mode .eq. "T MMS") then
        print *, "MANUFACTURED SOLUTION STARTED - TRANSPORT OPERATOR."
    else
        print *, "MANUFACTURED SOLUTION STARTED - BOLTZMANN OPERATOR."
    end if

    NQ = (NL+1)**2
    allocate(Sigmat(1,1))

    if (transport_mode .eq. "T MMS") then
        allocate(Sigma(1,1,NL+1,1))
        Sigma = 0.0
    end if
    Sigmat = MMS_attn
    if (.not. allocated(eltomat)) allocate(eltomat(NE))
    eltomat = 1

    ompstart = omp_get_wtime()
    if (transport_mode .eq. "T MMS") then
        if (MMS_angular .eq. "SN") then
            call SN_construct_T_MMS_vector(T_MMS_source, rglobal, mu, phiC, SN_sourcem)
            call SN_construct_T_MMS_vector(T_MMS_solution, rglobal, mu, phiC, SN_psim)
            if (FEM_type .eq. "T2") then
                allocate(e_SN_sourcem(NE,NKe,Nmu,Nphi,1))
                allocate(e_SN_psim(NE,NKe,Nmu,Nphi,1))
                do k = 1, NKe
                    do e = 1, NE
                        e_SN_sourcem(e,k,:,:,1) = SN_sourcem(Cekk(e,k),:,:,1)
                        e_SN_psim(e,k,:,:,1) = SN_psim(Cekk(e,k),:,:,1)
                    end do
                end do
                deallocate(SN_sourcem)
                deallocate(SN_psim)
            end if
            if (MMS_bdy) then
                call SN_boundary_source &
                    (Cekk, r, eprime, bdyel, normal, mu, phiC, khat, I1invI2, bdysrc)
            end if
        else if (MMS_angular .eq. "PN") then
            !call PN_construct_T_MMS_vector(T_MMS_source, rglobal, Cekk, mu, wL, phiC, RSH, PN_e_sourcem)
            !call PN_construct_T_MMS_vector(T_MMS_solution, rglobal, Cekk, mu, wL, phiC, RSH, PN_e_psim)
        end if
    else
        if (MMS_angular .eq. "SN") then
            call SN_construct_full_MMS_vectors(rglobal, mu, wL, phiC, RSH, Sigma, SN_psim, SN_sourcem)
            if (FEM_type .eq. "T2") then
                allocate(e_SN_sourcem(NE,NKe,Nmu,Nphi,1))
                allocate(e_SN_psim(NE,NKe,Nmu,Nphi,1))
                do k = 1, NKe
                    do e = 1, NE
                        e_SN_sourcem(e,k,:,:,1) = SN_sourcem(Cekk(e,k),:,:,1)
                        e_SN_psim(e,k,:,:,1) = SN_psim(Cekk(e,k),:,:,1)
                    end do
                end do
                deallocate(SN_sourcem)
                deallocate(SN_psim)
            end if
        else if (MMS_angular .eq. "PN") then
            print *, "WIP: code 1"
            stop
        end if
    end if
    ompfinish = omp_get_wtime()
    print *, "MODULE 8:", real(ompfinish - ompstart), "seconds"

    ompstart = omp_get_wtime()
    if (transport_mode .eq. "T MMS") then
        if (MMS_angular .eq. "SN") then
            allocate(M_Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

            call SN_MGXS_construct_inverted_transport_matrix &
                (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat(1,:), M_Tprimeinv)

            if (FEM_type .eq. "default") then
                allocate(SN_psimm, source = SN_sourcem)
                allocate(SN_sourcemm, source = SN_psim)
                call SN_MGXS_sweep &
                    (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
                    I1invI2f, bdysrc, M_Tprimeinv, SN_psimm(:,:,:,1))
                deallocate(M_Tprimeinv)
                call SN_MGXS_T &
                    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, &
                    Sigmat(1,:), bdysrc, SN_sourcemm(:,:,:,1))
            else if (FEM_type .eq. "T2") then
                allocate(e_SN_psimm, source = e_SN_sourcem)
                !allocate(e_SN_sourcemm, source = e_SN_psim)
                call t2_ele_SN_MGXS_sweep &
                    (Cekk, bdyel, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, varsigmadown, &
                    I1invI2f, bdysrc, M_Tprimeinv, e_SN_psimm(:,:,:,:,1))
                deallocate(M_Tprimeinv)
                !call ele_SN_MGXS_T &
                !    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, &
                !    Sigmat(1,:), bdysrc, e_SN_sourcemm(:,:,:,:,1))

                allocate(sweepnormcheck, source = e_SN_psimm)
                call ele_SN_MGXS_T &
                    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, &
                    Sigmat(1,:), bdysrc, sweepnormcheck(:,:,:,:,1))

                print *, "SWEEP NORM CHECK: |s - T*T^(-1)*s|_2"
                print *, norm2(e_SN_sourcem - sweepnormcheck)

                stop

                allocate(SN_sourcem(NK,Nmu,Nphi,1))
                allocate(SN_sourcemm(NK,Nmu,Nphi,1))
                allocate(SN_psim(NK,Nmu,Nphi,1))
                allocate(SN_psimm(NK,Nmu,Nphi,1))
                SN_sourcem = 0
                SN_sourcemm = 0
                SN_psim = 0
                SN_psimm = 0
                do k = 1, NKe
                    do e = 1, NE
                        SN_sourcem(Cekk(e,k),:,:,1) = SN_sourcem(Cekk(e,k),:,:,1) + &
                            e_SN_sourcem(e,k,:,:,1)/Nglobal(Cekk(e,k))
                        SN_sourcemm(Cekk(e,k),:,:,1) = SN_sourcemm(Cekk(e,k),:,:,1) + &
                            e_SN_sourcemm(e,k,:,:,1)/Nglobal(Cekk(e,k))
                        SN_psim(Cekk(e,k),:,:,1) = SN_psim(Cekk(e,k),:,:,1) + &
                            e_SN_psim(e,k,:,:,1)/Nglobal(Cekk(e,k))
                        SN_psimm(Cekk(e,k),:,:,1) = SN_psimm(Cekk(e,k),:,:,1) + &
                            e_SN_psimm(e,k,:,:,1)/Nglobal(Cekk(e,k))
                    end do
                end do
                deallocate(e_SN_sourcem)
                deallocate(e_SN_sourcemm)
                deallocate(e_SN_psim)
                deallocate(e_SN_psimm)
            end if

            !print *, "SOLUTION:"
            !do j = 1, Nphi
            !    do i = 1, Nmu
            !        print *, [i,j], norm2(SN_psim(:,i,j,1) - SN_psimm(:,i,j,1))
            !    end do
            !end do
            !print *, "SOURCE:"
            !do j = 1, Nphi
            !    do i = 1, Nmu
            !        print *, [i,j], norm2(SN_sourcem(:,i,j,1) - SN_sourcemm(:,i,j,1))
            !    end do
            !end do

        else if (MMS_angular .eq. "PN") then
            ! REMEMBER THESE!!
            allocate(intsgn, source = merge(1,-1,transpose(sgn)))
            I1invI2 = reshape(I1invI2, [NE,NKe,NFe,NKe], order = [4,2,3,1]) ! (e,kp,f,k)
            I1invI2f = reshape(I1invI2f, [NE,NKe,NFe,NKe], order = [2,4,3,1]) ! (e,kp,f,k)
            I1invI3vector = reshape(I1invI3vector, [NE,NKe,3,NKe], order = [4,2,3,1]) ! (e,kp,dir,k)
            stop
        end if
    else
        if (MMS_angular .eq. "SN") then
            allocate(M_Tprimeinv(NKe,NKe,NE,Nmu,Nphi))

            call SN_MGXS_construct_inverted_transport_matrix &
                (khat, I1invI2, I1invI3vector, enew, varsigmaup, Sigmat(1,:), M_Tprimeinv)

            if (FEM_type .eq. "default") then
                call SN_MGXS_EI_GMRESm &
                    (3, m_pmax, m_mmax, 1, Cekk, bdyel, wL, khat, cosmat, qPa, &
                    qfactorialmat, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
                    varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, Sigma, &
                    Sigmat, SN_sourcem, M_fEg, bdysrc, SN_psimm, phim)

                allocate(SN_sourcemm, source = SN_psim)
                call SN_MGXS_T &
                    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, &
                    Sigmat(1,:), bdysrc, SN_sourcemm(:,:,:,1))
                allocate(PN_psie(NK,Nmu,Nphi)) ! A temporary array
                PN_psie = SN_psim(:,:,:,1)
                call SN_MGXS_Kgg &
                    (wL, cosmat, qPa, qfactorialmat, Sigma(1,1,:,:), PN_psie)
                SN_sourcemm(:,:,:,1) = SN_sourcemm(:,:,:,1) - PN_psie
            else if (FEM_type .eq. "T2") then
                call t2_ele_SN_MGXS_EI_GMRESm &
                    (3, m_pmax, m_mmax, 1, Cekk, bdyel, wL, khat, cosmat, qPa, &
                    qfactorialmat, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
                    varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, Sigma, &
                    Sigmat, e_SN_sourcem, M_fEg, bdysrc, e_SN_psimm, phim)

                !allocate(e_SN_sourcemm, source = e_SN_psim)
                !call ele_SN_MGXS_T &
                !    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, &
                !    Sigmat(1,:), bdysrc, e_SN_sourcemm(:,:,:,:,1))
                !allocate(SN_psie(NE,NKe,Nmu,Nphi)) ! A temporary array
                !SN_psie = e_SN_psim(:,:,:,:,1)
                !call SN_MGXS_Kgg &
                !    (wL, cosmat, qPa, qfactorialmat, Sigma(1,1,:,:), SN_psie)
                !e_SN_sourcemm(:,:,:,1) = e_SN_sourcemm(:,:,:,1) - SN_psie

                allocate(sweepnormcheck, source = e_SN_psimm)
                call ele_SN_MGXS_T &
                    (Cekk, eprime, bdyel, normal, khat, varsigmaup, I1invI2, I1invI2f, I1invI3vector, &
                    Sigmat(1,:), bdysrc, sweepnormcheck(:,:,:,:,1))
                allocate(SN_psie(NE,NKe,Nmu,Nphi)) ! A temporary array
                SN_psie = e_SN_psimm(:,:,:,:,1)
                call t2_ele_SN_MGXS_Kgg &
                    (wL, cosmat, qPa, qfactorialmat, Sigma(1,1,:,:), SN_psie)
                sweepnormcheck(:,:,:,:,1) = sweepnormcheck(:,:,:,:,1) - SN_psie

                print *, "SWEEP NORM CHECK: |s - T*T^(-1)*s|_2"
                print *, norm2(e_SN_sourcem - sweepnormcheck)

                stop

            end if

            !print *, "SOLUTION:"
            !do j = 1, Nphi
            !    do i = 1, Nmu
            !        print *, [i,j], norm2(SN_psim(:,i,j,1) - SN_psimm(:,i,j,1))
            !    end do
            !end do
            !stop
            !print *, "SOLUTION:"
            !do j = 1, Nphi
            !    do i = 1, Nmu
            !        print *, [i,j], norm2(SN_psim(:,i,j,1))/norm2(SN_psimm(:,i,j,1))
            !    end do
            !end do
            !stop
        else if (MMS_angular .eq. "PN") then
            print *, "WIP"
            stop
        end if
    end if
    ompfinish = omp_get_wtime()
    print *, "MODULE 9:", real(ompfinish - ompstart), "seconds"
else
    if (energy_discretization_method .eq. "MGXS") then
        XSlogic = .false.
    else if (energy_discretization_method .eq. "FEXS") then
        XSlogic = .true.
    end if

    if (.not. solve_photons) go to 2
    print *, "PHOTON TRANSPORT STARTED."

    ! 7 - photons
    ompstart = omp_get_wtime()
    call construct_exceptions("photon", photon_exc)
    if (photon_energy_structure .eq. "logarithmic") then
        call logarithmic_energy_groups(Gp, Epmin, Epmax, photon_exc, Ep, dEp, Epmid)
    else if (photon_energy_structure .eq. "linear") then
        call linear_energy_groups(Gp, Epmin, Epmax, photon_exc, Ep, dEp, Epmid)
    else if (photon_energy_structure .eq. "exponential") then
        call exponential_energy_groups(Gp, Epmin, Epmax, photon_exc, Ep, dEp, Epmid)
    end if
    if (beam_energy_dist .eq. "polychromatic") then
        call uncollided_beam_energy_distribution(fEspec, fE)
    end if
    ompfinish = omp_get_wtime()
    print *, "MODULE 7 - PHOTONS:", real(ompfinish-ompstart), "seconds"

    if (customXS) then
        ! MGXS - photons
        ompstart = omp_get_wtime()
        call MGXS_read_custom_XS(1, M_Sigmapp, M_Sigmapt, MEDEP_p, MCDEP_p)
        print *, "MODULE MGXS - PHOTONS:", real(ompfinish-ompstart), "seconds"
        ompfinish = omp_get_wtime()
    else if (isMGXS) then
        ! MGXS - photons
        ompstart = omp_get_wtime()
        call MGXS_photon_photon(Ep, dEp, .true., 1.0E-6, 64, M_Sigmapp, MEDEP_p)
        call MGXS_photon_attn(Ep, dEp, .true., 1.0E-6, 64, M_Sigmapt, MEDEP_p)
        ompfinish = omp_get_wtime()
        print *, "MODULE MGXS - PHOTONS:", real(ompfinish-ompstart), "seconds"
    else if (isFEXS) then
        ! FEXS - photons
        ompstart = omp_get_wtime()
        call FEXS_photon_photon(Ep, dEp, .true., 1.0E-6, 64, F_Sigmapp)
        call FEXS_photon_attn(Ep, dEp, .true., 1.0E-6, 64, F_Sigmapt, i_Sigmapt)
        ompfinish = omp_get_wtime()
        print *, "MODULE FEXS - PHOTONS:", real(ompfinish-ompstart), "seconds"
    end if

    ! 8 - photons
    ompstart = omp_get_wtime()
    isSN = photon_angular .eq. "SN"
    isPN = .not. isSN
    if (.not. fresh_photon_run) go to 7 ! Having this come AFTER setting of isSN and isPN allows to use stored photons and fresh electrons, because fresh electron src. relies on these. Right???
    ! ### SN MGXS
    if      (isMGXS .and. isSN .and. (.not. FCS)) then
        call MGXS_discretize_energy_spectrum(Ep, fEspec, fE, M_fEg)
        call SN_external_beam_boundary_source &
            (Cekk, r, eprime, bdyel, normal, khat, RSH, RSH_unc, I1invI2, bdysrc)
    ! ### PN MGXS
    else if (isMGXS .and. isPN .and. (.not. FCS)) then
        print *, "WIP: PN external beam boundary source"
        stop
    ! ### SN MGXS FCS
    else if (isMGXS .and. isSN .and. FCS) then
        if (Nmats .eq. 1) then
            call single_material_ray_trace_mesh(rglobal, Cfkf, bdyfc, 1.0E-8, ells)
        else
            !call heterogeneous_medium_ray_trace_mesh(rglobal, Cfkf, matfc, 1.0E-8, ells)
            call TEMP_MANUAL_ray_tracing(rglobal, ells)
        end if
        call MGXS_discretize_energy_spectrum(Ep, fEspec, fE, M_fEg)
        call MGXS_construct_uncollided_fluence &
            ("photon", Cekk, r, ells, M_Sigmapt, M_fEg, phipu)
        if (FEM_type .eq. "default") then
            call SN_MGXS_construct_uncollided_source &
                (Cekk, Pa, factorialmat, M_Sigmapp, Pmlk, cmjk, phipu, SNM_sourcep)
        else if (FEM_type .eq. "T2") then
            call t2_ele_SN_MGXS_construct_uncollided_source &
                (Cekk, Pa, factorialmat, M_Sigmapp, Pmlk, cmjk, phipu, e_SNM_sourcep)
        end if
    ! ### PN MGXS FCS ! HAS WIP (T2)
    else if (isMGXS .and. isPN .and. FCS) then
        if (Nmats .eq. 1) then
            call single_material_ray_trace_mesh(rglobal, Cfkf, bdyfc, 1.0E-8, ells)
        else
            !call heterogeneous_medium_ray_trace_mesh(rglobal, Cfkf, matfc, 1.0E-8, ells)
            call TEMP_MANUAL_ray_tracing(rglobal, ells)
        end if
        call MGXS_discretize_energy_spectrum(Ep, fEspec, fE, M_fEg)
        call MGXS_construct_uncollided_fluence &
            ("photon", Cekk, r, ells, M_Sigmapt, M_fEg, phipu)
        if (FEM_type .eq. "default") then
            call PN_MGXS_construct_uncollided_source &
                (NL, Cekk, wL, RSH_unc, RSH, M_Sigmapp, phipu, PNM_sourcep)
        else if (FEM_type .eq. "T2") then
            print *, "WIP: PN T2 ele"
        end if
        call MGXS_discretize_energy_spectrum(Ep, fEspec, fE, M_fEg)
    ! ### SN FEXS ! HAS WIP (T1 and T2)
    else if (isFEXS .and. isSN) then
        call FEXS_discretize_energy_spectrum(Ep, fEspec, fE, F_fEg)
        call SN_external_beam_boundary_source &
            (Cekk, r, eprime, bdyel, normal, khat, RSH, RSH_unc, I1invI2, bdysrc)
    ! ### PN FEXS ! HAS WIP
    else if (isFEXS .and. isPN) then
        print *, "WIP: PN external beam boundary source"
        stop
    end if
    ompfinish = omp_get_wtime()
    print *, "MODULE 8 - PHOTONS:", real(ompfinish-ompstart), "seconds"

    ! 9 - photons
    outerompstart = omp_get_wtime()
    isSI = photon_solution_method .eq. "SI"
    isGMRESm = .not. isSI
    ! ### SN MGXS SI
    if      (isMGXS .and. isSN .and. isSI) then
        if (FEM_type .eq. "default") then
            call SN_MGXS_EI_SI &
                (1, p_pmax, Gp, Cekk, bdyel, eprime, normal, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, M_Sigmapp, M_Sigmapt, phipu, &
                SNM_sourcep, M_fEg, bdysrc, SN_psip, phip)
        else if (FEM_type .eq. "T1") then
            call t1_ele_SN_MGXS_EI_SI &
                (1, p_pmax, Gp, Cekk, bdyel, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
                varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, &
                M_Sigmapp, M_Sigmapt, phipu, SNM_sourcep, M_fEg, bdysrc, SN_psip, phip)
        else if (FEM_type .eq. "T2") then
            call t2_ele_SN_MGXS_EI_SI &
                (1, p_pmax, Gp, Cekk, bdyel, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, M_Sigmapp, M_Sigmapt, phipu, &
                e_SNM_sourcep, M_fEg, bdysrc, e_SN_psip, phip)
        end if
    ! ### SN MGXS GMRESm
    else if (isMGXS .and. isSN .and. isGMRESm) then
        if (FEM_type .eq. "default") then
            call SN_MGXS_EI_GMRESm &
                (1, p_pmax, p_mmax, Gp, Cekk, bdyel, wL, khat, cosmat, &
                qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, M_Sigmapp, M_Sigmapt, &
                SNM_sourcep, M_fEg, bdysrc, SN_psip, phip)
        else if (FEM_type .eq. "T1") then
            call t1_ele_SN_MGXS_EI_GMRESm &
                (1, p_pmax, p_mmax, Gp, Cekk, bdyel, wL, khat, &
                cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, I1invI2f, &
                I1invI3vector, M_Sigmapp, M_Sigmapt, SNM_sourcep, M_fEg, bdysrc, &
                SN_psip, phip)
        else if (FEM_type .eq. "T2") then
            call t2_ele_SN_MGXS_EI_GMRESm &
                (1, p_pmax, p_mmax, Gp, Cekk, bdyel, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, M_Sigmapp, M_Sigmapt, &
                e_SNM_sourcep, M_fEg, bdysrc, e_SN_psip, phip)
        end if
    ! ### PN MGXS GMRESm ! HAS WIP
    else if (isMGXS .and. isPN .and. isGMRESm) then
        print *, "WIP: PHOTON PN MGXS GMRESm"
    ! ### SN FEXS SI ! HAS WIP (T1 and T2)
    else if (isFEXS .and. isSN .and. isSI) then
        if (FEM_type .eq. "default") then
            call SN_FEXS_EI_SI &
                (1, p_pmax, Gp, Cekk, bdyel, wL, khat, cosmat, qPa, qfactorialmat, &
                esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
                varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, &
                F_Sigmapp, F_Sigmapt, SNF_sourcep, F_fEg, bdysrc, SN_psip, phip)
        else
            print *, "WIP: SN FEXS SI T1 and T2"
        end if
    ! ### SN FEXS GMRESm ! HAS WIP
    else if (isFEXS .and. isSN .and. isGMRESm) then
        print *, "WIP: PHOTON SN FEXS GMRESm"
    ! ### PN FEXS GMRESm ! HAS WIP
    else if (isFEXS .and. isPN .and. isGMRESm) then
        print *, "WIP: PHOTON PN FEXS GMRESm"
    end if
    go to 8

    7 continue
    print *, "NOTHING TO DO WITH MODULE 8 - PHOTONS"
    outerompstart = omp_get_wtime()
    ! ### SN
    if (isSN) then
        if (FEM_type .eq. "default") then
            call SN_read_field("photon", merge(Gp, Gp+1, isMGXS), phipu, phip, SN_psip)
        else if (FEM_type .eq. "T2") then
            call t2_ele_SN_read_field("photon", merge(Gp, Gp+1, isMGXS), phipu, phip, e_SN_psip)
        end if
    ! ### PN
    else if (isPN) then
        if (FEM_type .eq. "default") then
            call PN_read_field("photon", merge(Gp, Gp+1, isMGXS), NL, phipu, phip, PN_psip)
        else if (FEM_type .eq. "T2") then
            print *, "WIP: PN T2"
            stop
        end if
    end if
    print *, "PHOTON TRANSPORT READ."

    8 continue
    outerompfinish = omp_get_wtime()
    print *, "MODULE 9 - PHOTONS:", real(outerompfinish-outerompstart), "seconds"
    print *, "PHOTON TRANSPORT COMPLETE."

    if (allocated(SNM_sourcep)) deallocate(SNM_sourcep)
    if (allocated(SNF_sourcep)) deallocate(SNF_sourcep)
    if (allocated(PNM_sourcep)) deallocate(PNM_sourcep)
    if (allocated(PNF_sourcep)) deallocate(PNF_sourcep)
    if (allocated(bdysrc)) deallocate(bdysrc)
    if (allocated(M_fEg)) deallocate(M_fEg)
    if (allocated(F_fEg)) deallocate(F_fEg)

    if (.not. solve_electrons) go to 1

    2 continue
    print *, "ELECTRON TRANSPORT STARTED."

    ! 7 - electrons
    ompstart = omp_get_wtime()
    call construct_exceptions("electron", electron_exc)
    if (electron_energy_structure .eq. "logarithmic") then
        call logarithmic_energy_groups(Ge, Eemin, Eemax, electron_exc, Ee, dEe, Eemid)
    else if (electron_energy_structure .eq. "linear") then
        call linear_energy_groups(Ge, Eemin, Eemax, electron_exc, Ee, dEe, Eemid)
    else if (electron_energy_structure .eq. "exponential") then
        call exponential_energy_groups(Ge, Eemin, Eemax, electron_exc, Ee, dEe, Eemid)
    end if
    if (.not. solve_photons) then
        if (beam_energy_dist .eq. "polychromatic") then
            call uncollided_beam_energy_distribution(fEspec, fE)
        end if
    end if
    ompfinish = omp_get_wtime()
    print *, "MODULE 7 - ELECTRONS:", real(ompfinish-ompstart), "seconds"

    if (customXS) then
        ! MGXS - photons
        ompstart = omp_get_wtime()
        if (solve_photons .and. solve_electrons) then
            call MGXS_read_custom_XS(12, M_Sigmape, M_Sigmapt, MEDEP_p, MCDEP_p)
        end if
        call MGXS_read_custom_XS(2, MS_ee, MSt_e, MEDEP_e, MCDEP_e)
        print *, "MODULE MGXS - PHOTONS:", real(ompfinish-ompstart), "seconds"
        ompfinish = omp_get_wtime()
    else if (isMGXS) then
        ! MGXS - electrons
        ompstart = omp_get_wtime()
        if (solve_photons .and. solve_electrons) then
            call MGXS_photon_electron &
                (Ep, dEp, Ee, dEe, Eemid, .true., 1.0E-8, 64, M_Sigmape, MEDEP_p)
        end if
        call MGXS_electron_electron &
            (Ee, dEe, Eemid, .true., 1.0E-8, 64, MS_ee, MSt_e, MSa_e, MS2_e, MEDEP_e, MCDEP_e)
        ompfinish = omp_get_wtime()
        print *, "MODULE MGXS - ELECTRONS:", real(ompfinish-ompstart), "seconds"
    else if (isFEXS) then
        ! FEXS - electrons
        ompstart = omp_get_wtime()
        if (solve_photons .and. solve_electrons) then
            call FEXS_photon_electron(Ep, dEp, Ee, dEe, Eemid, .true., 1.0E-8, 64, F_Sigmape)
        end if
        call FEXS_electron_electron(Ee, dEe, Eemid, .false., 1.0E-2, 16, F_Sigmaee, Fp_Sigmaet, pi_Sigmaet)
        call FEXS_electron_attn(Ee, dEe, Fp_Sigmaet, pi_Sigmaet, .false., 1.0E-4, 16, F_Sigmaet, i_Sigmaet)
        ompfinish = omp_get_wtime()
        print *, "MODULE FEXS - ELECTRONS:", real(ompfinish-ompstart), "seconds"
    end if

    ! 8 - electrons
    ompstart = omp_get_wtime()
    t_isSN = isSN
    t_isPN = isPN
    isSN = electron_angular .eq. "SN"
    isPN = .not. isSN
    if (.not. fresh_electron_run) go to 9
    if (.not. solve_photons) then
        ! ### SN MGXS
        if      (isMGXS .and. isSN .and. (.not. FCS)) then
            call MGXS_discretize_energy_spectrum(Ee, fEspec, fE, M_fEg)
            call SN_external_beam_boundary_source &
                (Cekk, r, eprime, bdyel, normal, khat, RSH, RSH_unc, I1invI2, bdysrc)
        ! ### PN MGXS
        else if (isMGXS .and. isPN .and. (.not. FCS)) then
            print *, "WIP: PN external beam boundary source"
            stop
        ! ### SN MGXS FCS
        else if (isMGXS .and. isSN .and. FCS) then
            if (Nmats .eq. 1) then
                call single_material_ray_trace_mesh(rglobal, Cfkf, bdyfc, 1.0E-8, ells)
            else
                !call heterogeneous_medium_ray_trace_mesh(rglobal, Cfkf, matfc, 1.0E-8, ells)
                call TEMP_MANUAL_ray_tracing(rglobal, ells)
            end if
            call MGXS_discretize_energy_spectrum(Ee, fEspec, fE, M_fEg)
            call MGXS_construct_uncollided_fluence &
                ("electron", Cekk, r, ells, MSt_e, M_fEg, phieu)
            allocate(noquad_phieu, source = phieu)
            call MGXS_beam_quadrature &
                (Cekk, rglobal, Cfkf, bdyfc, matfc, MSt_e, ells, M_fEg, 4, phieu)
            if (inelastic_ETC .or. exact_RCSDA_angular) then
                allocate(onlyquad_phieu, source = phieu)
                call ele_SN_MGXS_ETC_uncollided_fluence &
                    (rglobal, Cekk, bdyel, wL, khat, Pa, factorialmat, Pmlk, cmjk, esweeplist, &
                    esweepbounds, enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                    I1invI2f, I1invI3vector, MS_ee(:,:,NL+2,:), MSt_e, phieu)
                noquad_phieu = noquad_phieu + phieu - onlyquad_phieu
                deallocate(onlyquad_phieu)
            end if

            ! If I unstrike these, keep MGXS_discretize_energy_spectrum unstruck
            !call SN_external_beam_boundary_source &
            !    (Cekk, r, eprime, bdyel, normal, khat, RSH, RSH_unc, I1invI2, bdysrc)
            !call temp_TESTING &
            !    (rglobal, Cekk, bdyel, wL, khat, Pa, factorialmat, Pmlk, cmjk, esweeplist, &
            !    esweepbounds, enew, fold, Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
            !    I1invI2f, I1invI3vector, MS_ee(:,:,NL+2,:), bdysrc, M_fEg, MSt_e, phieu, noquad_phieu)
            !
            !go to 123

            if (FEM_type .eq. "default") then
                call SN_MGXS_construct_uncollided_source &
                    (Cekk, Pa, factorialmat, MS_ee, Pmlk, cmjk, phieu, SNM_sourcee)
            else if (FEM_type .eq. "T2") then
                call t2_ele_SN_MGXS_construct_uncollided_source &
                    (Cekk, Pa, factorialmat, MS_ee, Pmlk, cmjk, phieu, e_SNM_sourcee)
            end if
        ! ### PN MGXS FCS ! HAS WIP
        else if (isMGXS .and. isPN .and. FCS) then
            print *, "WIP: PN MGXS ETC unc. fluence correction!! Actually, am i supposed to do this for PN?"
            print *, "may want to just make PN MGXS construct unc source do the delta down terms"
            stop
            if (Nmats .eq. 1) then
                call single_material_ray_trace_mesh(rglobal, Cfkf, bdyfc, 1.0E-8, ells)
            else
                !call heterogeneous_medium_ray_trace_mesh(rglobal, Cfkf, matfc, 1.0E-8, ells)
                call TEMP_MANUAL_ray_tracing(rglobal, ells)
            end if
            call MGXS_discretize_energy_spectrum(Ee, fEspec, fE, M_fEg)
            call MGXS_construct_uncollided_fluence &
                ("electron", Cekk, r, ells, MSt_e, M_fEg, phieu)
            allocate(noquad_phieu, source = phieu)
            call MGXS_beam_quadrature &
                (Cekk, rglobal, Cfkf, bdyfc, matfc, MSt_e, ells, M_fEg, 4, phieu)
            allocate(onlyquad_phieu, source = phieu)

            if (FEM_type .eq. "default") then
                call PN_MGXS_construct_uncollided_source &
                    (NL, Cekk, wL, RSH_unc, RSH, MS_ee, phieu, PNM_sourcee)
            else if (FEM_type .eq. "T2") then
                print *, "WIP: PN T2 ele"
            end if
            call MGXS_discretize_energy_spectrum(Ep, fEspec, fE, M_fEg)
        ! ### SN FEXS ! HAS WIP (T1 and T2)
        else if (isFEXS .and. isSN) then
            call FEXS_discretize_energy_spectrum(Ee, fEspec, fE, F_fEg)
            call SN_external_beam_boundary_source &
                (Cekk, r, eprime, bdyel, normal, khat, RSH, RSH_unc, I1invI2, bdysrc)
        ! ### PN FEXS ! HAS WIP
        else if (isFEXS .and. isPN) then
            print *, "WIP: PN external beam boundary source"
            stop
        end if
    else if (solve_photons .and. solve_electrons) then
        ! ### SN to SN MGXS ! HAS WIP (T1)
        if      (isMGXS .and. t_isSN .and. isSN) then
            if (FEM_type .eq. "default") then
                call SN_to_SN_MGXS_particle_particle_scattering_source &
                    (Cekk, wL, cosmat, Pa, factorialmat, qPa, qfactorialmat, M_Sigmape, &
                    Pmlk, cmjk, phipu, SN_psip, SNM_sourcee) ! EXTREMELY INEFFICIENT. WHY?
            else if (FEM_type .eq. "T1") then
                print *, "WIP: SN to SN MGXS T1"
                stop
            else if (FEM_type .eq. "T2") then
                call ele_SN_to_SN_MGXS_particle_particle_scattering_source &
                    (Cekk, wL, cosmat, Pa, factorialmat, qPa, qfactorialmat, M_Sigmape, &
                    Pmlk, cmjk, phipu, e_SN_psip, e_SNM_sourcee)
            end if
        ! ### SN to PN MGXS ! HAS WIP (T1 and T2)
        else if (isMGXS .and. t_isSN .and. isPN) then
            if (FEM_type .eq. "default") then
                call SN_to_PN_MGXS_particle_particle_scattering_source &
                    (NL, Cekk, wL, RSH_unc, RSH, M_Sigmape, phipu, SN_psip, PNM_sourcee)
            else
                print *, "WIP: SN to PN MGXS T1/T2"
                stop
            end if
        ! ### PN to SN MGXS ! HAS WIP
        else if (isMGXS .and. t_isPN .and. isSN) then
            print *, "WIP: PN to SN MGXS"
            stop
        ! ### PN to PN MGXS ! HAS WIP
        else if (isMGXS .and. t_isPN .and. isPN) then
            print *, "WIP: PN to PN MGXS"
            stop
        ! ### SN to SN FEXS ! HAS WIP (T1 and T2)
        else if (isFEXS .and. t_isSN .and. isSN) then
            if (FEM_type .eq. "default") then
                call SN_to_SN_FEXS_particle_particle_scattering_source &
                    (wL, cosmat, Pa, factorialmat, qPa, qfactorialmat, F_Sigmape, &
                    Pmlk, cmjk, SN_psip, SNF_sourcee)
            else
                print *, "WIP: SN to SN FEXS T1/T2"
                stop
            end if
        ! ### SN to PN FEXS ! HAS WIP (T1 and T2)
        else if (isFEXS .and. t_isSN .and. isPN) then
            if (FEM_type .eq. "default") then
                call SN_to_PN_FEXS_particle_particle_scattering_source &
                    (NL, wL, RSH_unc, RSH, F_Sigmape, SN_psip, PNF_sourcee)
            else
                print *, "WIP: SN to PN FEXS T1/T2"
                stop
            end if
        ! ### PN to SN FEXS ! HAS WIP
        else if (isFEXS .and. t_isPN .and. isSN) then
            print *, "WIP: PN to SN FEXS"
            stop
        ! ### PN to PN FEXS ! HAS WIP
        else if (isFEXS .and. t_isPN .and. isPN) then
            print *, "WIP: PN to PN FEXS"
            stop
        end if
    end if
    ompfinish = omp_get_wtime()
    print *, "MODULE 8 - ELECTRONS:", real(ompfinish-ompstart), "seconds"

    !go to 99

    ! 9 - electrons
    ompstart = omp_get_wtime()
    isSI = electron_solution_method .eq. "SI"
    isGMRESm = .not. isSI
    ! ### SN MGXS SI
    if      (isMGXS .and. isSN .and. isSI) then
        if (FEM_type .eq. "default") then
            call SN_MGXS_EI_SI &
                (2, e_pmax, Ge, Cekk, bdyel, eprime, normal, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, MS_ee, MSt_e, phieu, &
                SNM_sourcee, M_fEg, bdysrc, SN_psie, phie)
        else if (FEM_type .eq. "T1") then
            call t1_ele_SN_MGXS_EI_SI &
                (2, e_pmax, Ge, Cekk, bdyel, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
                varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, &
                MS_ee, MSt_e, phieu, SNM_sourcee, M_fEg, bdysrc, SN_psie, phie)
        else if (FEM_type .eq. "T2") then
            call t2_ele_SN_MGXS_EI_SI &
                (2, e_pmax, Ge, Cekk, bdyel, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, MS_ee, MSt_e, phieu, &
                e_SNM_sourcee, M_fEg, bdysrc, e_SN_psie, phie)
        end if
    ! ### SN MGXS GMRESm
    else if (isMGXS .and. isSN .and. isGMRESm) then
        if (FEM_type .eq. "default") then
            call SN_MGXS_EI_GMRESm &
                (2, e_pmax, e_mmax, Ge, Cekk, bdyel, wL, khat, cosmat, &
                qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, MS_ee, MSt_e, &
                SNM_sourcee, M_fEg, bdysrc, SN_psie, phie)
        else if (FEM_type .eq. "T1") then
            call t1_ele_SN_MGXS_EI_GMRESm &
                (2, e_pmax, e_mmax, Ge, Cekk, bdyel, wL, khat, &
                cosmat, qPa, qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, I1invI2f, &
                I1invI3vector, MS_ee, MSt_e, SNM_sourcee, M_fEg, bdysrc, &
                SN_psie, phie)
        else if (FEM_type .eq. "T2") then
            call t2_ele_SN_MGXS_EI_GMRESm &
                (2, e_pmax, e_mmax, Ge, Cekk, bdyel, wL, khat, cosmat, qPa, &
                qfactorialmat, esweeplist, esweepbounds, enew, fold, &
                Ia, Ja, noJa, varsigmaup, varsigmadown, I1invI2, &
                I1invI2f, I1invI3vector, MS_ee, MSt_e, &
                e_SNM_sourcee, M_fEg, bdysrc, e_SN_psie, phie)
        end if
    ! ### PN MGXS GMRESm ! HAS WIP
    else if (isMGXS .and. isPN .and. isGMRESm) then
        print *, "WIP: PHOTON PN MGXS GMRESm"
    ! ### SN FEXS SI ! HAS WIP (T1 and T2)
    else if (isFEXS .and. isSN .and. isSI) then
        if (FEM_type .eq. "default") then
            call SN_FEXS_EI_SI &
                (2, e_pmax, Ge, Cekk, bdyel, wL, khat, cosmat, qPa, qfactorialmat, &
                esweeplist, esweepbounds, enew, fold, Ia, Ja, noJa, &
                varsigmaup, varsigmadown, I1invI2, I1invI2f, I1invI3vector, &
                F_Sigmaee, F_Sigmaet, SNF_sourcee, F_fEg, bdysrc, SN_psie, phie)
        else
            print *, "WIP: SN FEXS SI T1 and T2"
        end if
    ! ### SN FEXS GMRESm ! HAS WIP
    else if (isFEXS .and. isSN .and. isGMRESm) then
        print *, "WIP: PHOTON SN FEXS GMRESm"
    ! ### PN FEXS GMRESm ! HAS WIP
    else if (isFEXS .and. isPN .and. isGMRESm) then
        print *, "WIP: PHOTON PN FEXS GMRESm"
    end if
    go to 10

    9 continue
    print *, "NOTHING TO DO WITH MODULE 8 - ELECTRONS"
    outerompstart = omp_get_wtime()
    ! ### SN
    if (isSN) then
        if (FEM_type .eq. "default") then
            call SN_read_field("electron", merge(Ge, Ge+1, isMGXS), phieu, phie, SN_psie)
        else if (FEM_type .eq. "T2") then
            call t2_ele_SN_read_field("electron", merge(Ge, Ge+1, isMGXS), phieu, phie, e_SN_psie)
        end if
    ! ### PN
    else if (isPN) then
        if (FEM_type .eq. "default") then
            call PN_read_field("electron", merge(Ge, Ge+1, isMGXS), NL, phieu, phie, PN_psie)
        else if (FEM_type .eq. "T2") then
            print *, "WIP: PN T2"
            stop
        end if
    end if

    !!!!!!!!!!KEEP AND USE FOR PN!!!!!!!!!!
    !I1invI2 = reshape(I1invI2, [NE,NKe,NFe,NKe], order = [4,2,3,1]) ! (e,kp,f,k)
    !I1invI2f = reshape(I1invI2f, [NE,NKe,NFe,NKe], order = [2,4,3,1]) ! (e,kp,f,k)
    !I1invI3vector = reshape(I1invI3vector, [NE,NKe,3,NKe], order = [4,2,3,1]) ! (e,kp,dir,k)

    print *, "ELECTRON TRANSPORT READ."

    10 continue
    ompfinish = omp_get_wtime()
    print *, "MODULE 9 - ELECTRONS:", real(ompfinish-ompstart), "seconds"
    print *, "ELECTRON TRANSPORT COMPLETE."

    1 continue
end if

if (gen_post_proc) then
    ! POST-PROCESSING
    ompstart = omp_get_wtime()

    if (transport_mode .eq. "T MMS" .or. transport_mode .eq. "full MMS") then ! Remember not to remove the doslices thing.
        call construct_mesh_slice_overlay("YZ", mshplane, mshslice, rglobal, Cfkf, slicek, slicerglobal, sliceCfkf)

        if (MMS_angular .eq. "SN") then
            call SN_visualize_mesh_slice_angular_fluence &
                ("SN_MMS_exact_source", Cekk, slicek, SN_sourcem)
            call SN_visualize_mesh_slice_angular_fluence &
                ("SN_MMS_exact_solution", Cekk, slicek, SN_psim)
            call SN_visualize_mesh_slice_angular_fluence &
                ("SN_MMS_manual_source", Cekk, slicek, SN_sourcemm)
            call SN_visualize_mesh_slice_angular_fluence &
                ("SN_MMS_manual_solution", Cekk, slicek, SN_psimm)

            call SN_construct_fluence &
                (wL, SN_sourcem, g_phim)
            call visualize_mesh_slice_fluence &
                ("SN_MMS_exact_source", Cekk, slicek, g_phim)
            deallocate(g_phim)

            call SN_construct_fluence &
                (wL, SN_psim, g_phim)
            call visualize_mesh_slice_fluence &
                ("SN_MMS_exact_solution", Cekk, slicek, g_phim)
            deallocate(g_phim)

            call SN_construct_fluence &
                (wL, SN_sourcemm, g_phim)
            call visualize_mesh_slice_fluence &
                ("SN_MMS_manual_source", Cekk, slicek, g_phim)
            deallocate(g_phim)

            call SN_construct_fluence &
                (wL, SN_psimm, g_phim)
            call visualize_mesh_slice_fluence &
                ("SN_MMS_manual_solution", Cekk, slicek, g_phim)
            deallocate(g_phim)

            call SN_MGXS_polar_angle_distribution &
                ("MMS exact", Cekk, vol, mu, wL, [0.0], SN_psim)
            call SN_MGXS_polar_angle_distribution &
                ("MMS manual", Cekk, vol, mu, wL, [0.0], SN_psimm)

            print *, "SOLUTION NORM:"
            print *, norm2(SN_psim - SN_psimm)
        else if (MMS_angular .eq. "PN") then
            call PN_visualize_mesh_slice_angular_fluence &
                ("PN_MMS_exact_source", Cekk, slicek, PN_sourcem)
            call PN_visualize_mesh_slice_angular_fluence &
                ("PN_MMS_exact_solution", Cekk, slicek, PN_psim)
            call PN_visualize_mesh_slice_angular_fluence &
                ("PN_MMS_manual_source", Cekk, slicek, PN_sourcemm)
            call PN_visualize_mesh_slice_angular_fluence &
                ("PN_MMS_manual_solution", Cekk, slicek, PN_psimm)
        end if
    else
        if (allocated(phipu)) phip = phip + phipu
        if (allocated(noquad_phieu)) phie = phie + noquad_phieu

        if (doslices) &
            call construct_mesh_slice_overlay &
                ("YZ", mshplane, mshslice, rglobal, Cfkf, slicek, slicerglobal, sliceCfkf)

        if (solve_photons) then ! IF T2 ele, must homogenize psi too. Slow??
            allocate(g_phip(NK,merge(Gp,Gp+1,isMGXS)))
            if (allocated(phipu)) allocate(g_phipu(NK,merge(Gp,Gp+1,isMGXS)))
            do g = 1, Gp
                if (allocated(phipu)) call homogenize_DFEM(Cekk, phipu(:,:,g), g_phipu(:,g))
                call homogenize_DFEM(Cekk, phip(:,:,g), g_phip(:,g))
            end do

            if (doslices) then
                if (allocated(g_phipu)) call visualize_mesh_slice_fluence("photon_uncollided", Cekk, slicek, g_phipu)
                call visualize_mesh_slice_fluence("photon", Cekk, slicek, g_phip)
                call SN_visualize_mesh_slice_angular_fluence("photon", Cekk, slicek, SN_psip)
            end if

            if (isMGXS) then
                if (allocated(g_phipu)) call MGXS_energy_spectrum("photon", "uncollided", Cekk, vol, Ep, g_phipu)
                call MGXS_energy_spectrum("photon", "NULL", Cekk, vol, Ep, g_phip)

                call SN_MGXS_polar_angle_distribution("photon", Cekk, vol, mu, wL, Ep, SN_psip)
            else
                if (allocated(g_phipu)) call FEXS_energy_spectrum("photon", "uncollided", Cekk, vol, Ep, g_phipu)
                call FEXS_energy_spectrum("photon", "NULL", Cekk, vol, Ep, g_phip)

                call SN_FEXS_polar_angle_distribution("photon", Cekk, vol, mu, wL, Ep, SN_psip)
            end if

            if (doEDEP) then
                if (isMGXS) call DFEM_MGXS_calculate_deposition(Cekk, MEDEP_p, phip, e_ENERGY_p)
            end if
            if (doCDEP) then
                if (isMGXS) call DFEM_MGXS_calculate_deposition(Cekk, MCDEP_p, phip, e_CHARGE_p)
            end if
            if (doDDEP) then
                do mat = 1, Nmats
                    MEDEP_p(:,mat) = MEDEP_p(:,mat)/mat_rho(mat)
                end do
                if (isMGXS) call DFEM_MGXS_calculate_deposition(Cekk, MEDEP_p, phip, e_DOSE_p)
            end if
        end if

        if (solve_electrons) then
            allocate(g_phie(NK,merge(Ge,Ge+1,isMGXS)))
            if (allocated(phieu)) allocate(g_phieu(NK,merge(Ge,Ge+1,isMGXS)))
            do g = 1, Ge
                if (allocated(phieu)) call homogenize_DFEM(Cekk, phieu(:,:,g), g_phieu(:,g))
                call homogenize_DFEM(Cekk, phie(:,:,g), g_phie(:,g))
            end do

            if (doslices) then
                if (allocated(g_phieu)) call visualize_mesh_slice_fluence("electron_uncollided", Cekk, slicek, g_phieu) ! This is after beam-quad though
                call visualize_mesh_slice_fluence("electron", Cekk, slicek, g_phie)
                call SN_visualize_mesh_slice_angular_fluence("electron", Cekk, slicek, SN_psie)
            end if

            if (isMGXS) then
                if (allocated(g_phieu)) call MGXS_energy_spectrum("electron", "uncollided", Cekk, vol, Ee, g_phieu)
                call MGXS_energy_spectrum("electron", "NULL", Cekk, vol, Ee, g_phie)

                call SN_MGXS_polar_angle_distribution("electron", Cekk, vol, mu, wL, Ee, SN_psie)
            else
                if (allocated(g_phieu)) call FEXS_energy_spectrum("electron", "uncollided", Cekk, vol, Ee, g_phieu)
                call FEXS_energy_spectrum("electron", "NULL", Cekk, vol, Ee, g_phie)

                call SN_FEXS_polar_angle_distribution("electron", Cekk, vol, mu, wL, Ee, SN_psie)
            end if

            if (doEDEP) then
                if (isMGXS) call DFEM_MGXS_calculate_deposition(Cekk, MEDEP_e, phie, e_ENERGY_e)
            end if
            if (doCDEP) then
                if (isMGXS) call DFEM_MGXS_calculate_deposition(Cekk, MCDEP_e, phie, e_CHARGE_e)
            end if
            if (doDDEP) then
                do mat = 1, Nmats
                    MEDEP_e(:,mat) = MEDEP_e(:,mat)/mat_rho(mat)
                end do
                if (isMGXS) call DFEM_MGXS_calculate_deposition(Cekk, MEDEP_e, phie, e_DOSE_e)
            end if
        end if

        if (doEDEP) then
            allocate(e_ENERGY(NE,NKe))
            e_ENERGY = 0
            if (solve_photons) e_ENERGY = e_ENERGY + e_ENERGY_p
            if (solve_electrons) e_ENERGY = e_ENERGY + e_ENERGY_e

            allocate(ENERGY(NK))
            call homogenize_DFEM(Cekk, e_ENERGY, ENERGY)

            call save_deposition("energy", ENERGY)

            if (doslices) &
                call visualize_mesh_slice_deposition("energy", Cekk, slicek, ENERGY)

            if (doline) &
                call construct_line("energy deposition", [1,2], [0.0,0.0], rglobal, ENERGY)

            if (doplanarints) &
                call planar_integrals &
                    (r, eprime, vol, a, bvector, 1.0E-8, planarx0, planarxN, planarNP, e_ENERGY)

            call DFEM_total_deposition("energy", vol, e_ENERGY)
        end if
        if (doCDEP) then
            allocate(e_CHARGE(NE,NKe))
            e_CHARGE = 0
            if (solve_photons) e_CHARGE = e_CHARGE + e_CHARGE_p
            if (solve_electrons) e_CHARGE = e_CHARGE + e_CHARGE_e

            allocate(CHARGE(NK))
            call homogenize_DFEM(Cekk, e_CHARGE, CHARGE)

            call save_deposition("charge", CHARGE)

            if (doslices) &
                call visualize_mesh_slice_deposition("charge", Cekk, slicek, CHARGE)

            if (doline) &
                call construct_line("charge deposition", [1,2], [0.0,0.0], rglobal, CHARGE)

            if (doplanarints) &
                call planar_integrals &
                    (r, eprime, vol, a, bvector, 1.0E-8, planarx0, planarxN, planarNP, e_CHARGE)

            call DFEM_total_deposition("energy", vol, e_CHARGE)
        end if
        if (doDDEP) then
            allocate(e_DOSE(NE,NKe))
            e_DOSE = 0
            if (solve_photons) e_DOSE = e_DOSE + e_DOSE_p
            if (solve_electrons) e_DOSE = e_DOSE + e_DOSE_e

            allocate(DOSE(NK))
            call homogenize_DFEM(Cekk, e_DOSE, DOSE)

            call save_deposition("dose", DOSE)

            if (doslices) &
                call visualize_mesh_slice_deposition("dose", Cekk, slicek, DOSE)

            if (doline) &
                call construct_line("dose deposition", [1,2], [0.0,0.0], rglobal, DOSE)

            if (doplanarints) &
                call planar_integrals &
                    (r, eprime, vol, a, bvector, 1.0E-8, planarx0, planarxN, planarNP, e_DOSE) !!!!!! INVESTIGATE SURFACE NANS!!!!

            call DFEM_total_deposition("energy", vol, e_DOSE)
        end if
    end if

    ompfinish = omp_get_wtime()
    print *, "POST_PROCESSING:", real(ompfinish-ompstart), "seconds"
end if

!if (any(isnan(e_DOSE))) print *, "1"
!if (any(isnan(DOSE))) print *, "2"
!if (any(isnan(phie))) print *, "3"
!if (any(isnan(g_phie))) print *, "4"
!!if (any(isnan(phieu))) print *, "5"
!!if (any(isnan(g_phieu))) print *, "6"
!!if (any(isnan(noquad_phieu))) print *, "7"
!if (any(isnan(SN_psie))) print *, "8"

lastompfinish = omp_get_wtime()
print *, "TOTAL:", real(lastompfinish-firstompstart), "seconds"

print *, "COMPUTATION COMPLETED SUCCESSFULLY."

go to 10000

1100 continue
print *, "wiscoslab: WIP, NOT FOR USERS."

! 2
firstompstart = omp_get_wtime()
call wiscoslab_all_transport_preproc(rglobal, Cekk, r, vol)
ompfinish = omp_get_wtime()
print *, "MODULE 2:", real(ompfinish-firstompstart), "seconds"

! 3
print *, "NOTHING TO DO WITH MODULE 3"

! 4
ompstart = omp_get_wtime()
call SN_angular_quadrature_abscissae_and_weights(Nmu, Nphi, mu, wL, phiC)
call construct_discrete_ordinates(mu, phiC, khat)
call populate_real_spherical_harmonics(NL, mu, phiC, RSH)
if ((photon_angular .eq. "PN" .and. transport_mode .eq. "external photon beam") .or. &
    (electron_angular .eq. "PN" .and. transport_mode .eq. "external photon beam coupled") .or. &
    (electron_angular .eq. "PN" .and. transport_mode .eq. "external electron beam") .or. &
    (MMS_angular .eq. "PN" .and. transport_mode .eq. "T MMS") .or. &
    (MMS_angular .eq. "PN" .and. transport_mode .eq. "full MMS")) then
    print *, "PN wiscoslab WIP"
    stop
end if
if (photon_angular .eq. "SN" .or. "electron_angular" .eq. "SN" .or. MMS_angular .eq. "SN") then
    call SN_scattering_angular_quadrature_terms &
        (mu, phiC, cosmat, Pa, factorialmat, qPa, qfactorialmat)
    call SN_uncollided_angular_quadrature_terms(rglobal, phiC, Pmlk, cmjk)
end if
NQ = (NL+1)**2
ompfinish = omp_get_wtime()
print *, "MODULE 4:", real(ompfinish-ompstart), "seconds"

! 5
print *, "NOTHING TO DO WITH MODULE 5"

! 6
ompstart = omp_get_wtime()
call wiscoslab_inner_products(vol, mu, GpFterm, Fsweep)
ompfinish = omp_get_wtime()
print *, "MODULE 6:", real(ompfinish-ompstart), "seconds"

print *, "PRE-PROCESSING COMPLETE."

if (energy_discretization_method .eq. "FEXS" .or. ccut .eq. 0) then
    XSlogic = .true.
else if (energy_discretization_method .eq. "MGXS") then
    XSlogic = .false.
end if

print *, "ELECTRON TRANSPORT STARTED."

! 7 - electrons
ompstart = omp_get_wtime()
call construct_exceptions("electron", electron_exc)
if (electron_energy_structure .eq. "logarithmic") then
    call logarithmic_energy_groups(Ge, Eemin, Eemax, electron_exc, Ee, dEe, Eemid)
else if (electron_energy_structure .eq. "linear") then
    call linear_energy_groups(Ge, Eemin, Eemax, electron_exc, Ee, dEe, Eemid)
else if (electron_energy_structure .eq. "exponential") then
    call exponential_energy_groups(Ge, Eemin, Eemax, electron_exc, Ee, dEe, Eemid)
end if
if (transport_mode .eq. "external electron beam") then
    if (beam_energy_dist .eq. "polychromatic") then
        call uncollided_beam_energy_distribution(fEspec, fE)
    end if
end if
ompfinish = omp_get_wtime()
print *, "MODULE 7 - ELECTRONS:", real(ompfinish-ompstart), "seconds"

if (energy_discretization_method .eq. "MGXS") then
    ! MGXS - electrons
    ompstart = omp_get_wtime()
    if (transport_mode .eq. "external photon beam coupled") then
        call MGXS_photon_electron &
            (Ep, dEp, Ee, dEe, Eemid, .true., 1.0E-8, 64, M_Sigmape, MEDEP_p)
    end if
    call MGXS_electron_electron &
        (Ee, dEe, Eemid, .true., 1.0E-8, 64, MS_ee, MSt_e, MSa_e, MS2_e, MEDEP_e, MCDEP_e)
    ompfinish = omp_get_wtime()
    print *, "MODULE MGXS - ELECTRONS:", real(ompfinish-ompstart), "seconds"
else if (energy_discretization_method .eq. "FEXS") then
    ! FEXS - electrons
    ompstart = omp_get_wtime()
    if (transport_mode .eq. "external photon beam coupled") then
        call FEXS_photon_electron(Ep, dEp, Ee, dEe, Eemid, .true., 1.0E-8, 64, F_Sigmape)
    end if
    call FEXS_electron_electron(Ee, dEe, Eemid, .false., 1.0E-2, 16, F_Sigmaee, Fp_Sigmaet, pi_Sigmaet)
    call FEXS_electron_attn(Ee, dEe, Fp_Sigmaet, pi_Sigmaet, .false., 1.0E-4, 16, F_Sigmaet, i_Sigmaet)
    ompfinish = omp_get_wtime()
    print *, "MODULE FEXS - ELECTRONS:", real(ompfinish-ompstart), "seconds"
end if

! 8 - electrons
ompstart = omp_get_wtime()
if (energy_discretization_method .eq. "MGXS") then
    if (FCS) then
        call wiscoslab_MGXS_uncollided_fluence &
            ("electron", rglobal, Ee, fEspec, fE, MSt_e, phieu)
        !allocate(noquad_phieu, source = phieu)
        call wiscoslab_MGXS_beam_quadrature(rglobal, vol, MSt_e, 4, phieu)
        !allocate(onlyquad_phieu, source = phieu)
        !go to 12123
        if (inelastic_ETC .or. exact_RCSDA_angular) then
            call wiscoslab_SN_MGXS_ETC_uncollided_fluence &
                (rglobal, mu, wL, GpFterm, Fsweep, MS_ee(:,:,NL+2,:), MSt_e, phieu)
            !noquad_phieu = noquad_phieu + phieu - onlyquad_phieu
        end if
        call wiscoslab_SN_MGXS_uncollided_source(mu, MS_ee, phieu, SNM_sourcee)
    else
        call wiscoslab_SN_MGXS_boundary_source_term(mu, Fsweep, Ee, fEspec, fE, bdy, M_fEg)
    end if
else
    call wiscoslab_SN_FEXS_boundary_source_term(mu, Fsweep, Ee, fEspec, fE, bdy, F_fEg)
end if
ompfinish = omp_get_wtime()
print *, "MODULE 8 - ELECTRONS:", real(ompfinish-ompstart), "seconds"

! 9 - electrons
ompstart = omp_get_wtime()
if (energy_discretization_method .eq. "MGXS") then
    call wiscoslab_SN_MGXS_EI_GMRESm &
        (2, e_pmax, e_mmax, Ge, mu, wL, GpFterm, Fsweep, MS_ee, MSt_e, &
        SNM_sourcee, bdy, M_fEg, SN_psie, phie)

    if (allocated(bdy)) deallocate(bdy)
else
    call wiscoslab_SN_FEXS_EI_GMRESm &
        (2, e_pmax, e_mmax, Ge, mu, wL, GpFterm, Fsweep, F_Sigmaee, F_Sigmaet, &
        SNF_sourcee, bdy, F_fEg, e_SN_psie, SN_psie)
    if (allocated(bdy)) deallocate(bdy)
end if
ompfinish = omp_get_wtime()
print *, "MODULE 9 - ELECTRONS:", real(ompfinish-ompstart), "seconds"
print *, "ELECTRON TRANSPORT COMPLETE."

if (energy_discretization_method .eq. "MGXS") then
    if (FCS) phie = phie + phieu ! I've turned off beam quadrature for the moment !!! 1D beam quad has issues right now!!!
    call DFEM_MGXS_calculate_deposition(Cekk, MEDEP_e, phie, e_ENERGY)
    call DFEM_MGXS_calculate_deposition(Cekk, MCDEP_e, phie, e_CHARGE)

    do mat = 1, Nmats
        MEDEP_e(:,mat) = MEDEP_e(:,mat)/mat_rho(mat)
    end do
    call DFEM_MGXS_calculate_deposition(Cekk, MEDEP_e, phie, e_DOSE)
else
    call FEXS_stopping_power(Ee, dEe, F_Scol)
    call DFEM_FEXS_calculate_deposition(Cekk, F_Scol, SN_psie, e_DOSE) ! USING SN_psi AS DUMMY ARRAY for fluence
end if

!! ### DOSE, CHARGE, SPACE
!do e = 1, NE
!    do k = 1, NKe
!        print *, r(1,k,e)
!    end do
!end do
!print *, "BREAK"
!do e = 1, NE
!    do k = 1, NKe
!        print *, e_DOSE(e,k)
!    end do
!end do
!print *, "BREAK"
!do e = 1, NE
!    do k = 1, NKe
!        print *, e_CHARGE(e,k)
!    end do
!end do
!print *, "BREAK"
!do e = 1, NE
!    do k = 1, NKe
!        print *, sum(phie(e,k,:))
!    end do
!end do

!! ### ANGLE
!do i = 1, Nmu
!    print *, mu(i)
!end do
!print *, "BREAK" ! Does not include phiu. Means its a little biased across FCS and non-FCS
!do i = 1, Nmu
!    S1 = 0
!    do g = 1, Ge
!        S1 = S1 + sum(matmul(transpose(SN_psie(1:NE,1:NKe,i,g)),vol))/NKe
!    end do
!    print *, S1
!end do

!! ### ENERGY
!do g = 1, Ge
!    do n = 1, 2
!        print *, Ee(g+n-1)
!    end do
!end do
!print *, "BREAK" ! Separate collided and uncollided?? Would not be good for cross-comparing FCS and non-FCS
!do g = 1, Ge
!    do n = 1,2
!        print *, sum(matmul(transpose(phie(1:NE,1:NKe,g)),vol))/NKe
!    end do
!end do

allocate(DOSE(NK))
DOSE = 0
do k = 1, NKe
    do e = 1, NE
        DOSE(e+k-1) = DOSE(e+k-1) + &
            e_DOSE(e,k)/merge(2, 1, e+k-1 .ne. 1 .and. e+k-1 .ne. NK)
    end do
end do

do k = 1, NK
    print *, rglobal(1,k)
end do
print *, "BREAK"
do k = 1, NK
    print *, DOSE(k)
end do

!call DFEM_total_deposition("energy", vol, e_ENERGY)
!call DFEM_total_deposition("charge", vol, e_CHARGE)
!
!do mat = 1, Nmats
!    MEDEP_e(:,mat) = MEDEP_e(:,mat)*mat_rho(mat)
!end do
!call wiscoslab_DFEM_SN_MGXS_check_energy_conservation &
!    (vol, mu, wL, Eemid, MEDEP_e, M_fEg, phieu, phie, SN_psie)

!12123 continue
!allocate(DOSE(NK))
!DOSE = 0
!do k = 1, NKe
!    do e = 1, NE
!        DOSE(e+k-1) = DOSE(e+k-1) + &
!            sum(phieu(e,k,:))/merge(2, 1, e+k-1 .ne. 1 .and. e+k-1 .ne. NK)
!    end do
!end do
!do k = 1, NK
!    print *, DOSE(k)
!end do

10000 continue

!123 continue

!call DFEM_construct_line("numerical phiu alone", [1,2], [0.0,0.0], r, sum(noquad_phieu,3))
!!call DFEM_construct_line("onlyquad_phieu", [1,2], [0.0,0.0], r, sum(onlyquad_phieu,3))
!call DFEM_construct_line("full", [1,2], [0.0,0.0], r, sum(phieu,3))

end program solver