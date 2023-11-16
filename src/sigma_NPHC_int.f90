subroutine sigma_NPHC_int
    !> Calculate the intrinsic nonlinear planar Hall conductivity, the xyyy and yxxx elements
    !
    !> usage: sigma_NPHC_int_calc = T
    !
    !> ref : 10.1103/PhysRevLett.130.126303
    !
    !> 2023/10/31 Fan Yang
    !
    !> The Lande g-factor for spin degree is fixed to 2
    !> The input hr.dat must be from wannier90 v2.0 or newer versions

    use wmpi
    use para
    use nonlinear_transport
    use magnetic_moments
    implicit none

    integer :: NumberofEta              !> NumT

    integer :: ik, ikx, iky, ikz
    integer :: ie, ieta
    integer :: ierr, knv3
    real(dp) :: k(3)

    character*40 :: ahcfilename, etaname
    real(dp), allocatable :: Eta_array(:)

    real(dp), allocatable :: energy(:)!> Fermi energy, dim= OmegaNum

    real(dp) :: time_start, time_end

    complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
    complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments

    real(dp), allocatable :: Chi_xyyy_tensor    (:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi(:,:)
    real(dp), allocatable :: Chi_yxxx_tensor    (:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi(:,:)

    allocate( energy (OmegaNum))

    allocate( M_S(Num_wann, Num_wann,3) )
    allocate( M_L(Num_wann, Num_wann,3) )

    !> temperature, Eta = 1/k_B/T
    NumberofEta = NumT
    allocate( Eta_array(NumberofEta) )
    Eta_array(1) = Eta_Arc !> from wt.in
    if ((NumberofEta>1) .and. (NumberofEta<8)) then ! 1+6=7
        Eta_array(2:NumberofEta) = Eta_array_all(1:NumberofEta-1)
    else if (NumberofEta>7) then
        stop "The NumT should not more than 7"
    endif

    allocate( Chi_xyyy_tensor    (OmegaNum, NumberofEta))
    allocate( Chi_xyyy_tensor_mpi(OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor    (OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor_mpi(OmegaNum, NumberofEta))

    Chi_xyyy_tensor     = 0d0
    Chi_xyyy_tensor_mpi = 0d0
    Chi_yxxx_tensor     = 0d0
    Chi_yxxx_tensor_mpi = 0d0

    !> Fermi energy in Hatree energy, not eV
    do ie=1, OmegaNum
        if (OmegaNum>1) then
            energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie- 1d0)/dble(OmegaNum- 1)
        else
            energy= OmegaMin
        endif
    enddo ! ie

    !> difference on magnetic field
    if (cpuid .eq. 0) then
        if ((Bx==0d0) .or. (By==0d0)) then
            write(stdout, '("We did not find Bx or By in wt.in, so we use the default value Bx=By=2 Tesla")')
            Bx = 2d0
            By = 2d0
            Bz = 2d0
        endif
    endif

    !> generate Pauli matrix
    call spin_magnetic_moments(M_S)
    M_S = M_S * Lande_g_S * mu_B * 0.5d0
    M_S(:,:,1) = M_S(:,:,1) * Bx
    M_S(:,:,2) = M_S(:,:,2) * By
    M_S(:,:,3) = M_S(:,:,3) * Bz

    knv3= Nk1*Nk2*Nk3

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 1000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/1000d0/60d0
            time_start= time_end
        endif

        ikx= (ik-1)/(Nk2*Nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube  &
            + K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(Nk3)

        call sigma_NPHC_int_single_k(M_S, &
            k, energy, NumberofEta, Eta_array, Chi_xyyy_k, Chi_yxxx_k)

        Chi_xyyy_tensor_mpi = Chi_xyyy_tensor_mpi+ Chi_xyyy_k
        Chi_yxxx_tensor_mpi = Chi_yxxx_tensor_mpi+ Chi_yxxx_k
    enddo ! ik

#if defined (MPI)
    call mpi_allreduce(Chi_xyyy_tensor_mpi,Chi_xyyy_tensor,size(Chi_xyyy_tensor),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Chi_yxxx_tensor_mpi,Chi_yxxx_tensor,size(Chi_yxxx_tensor),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
    Chi_xyyy_tensor= Chi_xyyy_tensor_mpi
    Chi_yxxx_tensor= Chi_yxxx_tensor_mpi
#endif
    Chi_xyyy_tensor= Chi_xyyy_tensor * SOAHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    Chi_yxxx_tensor= Chi_yxxx_tensor * SOAHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, NumberofEta
            write(etaname, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree
            write(ahcfilename, '(7a)')'sigma_NPHC_int_eta', trim(adjustl(etaname)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1'
            write(outfileindex, '("#",a)')' For 2D cases, you need to multiply the 3rd lattice vector in SI unit'
            write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy', '\sigma_yxxx'
            do ie=1, OmegaNum
                write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, Chi_xyyy_tensor(ie,ieta), &
                    Chi_yxxx_tensor(ie,ieta)
            enddo
            close(outfileindex)
        enddo
    endif

    deallocate(energy, Eta_array)
    deallocate(M_S, M_L)
    deallocate(Chi_xyyy_tensor, Chi_xyyy_tensor_mpi, Chi_yxxx_tensor, Chi_yxxx_tensor_mpi)

    return

end subroutine sigma_NPHC_int


subroutine sigma_NPHC_int_single_k(M_S, &
    k, energy, NumberofEta, Eta_array, Chi_xyyy_k, Chi_yxxx_k)

    use nonlinear_transport
    use para
    implicit none

    complex(dp),intent(in):: M_S(Num_wann, Num_wann,3) 

    real(dp), intent(in)  :: k(3)
    real(dp), intent(in)  :: energy(OmegaNum)
    integer , intent(in)  :: NumberofEta
    real(dp), intent(in)  :: Eta_array(NumberofEta)
    real(dp), intent(out) :: Chi_xyyy_k(OmegaNum, NumberofEta)
    real(dp), intent(out) :: Chi_yxxx_k(OmegaNum, NumberofEta)

    ! eigen value of H
    real(dp),    allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk   (:, :)
    complex(dp), allocatable :: Hamk_bulk_dx(:, :)
    complex(dp), allocatable :: Hamk_bulk_dy(:, :)
    complex(dp), allocatable :: UU(:, :)
    !> velocities
    complex(dp), allocatable :: vx(:, :), vy(:, :)

    real(dp), allocatable :: sigma_xyy_k   (:,:)
    real(dp), allocatable :: sigma_yxx_k   (:,:)
    real(dp), allocatable :: sigma_xyy_k_dx(:,:)
    real(dp), allocatable :: sigma_yxx_k_dx(:,:)
    real(dp), allocatable :: sigma_xyy_k_dy(:,:)
    real(dp), allocatable :: sigma_yxx_k_dy(:,:)

    allocate( sigma_xyy_k    (OmegaNum, NumberofEta))
    allocate( sigma_yxx_k    (OmegaNum, NumberofEta))
    allocate( sigma_xyy_k_dx (OmegaNum, NumberofEta))
    allocate( sigma_yxx_k_dx (OmegaNum, NumberofEta))
    allocate( sigma_xyy_k_dy (OmegaNum, NumberofEta))
    allocate( sigma_yxx_k_dy (OmegaNum, NumberofEta))

    allocate( W(Num_wann))
    allocate( Hamk_bulk   (Num_wann, Num_wann))
    allocate( Hamk_bulk_dx(Num_wann, Num_wann))
    allocate( Hamk_bulk_dy(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))

    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))

    Hamk_bulk= 0d0
    UU= 0d0

    !> ground state, B=0 ---------------------------------------------------------
    call ham_bulk_latticegauge(k, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    call velocity_latticegauge_simple(k, UU, vx, vy)
    call Lambda_abc_sum(energy, W, vx, vy, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)

    Hamk_bulk_dx = Hamk_bulk
    Hamk_bulk_dy = Hamk_bulk
    if (include_m_spin) then
        Hamk_bulk_dx = Hamk_bulk_dx + M_S(:,:,1)
        Hamk_bulk_dy = Hamk_bulk_dy + M_S(:,:,2)
    endif
    if (include_m_orb) then
        !Hamk_bulk_dx = Hamk_bulk_dx + mx_L
        !Hamk_bulk_dy = Hamk_bulk_dy + my_L
    endif

    UU=Hamk_bulk_dx
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    call velocity_latticegauge_simple(k, UU, vx, vy)
    call Lambda_abc_sum(energy, W, vx, vy, NumberofEta, Eta_array, sigma_xyy_k_dx, sigma_yxx_k_dx)

    UU=Hamk_bulk_dy
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    call velocity_latticegauge_simple(k, UU, vx, vy)
    call Lambda_abc_sum(energy, W, vx, vy, NumberofEta, Eta_array, sigma_xyy_k_dy, sigma_yxx_k_dy)

    Chi_xyyy_k = Chi_xyyy_k + (sigma_xyy_k_dy - sigma_xyy_k)/By
    Chi_yxxx_k = Chi_yxxx_k + (sigma_yxx_k_dx - sigma_yxx_k)/Bx

    deallocate( W, Hamk_bulk, Hamk_bulk_dx, Hamk_bulk_dy, vx, vy, UU)
    return
end subroutine sigma_NPHC_int_single_k

