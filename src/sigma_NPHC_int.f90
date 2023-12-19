subroutine sigma_NPHC_int
    !> Calculate the intrinsic nonlinear planar Hall conductivity, the xyyy and yxxx elements
    !
    !> usage: sigma_NPHC_int_calc = T
    !
    !> ref : 10.1103/PhysRevLett.130.126303
    !
    !> 2023/10/31 Fan Yang
    !

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    !> magnetic moments in nonlinear planar Hall
    logical               :: include_m_spin = .false.
    logical               :: include_m_orb  = .true.

    integer :: NumberofEta              !> NumT

    integer :: ik, ikx, iky, ikz
    integer :: ie, ieta
    integer :: ierr, knv3
    real(dp) :: k(3)

    character*40 :: ahcfilename, etaname
    real(dp), allocatable :: Eta_array(:)

    real(dp), allocatable :: energy(:)!> Fermi energy, dim= OmegaNum

    real(dp) :: time_start, time_end

    real(dp), allocatable :: Chi_xyyy_k         (:,:)
    real(dp), allocatable :: Chi_xyyy_tensor    (:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi(:,:)
    real(dp), allocatable :: Chi_yxxx_k         (:,:)
    real(dp), allocatable :: Chi_yxxx_tensor    (:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi(:,:)

    allocate( energy (OmegaNum))

    !> temperature, Eta = 1/k_B/T
    NumberofEta = NumT
    allocate( Eta_array(NumberofEta) )
    Eta_array(1) = Eta_Arc !> from wt.in
    if ((NumberofEta>1) .and. (NumberofEta<8)) then ! 1+6=7
        Eta_array(2:NumberofEta) = Eta_array_all(1:NumberofEta-1)
    else if (NumberofEta>7) then
        stop "The NumT should not more than 7"
    endif

    allocate( Chi_xyyy_k    (OmegaNum, NumberofEta))
    allocate( Chi_yxxx_k    (OmegaNum, NumberofEta))
    allocate( Chi_xyyy_tensor    (OmegaNum, NumberofEta))
    allocate( Chi_xyyy_tensor_mpi(OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor    (OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor_mpi(OmegaNum, NumberofEta))

    Chi_xyyy_tensor     = 0d0
    Chi_xyyy_tensor_mpi = 0d0
    Chi_yxxx_tensor     = 0d0
    Chi_yxxx_tensor_mpi = 0d0

    if (cpuid .eq. 0) then
        if (include_m_spin) then
            write(stdout, '("You have considered the spin    magnetic moments.")')
        endif
        if (include_m_orb) then
            write(stdout, '("You have considered the orbital magnetic moments.")')
        endif
    endif

    !> Fermi energy in Hatree energy, not eV
    do ie=1, OmegaNum
        if (OmegaNum>1) then
            energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie- 1d0)/dble(OmegaNum- 1)
        else
            energy= OmegaMin
        endif
    enddo ! ie

    knv3= Nk1*Nk2*Nk3

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if ( cpuid.eq.0 .and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif

        if ( (real(ik)/real(knv3))>0.95 .and. mod( (ik-cpuid)/num_cpu, 2000 ).eq.0) then
            call now(time_end)
            write(*     , '(a, i5, a, f10.2, "min")') 'cpuid =', cpuid, &
                '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif

        ikx= (ik-1)/(Nk2*Nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube  &
            + K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(Nk3)

        call sigma_INPHC_single_k(include_m_spin, include_m_orb, k, energy, NumberofEta, Eta_array, Chi_xyyy_k, Chi_yxxx_k)

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
    Chi_xyyy_tensor= Chi_xyyy_tensor * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    Chi_yxxx_tensor= Chi_yxxx_tensor * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, NumberofEta
            write(etaname, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree
            if ((include_m_spin) .and. (include_m_orb)) then
                write(ahcfilename, '(7a)')'sigma_INPHC_LS_eta', trim(adjustl(etaname)), 'meV.dat'
            else if (include_m_spin) then
                write(ahcfilename, '(7a)')'sigma_INPHC_S_eta', trim(adjustl(etaname)), 'meV.dat'
            else if (include_m_orb ) then
                write(ahcfilename, '(7a)')'sigma_INPHC_L_eta', trim(adjustl(etaname)), 'meV.dat'
            else
                write(ahcfilename, '(7a)')'sigma_INPHC_eta', trim(adjustl(etaname)), 'meV.dat'
            endif

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
    deallocate(Chi_xyyy_k, Chi_yxxx_k, Chi_xyyy_tensor, Chi_xyyy_tensor_mpi, Chi_yxxx_tensor, Chi_yxxx_tensor_mpi)

    return

end subroutine


subroutine sigma_INPHC_single_k(include_m_spin, include_m_orb, k, energy, NumberofEta, Eta_array, Chi_xyyy_k, Chi_yxxx_k)

    use nonlinear_transport, only: velocity_latticegauge_simple
    use magnetic_moments
    use para
    implicit none

    logical, intent(in)   :: include_m_spin
    logical, intent(in)   :: include_m_orb 

    real(dp), intent(in)  :: k(3)
    real(dp), intent(in)  :: energy(OmegaNum)
    integer , intent(in)  :: NumberofEta
    real(dp), intent(in)  :: Eta_array(NumberofEta)
    real(dp), intent(out) :: Chi_xyyy_k(OmegaNum, NumberofEta)
    real(dp), intent(out) :: Chi_yxxx_k(OmegaNum, NumberofEta)

    complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
    complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments

    integer :: n, m, l, ie, ieta
    
    real(dp) :: k_dx(3)
    real(dp) :: k_dy(3)

    !> Fermi-Dirac distribution
    real(dp) :: mu, diffFermi
    
    ! eigen value of H
    real(dp),    allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

    real(dp), allocatable :: W_dx(:)
    real(dp), allocatable :: W_dy(:)

    complex(dp), allocatable :: sx(:, :), sy(:, :), sz(:, :)
    complex(dp), allocatable :: vx(:, :), vy(:, :)
    complex(dp), allocatable :: velocities(:,:,:)

    complex(dp), allocatable :: vx_dx(:, :), vy_dx(:, :)  
    complex(dp), allocatable :: vx_dy(:, :), vy_dy(:, :) 

    real(dp) :: Lambda_xyy, Lambda_yyy, Lambda_yxx, Lambda_xxx
    real(dp) :: G_xx, G_xy, G_yx, G_yy, G_yy_dx, G_xy_dy, G_xx_dy, G_yx_dx
    real(dp) :: dEnm, dEnm3, dEml, dEnl
    
    !===========================================================================
    !> original kpoints
    allocate( W (Num_wann))
    allocate( Hamk_bulk (Num_wann, Num_wann))
    allocate( Amat (Num_wann, Num_wann))
    allocate( UU (Num_wann, Num_wann))
    allocate( UU_dag (Num_wann, Num_wann))

    allocate( velocities(Num_wann, Num_wann, 3))
    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
    allocate( sx(Num_wann, Num_wann), sy(Num_wann, Num_wann), sz(Num_wann, Num_wann))

    call ham_bulk_latticegauge(k, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    UU_dag= conjg(transpose(UU))
    call velocity_latticegauge_simple(k, UU, velocities)
    vx = velocities(:,:,1)
    vy = velocities(:,:,2)

    !===========================================================================
    !> magnetic operators 
    allocate( M_S(Num_wann, Num_wann,3) )
    allocate( M_L(Num_wann, Num_wann,3) )

    M_S = 0d0
    M_L = 0d0
    if (include_m_spin) then
        call spin_magnetic_moments(M_S)

        call mat_mul(Num_wann, M_S(:,:,1), UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, M_S(:,:,1)) 
        call mat_mul(Num_wann, M_S(:,:,2), UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, M_S(:,:,2))
        ! call mat_mul(Num_wann, M_S(:,:,3), UU, Amat)
        ! call mat_mul(Num_wann, UU_dag, Amat, M_S(:,:,3))
    endif
    if (include_m_orb) then
        call orbital_magnetic_moments(W, velocities, M_L)
    endif

    sx = -0.5d0 * Lande_g_S * M_S(:,:,1) + Lande_g_L * M_L(:,:,1)
    sy = -0.5d0 * Lande_g_S * M_S(:,:,2) + Lande_g_L * M_L(:,:,2)
    sz = -0.5d0 * Lande_g_S * M_S(:,:,3) + Lande_g_L * M_L(:,:,3) 


    !===========================================================================
    !> k + dk_x
    allocate( W_dx (Num_wann))   
    allocate( vx_dx(Num_wann, Num_wann), vy_dx(Num_wann, Num_wann))

    k_dx = k+(/Origin_cell%Rua(1)*dx , Origin_cell%Rub(1)*dx , Origin_cell%Ruc(1)*dx/)/twopi

    call ham_bulk_latticegauge(k_dx, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W_dx)
    UU_dag= conjg(transpose(UU))
    call velocity_latticegauge_simple(k_dx, UU, velocities)
    vx_dx = velocities(:,:,1)
    vy_dx = velocities(:,:,2)
    !===========================================================================

    !===========================================================================
    !> k + dk_y
    allocate( W_dy (Num_wann))
    allocate( vx_dy(Num_wann, Num_wann), vy_dy(Num_wann, Num_wann))

    k_dy = k+(/Origin_cell%Rua(2)*dy , Origin_cell%Rub(2)*dy , Origin_cell%Ruc(2)*dy/)/twopi

    call ham_bulk_latticegauge(k_dy, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W_dy)
    UU_dag= conjg(transpose(UU))
    call velocity_latticegauge_simple(k_dy, UU, velocities)
    vx_dy = velocities(:,:,1)
    vy_dy = velocities(:,:,2)
    !===========================================================================

    Chi_xyyy_k = 0d0
    Chi_yxxx_k = 0d0

    do n= 1, Num_wann
        if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle 
        G_xx= 0d0
        G_xy= 0d0
        G_yx= 0d0
        G_yy= 0d0
        G_yy_dx= 0d0 
        G_xy_dy= 0d0 
        G_xx_dy= 0d0 
        G_yx_dx= 0d0        
        Lambda_xyy= 0d0
        Lambda_yyy= 0d0
        Lambda_yxx= 0d0
        Lambda_xxx= 0d0

        do m= 1, Num_wann
            dEnm= W(n) - W(m)           
            if (ABS(dEnm) < band_degeneracy_threshold) cycle

            dEnm3= dEnm**3
            G_xx= G_xx+ 2.d0*real( vx(n, m)*vx(m, n) )/dEnm3
            G_xy= G_xy+ 2.d0*real( vx(n, m)*vy(m, n) )/dEnm3
            G_yx= G_yx+ 2.d0*real( vy(n, m)*vx(m, n) )/dEnm3
            G_yy= G_yy+ 2.d0*real( vy(n, m)*vy(m, n) )/dEnm3

            G_yy_dx= G_yy_dx + 2.d0*real( vy_dx(n, m)*vy_dx(m, n) )/(W_dx(n) - W_dx(m))**3
            G_yx_dx= G_yx_dx + 2.d0*real( vy_dx(n, m)*vx_dx(m, n) )/(W_dx(n) - W_dx(m))**3
            
            G_xy_dy= G_xy_dy + 2.d0*real( vx_dy(n, m)*vy_dy(m, n) )/(W_dy(n) - W_dy(m))**3
            G_xx_dy= G_xx_dy + 2.d0*real( vx_dy(n, m)*vx_dy(m, n) )/(W_dy(n) - W_dy(m))**3

            Lambda_xyy= Lambda_xyy + 6.d0* real( vx(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
            Lambda_yyy= Lambda_yyy + 6.d0* real( vy(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
            Lambda_yxx= Lambda_yxx + 6.d0* real( vy(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
            Lambda_xxx= Lambda_xxx + 6.d0* real( vx(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
            
            do l= 1, Num_wann
                dEnl= W(n)-W(l)
                dEml= W(m)-W(l)                    
                if (ABS(dEnl) > band_degeneracy_threshold) then
                    Lambda_xyy= Lambda_xyy - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *sy(n, l)) /dEnm3/dEnl
                    Lambda_yyy= Lambda_yyy - 2.d0* real( (vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)) *sy(n, l)) /dEnm3/dEnl
                    Lambda_yxx= Lambda_yxx - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *sx(n, l)) /dEnm3/dEnl
                    Lambda_xxx= Lambda_xxx - 2.d0* real( (vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)) *sx(n, l)) /dEnm3/dEnl
                endif
                if (ABS(dEml) > band_degeneracy_threshold) then
                    Lambda_xyy= Lambda_xyy - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *sy(m, l)) /dEnm3/dEml
                    Lambda_yyy= Lambda_yyy - 2.d0* real( (vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)) *sy(m, l)) /dEnm3/dEml
                    Lambda_yxx= Lambda_yxx - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *sx(m, l)) /dEnm3/dEml
                    Lambda_xxx= Lambda_xxx - 2.d0* real( (vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)) *sx(m, l)) /dEnm3/dEml
                endif
            enddo ! l

        enddo ! m

        do ieta=1, NumberofEta
            do ie=1, OmegaNum
                mu = energy(ie)

                if ((W(n)-mu)/Eta_array(ieta)<50) then
                    diffFermi= -Exp((W(n)-mu)/Eta_array(ieta))/(Exp((W(n)-mu)/Eta_array(ieta))+1d0)**2 /Eta_array(ieta)
                else
                    diffFermi=0.d0
                endif

                Chi_xyyy_k(ie,ieta) = Chi_xyyy_k(ie,ieta) + real( (vx(n,n)*Lambda_yyy - vy(n,n)*Lambda_xyy) &
                    + ((G_yy_dx - G_yy)/dx - (G_xy_dy - G_xy)/dy)*sy(n,n) ) * diffFermi
                Chi_yxxx_k(ie,ieta) = Chi_yxxx_k(ie,ieta) + real( (vy(n,n)*Lambda_xxx - vx(n,n)*Lambda_yxx) &
                    + ((G_xx_dy - G_xx)/dy - (G_yx_dx - G_yx)/dx)*sx(n,n) ) * diffFermi              

            enddo ! ie
        enddo ! ieta
    enddo ! n

    deallocate(W, vx, vy, Hamk_bulk, Amat, UU, UU_dag, velocities)
    deallocate(W_dx, W_dy, vx_dx, vx_dy, vy_dx, vy_dy)
    deallocate(sx, sy, sz)
    deallocate(M_S, M_L)
    return

end subroutine sigma_INPHC_single_k

