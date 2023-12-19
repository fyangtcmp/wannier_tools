module nonlinear_transport
    use para, only: dp, eV2Hartree, Echarge, mu_B, Hartree2J, hbar, Num_wann, OmegaNum, zi
    implicit none

    !> adaptive k-meshes method
    logical :: use_adaptive_method = .false.
    real(dp):: adaptive_threshold  = 1e-10
    integer :: Nk_local = 6          !> local k-points along each directions of the dense k-meshes

    ! The energy threshold to judge the degenerate bands
    ! The intrinsic SOAHC is usually calculated in PT-symmetric materials, with double degenerate bands.
    ! And for normal hr.dat, the numerical error between two physical degenerate bands is about 0.1~ meV.
    ! So we set it to 1meV = 3.6749d-5 Hartree, a relatively large value. It will lead to smoother curves of the conductivities.
    ! If you want larger peaks of the conductivities, or you have symmetrized your hr.dat to reduce the numerical error,
    ! you can set the threshold to a smaller value like 0.01meV = 3.6749d-7.
    real(dp):: band_degeneracy_threshold  = 3.6749d-5

    !> temperature lists:           10K        20K       70K      100K      200K      300K
    real(dp):: Eta_array_all(6) = (/0.00086d0, 0.0017d0, 0.006d0, 0.0086d0, 0.0172d0, 0.0259d0/)*eV2Hartree
    
    real(dp), parameter :: SOAHC_unit_factor = Echarge**3/hbar/Hartree2J
    real(dp), parameter :: INPHC_unit_factor = Echarge**3/hbar/Hartree2J * mu_B

    !> Differential step size on G and Lambda, in unit of [Length]^-1
    !> Too small value may lead to large error, our tests show that 1e-5 is better than 1e-6 and 1e-7 
    real(dp) :: dx = 1d-5
    real(dp) :: dy = 1d-5

contains
    subroutine velocity_latticegauge_simple(k, UU, velocities) !> dH_dk, without 1/hbar
        use para, only: irvec, crvec, HmnR, pi2zi, ndegen, Nrpts
        implicit none

        real(dp),    intent(in)  :: k(3)
        complex(dp), intent(in)  :: UU(Num_wann, Num_wann)
        complex(dp), intent(out) :: velocities(Num_wann, Num_wann, 3)

        real(dp):: kdotr
        integer :: iR
        complex(dp), allocatable :: Amat(:, :), UU_dag(:,:)
        allocate( Amat(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann))

        velocities= 0d0
        do iR= 1, Nrpts
            kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
            velocities(:,:,1)= velocities(:,:,1) + zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            velocities(:,:,2)= velocities(:,:,2) + zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            velocities(:,:,3)= velocities(:,:,3) + zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        enddo ! iR

        UU_dag= conjg(transpose(UU))
        !> unitility rotate velocity
        call mat_mul(Num_wann, velocities(:,:,1), UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, velocities(:,:,1))
        call mat_mul(Num_wann, velocities(:,:,2), UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, velocities(:,:,2))
        call mat_mul(Num_wann, velocities(:,:,3), UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, velocities(:,:,3))

        deallocate(Amat, UU_dag)
        return
    end subroutine velocity_latticegauge_simple

    
    subroutine Lambda_abc_df(k, energy, W, velocities, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)
        use para, only:  OmegaMin, OmegaMax
        implicit none

        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
        real(dp), intent(in)  :: W(Num_wann)
        complex(dp),intent(in):: velocities(Num_wann, Num_wann, 3)
        integer , intent(in)  :: NumberofEta
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: sigma_xyy_k(OmegaNum, NumberofEta)
        real(dp), intent(out) :: sigma_yxx_k(OmegaNum, NumberofEta)

        integer :: m, n, ie, ieta

        real(dp) :: mu, diffFermi
        real(dp) :: G_xy, G_yx, G_xx, G_yy
        complex(dp), allocatable :: vx(:, :), vy(:, :)

        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann) )
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)

        sigma_xyy_k        = 0d0
        sigma_yxx_k        = 0d0

        do m= 1, Num_wann
            !> At room temperature, the derivatives of the Fermi-Dirac distribution
            !> at (E-Ef)=0.5 is seven orders smaller than (E-Ef)=0.1.
            !> So we choose the energy truncation as [-0.5eV 0.5eV]*eV2Hartree,
            !> it will not affect the precsion and will accelerate the calculations
            if (W(m)<OmegaMin- 2.d-2 .or. W(m)>OmegaMax+ 2.d-2) cycle

            !> calculate G for each band
            G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0

            do n= 1, Num_wann
                if (ABS(W(m)-W(n)) < band_degeneracy_threshold) cycle
                G_xx= G_xx+ 2.d0*real(vx(m, n)*vx(n, m)/((W(m)-W(n))**3))
                G_xy= G_xy+ 2.d0*real(vx(m, n)*vy(n, m)/((W(m)-W(n))**3))
                G_yx= G_yx+ 2.d0*real(vy(m, n)*vx(n, m)/((W(m)-W(n))**3))
                G_yy= G_yy+ 2.d0*real(vy(m, n)*vy(n, m)/((W(m)-W(n))**3))
            enddo ! n

            !> consider the Fermi-distribution according to the brodening Earc_eta
            do ieta=1, NumberofEta
                do ie=1, OmegaNum
                    mu = energy(ie)

                    !> the if...else... statement here is to avoid infinite values,
                    !> generally, the 'Beta_fake*(W(m)-mu)<50' condition will be satisfied.
                    if ((W(m)-mu)/Eta_array(ieta)<50) then
                        diffFermi= -Exp((W(m)-mu)/Eta_array(ieta))/(Exp((W(m)-mu)/Eta_array(ieta))+1d0)**2 /Eta_array(ieta)
                    else
                        diffFermi=0.d0
                    endif

                    sigma_xyy_k(ie,ieta)= sigma_xyy_k(ie,ieta) &
                        + (G_yy*real(vx(m,m))-G_xy*real(vy(m,m)))*diffFermi
                    sigma_yxx_k(ie,ieta)= sigma_yxx_k(ie,ieta) &
                        + (G_xx*real(vy(m,m))-G_yx*real(vx(m,m)))*diffFermi
                enddo ! ie
            enddo ! ieta
        enddo ! m

        deallocate(vx, vy)
        return
    end subroutine Lambda_abc_df


    subroutine Lambda_abc_dG(k, energy, W, velocities, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)
        use para
        implicit none

        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
        real(dp), intent(in)  :: W(Num_wann)
        complex(dp),intent(inout):: velocities(Num_wann, Num_wann, 3)
        integer , intent(in)  :: NumberofEta
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: sigma_xyy_k(OmegaNum, NumberofEta)
        real(dp), intent(out) :: sigma_yxx_k(OmegaNum, NumberofEta)

        integer :: m, n, ie, ieta

        real(dp) :: mu, Fermi

        real(dp) :: k_dx(3)
        real(dp) :: k_dy(3)
        
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: Amat(:, :)
        complex(dp), allocatable :: UU(:, :)
        complex(dp), allocatable :: UU_dag(:, :)

        real(dp), allocatable :: W_dx(:)
        real(dp), allocatable :: W_dy(:)

        complex(dp), allocatable :: vx(:, :),    vy(:, :)
        complex(dp), allocatable :: vx_dx(:, :), vy_dx(:, :)  
        complex(dp), allocatable :: vx_dy(:, :), vy_dy(:, :) 

        real(dp) :: G_xx, G_xy, G_yx, G_yy, G_yy_dx, G_xy_dy, G_xx_dy, G_yx_dx

        allocate( Hamk_bulk (Num_wann, Num_wann))
        allocate( Amat (Num_wann, Num_wann))
        allocate( UU (Num_wann, Num_wann))
        allocate( UU_dag (Num_wann, Num_wann))

        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)

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

        sigma_xyy_k        = 0d0
        sigma_yxx_k        = 0d0

        do m= 1, Num_wann

            !> calculate G for each band
            G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0
            G_yy_dx=0d0; G_xy_dy=0d0; G_xx_dy=0d0; G_yx_dx=0d0

            do n= 1, Num_wann
                if (ABS(W(m)-W(n)) < band_degeneracy_threshold) cycle
                G_xx= G_xx+ 2.d0*real(vx(m, n)*vx(n, m)/((W(m)-W(n))**3))
                G_xy= G_xy+ 2.d0*real(vx(m, n)*vy(n, m)/((W(m)-W(n))**3))
                G_yx= G_yx+ 2.d0*real(vy(m, n)*vx(n, m)/((W(m)-W(n))**3))
                G_yy= G_yy+ 2.d0*real(vy(m, n)*vy(n, m)/((W(m)-W(n))**3))

                G_yy_dx= G_yy_dx + 2.d0*real( vy_dx(n, m)*vy_dx(m, n) )/(W_dx(m) - W_dx(n))**3
                G_yx_dx= G_yx_dx + 2.d0*real( vy_dx(n, m)*vx_dx(m, n) )/(W_dx(m) - W_dx(n))**3
                
                G_xy_dy= G_xy_dy + 2.d0*real( vx_dy(n, m)*vy_dy(m, n) )/(W_dy(m) - W_dy(n))**3
                G_xx_dy= G_xx_dy + 2.d0*real( vx_dy(n, m)*vx_dy(m, n) )/(W_dy(m) - W_dy(n))**3
            enddo ! n

            !> consider the Fermi-distribution according to the brodening Earc_eta
            do ieta=1, NumberofEta
                do ie=1, OmegaNum
                    mu = energy(ie)
                    Fermi= 1d0/(Exp((W(m)-mu)/Eta_array(ieta))+1d0)
                    
                    sigma_xyy_k(ie,ieta)= sigma_xyy_k(ie,ieta) &
                        - ((G_yy_dx - G_yy)/dx - (G_xy_dy - G_xy)/dy)*Fermi
                    sigma_yxx_k(ie,ieta)= sigma_yxx_k(ie,ieta) &
                        - ((G_xx_dy - G_xx)/dy - (G_yx_dx - G_yx)/dx)*Fermi
                enddo ! ie
            enddo ! ieta
        enddo ! m

        deallocate(vx, vy, Hamk_bulk, Amat, UU, UU_dag)
        deallocate(W_dx, W_dy, vx_dx, vx_dy, vy_dx, vy_dy)
        return
    end subroutine Lambda_abc_dG


    subroutine sigma_SOAHC_int_single_k(k, energy, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)

        implicit none
    
        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
        integer , intent(in)  :: NumberofEta
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: sigma_xyy_k(OmegaNum, NumberofEta)
        real(dp), intent(out) :: sigma_yxx_k(OmegaNum, NumberofEta)
    
        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: UU(:, :)
    
        !> velocities
        complex(dp), allocatable :: velocities(:,:,:)
    
        allocate( W(Num_wann))
        allocate( Hamk_bulk(Num_wann, Num_wann))
        allocate( UU(Num_wann, Num_wann))
        allocate( velocities(Num_wann, Num_wann, 3))
    
        Hamk_bulk= 0d0
        UU= 0d0
    
        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_latticegauge(k, Hamk_bulk)
    
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    
        call velocity_latticegauge_simple(k, UU, velocities)
    
        call Lambda_abc_df(k, energy, W, velocities, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)
    
        deallocate( W, Hamk_bulk, UU, velocities)
        return
    end subroutine sigma_SOAHC_int_single_k

end module

!> Calculate the intrinsic second order anomalous hall conductivity, the xyy and yxx elements
!
!> usage: sigma_SOAHC_int_calc = T
!
!> ref1 : 10.1103/PhysRevLett.127.277201
!> ref2 : 10.1103/PhysRevLett.127.277202
!
!> Original developed by Huiying Liu
!> 2022/07/15 Fan Yang, correct the units
!> 2023/10/30 Fan Yang, update to wannier tools 2.7.0
!> 2023/11/06 Fan Yang, adaptive k-meshes methods to accelerate speed
!

subroutine sigma_SOAHC_int
    use wmpi
    use para
    use nonlinear_transport
    implicit none

    integer :: NumberofEta           !> NumT
    real(dp):: adaptive_threshold_re !> including the SI units and the volume of the cell

    !> local k-points along each directions of the dense k-meshes
    integer :: Nk1_local
    integer :: Nk2_local
    integer :: Nk3_local
    integer :: knv3_local
    real(dp), allocatable :: klist_local(:,:)

    !> to decide the number of k-points which will be implemented with local dense k-meshes
    integer, allocatable :: displacement(:)
    integer ::              Nk_adaptive_mpi
    integer, allocatable :: Nk_adaptive(:) ! Nk_adaptive on every cores
    integer ::              Nk_adaptive_tol
    integer, allocatable :: klist_adaptive_mpi(:)
    integer, allocatable :: klist_adaptive    (:)

    !> estimate the error from the excluded k-points
    real(dp), allocatable :: sigma_included(:,:)
    real(dp), allocatable :: sigma_excluded(:,:)
    real(dp), allocatable :: sigma_included_mpi(:,:)
    real(dp), allocatable :: sigma_excluded_mpi(:,:)

    character*40 :: ahcfilename, etaname
    real(dp), allocatable :: Eta_array(:)

    integer :: ik, ik_local, ikx, iky, ikz, knv3
    integer :: ik_index
    integer :: ie, icore, ieta
    integer :: ierr

    real(dp) :: k(3)

    real(dp) :: time_start, time_end

    !> conductivity  dim = OmegaNum
    real(dp), allocatable :: energy(:)

    real(dp), allocatable :: sigma_xyy    (:,:)
    real(dp), allocatable :: sigma_yxx    (:,:)
    real(dp), allocatable :: sigma_xyy_k  (:,:)
    real(dp), allocatable :: sigma_yxx_k  (:,:)
    real(dp), allocatable :: sigma_xyy_mpi(:,:)
    real(dp), allocatable :: sigma_yxx_mpi(:,:)

    allocate( energy(OmegaNum))

    !> temperature, Eta = 1/k_B/T
    NumberofEta = NumT
    allocate( Eta_array(NumberofEta) )
    Eta_array(1) = Eta_Arc !> from wt.in
    if ((NumberofEta>1) .and. (NumberofEta<8)) then ! 1+6=7
        Eta_array(2:NumberofEta) = Eta_array_all(1:NumberofEta-1)
    else if (NumberofEta>7) then
        stop "The NumT should not more than 7"
    endif

    allocate( sigma_xyy        (OmegaNum, NumberofEta))
    allocate( sigma_yxx        (OmegaNum, NumberofEta))
    allocate( sigma_xyy_k      (OmegaNum, NumberofEta))
    allocate( sigma_yxx_k      (OmegaNum, NumberofEta))
    allocate( sigma_xyy_mpi    (OmegaNum, NumberofEta))
    allocate( sigma_yxx_mpi    (OmegaNum, NumberofEta))

    if (cpuid .eq. 0) then
        if (use_adaptive_method) then
            write(stdout, '("You have turned on the adaptive k-meshes feature.")')
            write(stdout, '("We advise to use a coarse k-mesh with a spacing of 0.005 A^(-1) in this case.")')
            write(stdout, '("(1/2): Testing the coarse k-mesh")')
        else
            write(stdout, '("You have turned off the adaptive k-meshes feature.")')
            write(stdout, '("We advise to use a dense k-mesh with a spacing of 0.001 A^(-1) in this case.")')
        endif
    endif

    !=============================================
    if (use_adaptive_method) then
        allocate( Nk_adaptive(num_cpu), displacement(num_cpu))
        allocate( klist_adaptive_mpi(Nk1*Nk2*Nk3), klist_adaptive(Nk1*Nk2*Nk3*num_cpu))

        allocate( sigma_included(OmegaNum,2))
        allocate( sigma_excluded(OmegaNum,2))
        allocate( sigma_included_mpi(OmegaNum,2))
        allocate( sigma_excluded_mpi(OmegaNum,2))

        !> generate the local k-meshes
        Nk1_local = Nk_local
        Nk2_local = Nk_local
        if (Nk3 == 1) then !> 2D system
            Nk3_local = 1
        else
            Nk3_local = Nk_local
        endif
        knv3_local = Nk1_local*Nk2_local*Nk3_local
        allocate( klist_local(knv3_local, 3) )

        do ik= 1, knv3_local
            ikx= (ik-1)/(Nk2_local*Nk3_local)+1
            iky= ((ik-1-(ikx-1)*Nk2_local*Nk3_local)/Nk3_local)+1
            ikz= (ik-(iky-1)*Nk3_local- (ikx-1)*Nk2_local*Nk3_local)
            klist_local(ik,:) = K3D_vec1_cube*(ikx-1)/dble(Nk1_local*Nk1)  &
                + K3D_vec2_cube*(iky-1)/dble(Nk2_local*Nk2)  &
                + K3D_vec3_cube*(ikz-1)/dble(Nk3_local*Nk3)
        enddo

        Nk_adaptive_mpi      = 0
        klist_adaptive_mpi   = 0

        adaptive_threshold_re = adaptive_threshold/(SOAHC_unit_factor/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume)

        sigma_included     = 0.d0
        sigma_excluded     = 0.d0
        sigma_included_mpi = 0.d0
        sigma_excluded_mpi = 0.d0
    else
        !> they must be allocated, so give the smallest szie
        allocate( Nk_adaptive(1), displacement(1))
        allocate( klist_adaptive_mpi(1), klist_adaptive(1))

        allocate( sigma_included(1,1))
        allocate( sigma_excluded(1,1))
        allocate( sigma_included_mpi(1,1))
        allocate( sigma_excluded_mpi(1,1))
        allocate( klist_local(1,1))
    endif
    !=====================================================

    !> energy
    do ie=1, OmegaNum
        if (OmegaNum>1) then
            energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
            energy= OmegaMin
        endif
    enddo ! ie

    knv3= Nk1*Nk2*Nk3

    sigma_xyy_mpi    = 0.d0
    sigma_yxx_mpi    = 0.d0

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
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(Nk3)

        call sigma_SOAHC_int_single_k(k, energy, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)

        if (.not. use_adaptive_method) then
            sigma_xyy_mpi = sigma_xyy_mpi + sigma_xyy_k
            sigma_yxx_mpi = sigma_yxx_mpi + sigma_yxx_k

            !> only use the block from the input Eta_Arc
        else if ((maxval(abs(sigma_xyy_k(:,1))) > adaptive_threshold_re) .or. (maxval(abs(sigma_yxx_k(:,1)))>adaptive_threshold_re)) then
            Nk_adaptive_mpi  = Nk_adaptive_mpi  + 1
            klist_adaptive_mpi(Nk_adaptive_mpi) = ik

            sigma_included_mpi(:,1) = sigma_included_mpi(:,1) + sigma_xyy_k(:,1)
            sigma_included_mpi(:,2) = sigma_included_mpi(:,2) + sigma_yxx_k(:,1)
        else
            sigma_excluded_mpi(:,1) = sigma_excluded_mpi(:,1) + sigma_xyy_k(:,1)
            sigma_excluded_mpi(:,2) = sigma_excluded_mpi(:,2) + sigma_yxx_k(:,1)
        endif
    enddo ! ik

    if (use_adaptive_method) then
        Nk_adaptive        = 0
        klist_adaptive     = 0
        displacement      = 0
#if defined (MPI)
        call mpi_barrier(mpi_cmw,ierr)
        call mpi_allreduce(sigma_included_mpi,sigma_included,size(sigma_included),&
            mpi_dp,mpi_sum,mpi_cmw,ierr)
        call mpi_allreduce(sigma_excluded_mpi,sigma_excluded,size(sigma_excluded),&
            mpi_dp,mpi_sum,mpi_cmw,ierr)
        call mpi_allgather(Nk_adaptive_mpi, 1, mpi_in, Nk_adaptive, 1, mpi_in, mpi_cmw,ierr)
        do icore=2, size(Nk_adaptive)
            displacement(icore)=sum(Nk_adaptive(1:icore-1))
        enddo
        call mpi_allgatherv(klist_adaptive_mpi, Nk_adaptive_mpi, mpi_in, klist_adaptive, Nk_adaptive, &
            displacement, mpi_in, mpi_cmw,ierr)
#else
        Nk_adaptive = Nk_adaptive_mpi
        klist_adaptive = klist_adaptive_mpi
        sigma_included = sigma_included_mpi
        sigma_excluded = sigma_excluded_mpi
#endif

        Nk_adaptive_tol = sum(Nk_adaptive)
        if (cpuid .eq. 0) then
            write(stdout, '(" ")')
            write(stdout, '("There are ", i15, "/", i18, "  k-points hit the threshold")') Nk_adaptive_tol, knv3
            write(stdout, '("The error from the excluded k-points is roughly ", f8.3, ", please check whether it is acceptable.")') &
                maxval(abs(sigma_excluded/sigma_included))
            write(stdout, '(" ")')
            write(stdout, '("(2/2): Scanning the local dense k-mesh")')
        endif

        call now(time_start)
        do ik_index = 1+ cpuid, Nk_adaptive_tol, num_cpu
            if (cpuid.eq.0.and. mod(ik_index/num_cpu, 100).eq.0) then
                call now(time_end)
                write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                    ik_index, Nk_adaptive_tol, '  time left', (Nk_adaptive_tol - ik_index)*(time_end-time_start)/num_cpu/100d0/60d0
                time_start= time_end
            endif

            ik = klist_adaptive(ik_index)
            ikx= (ik-1)/(Nk2*Nk3)+1
            iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
            ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
            k = K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
                + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
                + K3D_vec3_cube*(ikz-1)/dble(Nk3)

            do ik_local = 1, knv3_local
                call sigma_SOAHC_int_single_k(k+klist_local(ik_local,:), energy, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)

                sigma_xyy_mpi = sigma_xyy_mpi + sigma_xyy_k/dble(knv3_local)
                sigma_yxx_mpi = sigma_yxx_mpi + sigma_yxx_k/dble(knv3_local)
            enddo ! ik_local
        enddo ! ik
    endif ! use_adaptive_method


#if defined (MPI)
    call mpi_allreduce(sigma_xyy_mpi,sigma_xyy,size(sigma_xyy),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)

    call mpi_allreduce(sigma_yxx_mpi,sigma_yxx,size(sigma_yxx),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
    sigma_xyy = sigma_xyy_mpi
    sigma_yxx = sigma_yxx_mpi
#endif

    !> the sigma_xyy contains an additional [energy]^-1 dimension, so besides e^3/hbar, we need to convert hartree to joule
    sigma_xyy= sigma_xyy * SOAHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    sigma_yxx= sigma_yxx * SOAHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, NumberofEta
            write(etaname, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree
            write(ahcfilename, '(7a)')'sigma_ISOAHC_eta', trim(adjustl(etaname)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Intrinsic 2nd anomalous hall conductivity, in unit of A.V^-2 for 3D cases.'
            write(outfileindex, '("#",a)')' For 2D cases, you need to multiply the 3rd lattice vector in SI unit'
            write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyy', '\sigma_yxx'
            do ie=1, OmegaNum
                write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, sigma_xyy(ie,ieta), &
                    sigma_yxx(ie,ieta)
            enddo
            close(outfileindex)
        enddo
    endif

    deallocate( sigma_included, sigma_excluded, sigma_included_mpi, sigma_excluded_mpi)
    deallocate( Nk_adaptive, displacement, klist_adaptive_mpi, klist_adaptive, klist_local)
    deallocate( energy, Eta_array)
    deallocate( sigma_xyy, sigma_yxx, sigma_xyy_mpi, sigma_yxx_mpi )

    return
end subroutine sigma_SOAHC_int





