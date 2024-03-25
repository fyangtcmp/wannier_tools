!> Do not parallel these codes over multiple nodes, cause it does not support the OPENMP.
!> the calculations on nonlinear transport are heavy, so the parallel version of wt.x is needed.

module nonlinear_transport
    use para, only: dp, eV2Hartree, Echarge, mu_B, Hartree2J, hbar, Num_wann, OmegaNum, zi, band_degeneracy_threshold
    implicit none

    !> magnetic moments in nonlinear planar Hall
    ! logical               :: include_m_spin = .false.
    ! logical               :: include_m_orb  = .true.

    !> adaptive k-meshes method
    logical :: use_adaptive_method = .false.
    real(dp):: adaptive_threshold  = 1e-10
    integer :: Nk_local = 6          !> local k-points along each directions of the dense k-meshes

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

        real(dp), intent(in)  :: k(3)  !> not used, in order to keep the inputs the same as Lambda_abc_dG
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


    subroutine sigma_ISOAHC_single_k(k, energy, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)

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
    end subroutine sigma_ISOAHC_single_k


    subroutine sigma_INPHC_single_k(k, energy, NumberofEta, Eta_array, Chi_xyyy_k_S, Chi_xyyy_k_L, Chi_yxxx_k_S, Chi_yxxx_k_L)
        use magnetic_moments
        use para
        implicit none
    
        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
        integer , intent(in)  :: NumberofEta
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: Chi_xyyy_k_S(OmegaNum, NumberofEta)
        real(dp), intent(out) :: Chi_yxxx_k_S(OmegaNum, NumberofEta)
        real(dp), intent(out) :: Chi_xyyy_k_L(OmegaNum, NumberofEta)
        real(dp), intent(out) :: Chi_yxxx_k_L(OmegaNum, NumberofEta)
    
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
        complex(dp), allocatable :: lx(:, :), ly(:, :), lz(:, :)
        complex(dp), allocatable :: vx(:, :), vy(:, :)
        complex(dp), allocatable :: velocities(:,:,:)
    
        complex(dp), allocatable :: vx_dx(:, :), vy_dx(:, :)  
        complex(dp), allocatable :: vx_dy(:, :), vy_dy(:, :) 
    
        real(dp) :: Lambda_xyy_S, Lambda_yyy_S, Lambda_yxx_S, Lambda_xxx_S
        real(dp) :: Lambda_xyy_L, Lambda_yyy_L, Lambda_yxx_L, Lambda_xxx_L
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
        allocate( lx(Num_wann, Num_wann), ly(Num_wann, Num_wann), lz(Num_wann, Num_wann))
    
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
            sx = -0.5d0 * Lande_g_S * M_S(:,:,1)
            sy = -0.5d0 * Lande_g_S * M_S(:,:,2)
            sz = -0.5d0 * Lande_g_S * M_S(:,:,3)
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W, velocities, M_L)
            lx = Lande_g_L * M_L(:,:,1)
            ly = Lande_g_L * M_L(:,:,2)
            lz = Lande_g_L * M_L(:,:,3) 
        endif    
    
        !> k + dk_x <===============================================================
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
    
        !> k + dk_y <===============================================================
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
    
        Chi_xyyy_k_S = 0d0
        Chi_yxxx_k_S = 0d0
        Chi_xyyy_k_L = 0d0
        Chi_yxxx_k_L = 0d0
    
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
            Lambda_xyy_S = 0d0
            Lambda_yyy_S = 0d0
            Lambda_yxx_S = 0d0
            Lambda_xxx_S = 0d0
            Lambda_xyy_L = 0d0
            Lambda_yyy_L = 0d0
            Lambda_yxx_L = 0d0
            Lambda_xxx_L = 0d0
    
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
    
                if (include_m_spin) then
                    Lambda_xyy_S = Lambda_xyy_S + 6.d0* real( vx(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                    Lambda_yyy_S = Lambda_yyy_S + 6.d0* real( vy(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                    Lambda_yxx_S = Lambda_yxx_S + 6.d0* real( vy(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                    Lambda_xxx_S = Lambda_xxx_S + 6.d0* real( vx(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                    
                    do l= 1, Num_wann
                        dEnl= W(n)-W(l)
                        dEml= W(m)-W(l)                    
                        if (ABS(dEnl) > band_degeneracy_threshold) then
                            Lambda_xyy_S = Lambda_xyy_S - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *sy(n, l)) /dEnm3/dEnl
                            Lambda_yyy_S = Lambda_yyy_S - 2.d0* real( (vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)) *sy(n, l)) /dEnm3/dEnl
                            Lambda_yxx_S = Lambda_yxx_S - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *sx(n, l)) /dEnm3/dEnl
                            Lambda_xxx_S = Lambda_xxx_S - 2.d0* real( (vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)) *sx(n, l)) /dEnm3/dEnl
                        endif
                        if (ABS(dEml) > band_degeneracy_threshold) then
                            Lambda_xyy_S = Lambda_xyy_S - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *sy(m, l)) /dEnm3/dEml
                            Lambda_yyy_S = Lambda_yyy_S - 2.d0* real( (vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)) *sy(m, l)) /dEnm3/dEml
                            Lambda_yxx_S = Lambda_yxx_S - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *sx(m, l)) /dEnm3/dEml
                            Lambda_xxx_S = Lambda_xxx_S - 2.d0* real( (vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)) *sx(m, l)) /dEnm3/dEml
                        endif
                    enddo ! l
                endif

                if (include_m_orb) then
                    Lambda_xyy_L = Lambda_xyy_L + 6.d0* real( vx(n, m)*vy(m, n)*(ly(n, n)-ly(m, m)) )/dEnm3/dEnm
                    Lambda_yyy_L = Lambda_yyy_L + 6.d0* real( vy(n, m)*vy(m, n)*(ly(n, n)-ly(m, m)) )/dEnm3/dEnm
                    Lambda_yxx_L = Lambda_yxx_L + 6.d0* real( vy(n, m)*vx(m, n)*(lx(n, n)-lx(m, m)) )/dEnm3/dEnm
                    Lambda_xxx_L = Lambda_xxx_L + 6.d0* real( vx(n, m)*vx(m, n)*(lx(n, n)-lx(m, m)) )/dEnm3/dEnm
                    
                    do l= 1, Num_wann
                        dEnl= W(n)-W(l)
                        dEml= W(m)-W(l)                    
                        if (ABS(dEnl) > band_degeneracy_threshold) then
                            Lambda_xyy_L = Lambda_xyy_L - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *ly(n, l)) /dEnm3/dEnl
                            Lambda_yyy_L = Lambda_yyy_L - 2.d0* real( (vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)) *ly(n, l)) /dEnm3/dEnl
                            Lambda_yxx_L = Lambda_yxx_L - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *lx(n, l)) /dEnm3/dEnl
                            Lambda_xxx_L = Lambda_xxx_L - 2.d0* real( (vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)) *lx(n, l)) /dEnm3/dEnl
                        endif
                        if (ABS(dEml) > band_degeneracy_threshold) then
                            Lambda_xyy_L = Lambda_xyy_L - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *ly(m, l)) /dEnm3/dEml
                            Lambda_yyy_L = Lambda_yyy_L - 2.d0* real( (vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)) *ly(m, l)) /dEnm3/dEml
                            Lambda_yxx_L = Lambda_yxx_L - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *lx(m, l)) /dEnm3/dEml
                            Lambda_xxx_L = Lambda_xxx_L - 2.d0* real( (vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)) *lx(m, l)) /dEnm3/dEml
                        endif
                    enddo ! l
                endif
    
            enddo ! m
    
            do ieta=1, NumberofEta
                do ie=1, OmegaNum
                    mu = energy(ie)
    
                    if ((W(n)-mu)/Eta_array(ieta)<50) then
                        diffFermi= -Exp((W(n)-mu)/Eta_array(ieta))/(Exp((W(n)-mu)/Eta_array(ieta))+1d0)**2 /Eta_array(ieta)
                    else
                        diffFermi=0.d0
                    endif
    
                    if (include_m_spin) then
                        Chi_xyyy_k_S(ie,ieta) = Chi_xyyy_k_S(ie,ieta) + real( (vx(n,n)*Lambda_yyy_S - vy(n,n)*Lambda_xyy_S) &
                            + ((G_yy_dx - G_yy)/dx - (G_xy_dy - G_xy)/dy)*sy(n,n) ) * diffFermi
                        Chi_yxxx_k_S(ie,ieta) = Chi_yxxx_k_S(ie,ieta) + real( (vy(n,n)*Lambda_xxx_S - vx(n,n)*Lambda_yxx_S) &
                            + ((G_xx_dy - G_xx)/dy - (G_yx_dx - G_yx)/dx)*sx(n,n) ) * diffFermi
                    endif
                    if (include_m_orb) then
                        Chi_xyyy_k_L(ie,ieta) = Chi_xyyy_k_L(ie,ieta) + real( (vx(n,n)*Lambda_yyy_L - vy(n,n)*Lambda_xyy_L) &
                            + ((G_yy_dx - G_yy)/dx - (G_xy_dy - G_xy)/dy)*ly(n,n) ) * diffFermi
                        Chi_yxxx_k_L(ie,ieta) = Chi_yxxx_k_L(ie,ieta) + real( (vy(n,n)*Lambda_xxx_L - vx(n,n)*Lambda_yxx_L) &
                            + ((G_xx_dy - G_xx)/dy - (G_yx_dx - G_yx)/dx)*lx(n,n) ) * diffFermi
                    endif    
    
                enddo ! ie
            enddo ! ieta
        enddo ! n
    
        deallocate(W, vx, vy, Hamk_bulk, Amat, UU, UU_dag, velocities)
        deallocate(W_dx, W_dy, vx_dx, vx_dy, vy_dx, vy_dy)
        deallocate(sx, sy, sz, lx, ly, lz)
        deallocate(M_S, M_L)
        return
    end subroutine sigma_INPHC_single_k

end module


#if defined (MPI)

subroutine sigma_ISOAHC

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

    integer :: ik, ik_local, ikx, iky, ikz 
    integer(8) :: knv3
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
        if (cpuid.eq.0.and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif

        ikx= (ik-1)/(Nk2*Nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(Nk3)

        call sigma_ISOAHC_single_k(k, energy, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)

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
                call sigma_ISOAHC_single_k(k+klist_local(ik_local,:), energy, NumberofEta, Eta_array, sigma_xyy_k, sigma_yxx_k)

                sigma_xyy_mpi = sigma_xyy_mpi + sigma_xyy_k/dble(knv3_local)
                sigma_yxx_mpi = sigma_yxx_mpi + sigma_yxx_k/dble(knv3_local)
            enddo ! ik_local
        enddo ! ik
    endif ! use_adaptive_method

    call mpi_allreduce(sigma_xyy_mpi,sigma_xyy,size(sigma_xyy),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)

    call mpi_allreduce(sigma_yxx_mpi,sigma_yxx,size(sigma_yxx),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)

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
end subroutine sigma_ISOAHC


subroutine sigma_INPHC
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

    integer :: NumberofEta              !> NumT

    integer :: ik, ikx, iky, ikz
    integer :: ie, ieta
    integer :: ierr
    integer(8) :: knv3
    real(dp) :: k(3)

    character*40 :: ahcfilename, etaname
    real(dp), allocatable :: Eta_array(:)

    real(dp), allocatable :: energy(:)  !> Fermi energy, dim= OmegaNum

    real(dp) :: time_start, time_end

    real(dp), allocatable :: Chi_xyyy_k_S         (:,:)
    real(dp), allocatable :: Chi_xyyy_k_L         (:,:)
    real(dp), allocatable :: Chi_yxxx_k_S         (:,:)
    real(dp), allocatable :: Chi_yxxx_k_L         (:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_S    (:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_L    (:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_S    (:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_L    (:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi_S(:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi_L(:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi_S(:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi_L(:,:)

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

    allocate( Chi_xyyy_k_S    (OmegaNum, NumberofEta))
    allocate( Chi_yxxx_k_S    (OmegaNum, NumberofEta))
    allocate( Chi_xyyy_k_L    (OmegaNum, NumberofEta))
    allocate( Chi_yxxx_k_L    (OmegaNum, NumberofEta))
    allocate( Chi_xyyy_tensor_S    (OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor_S    (OmegaNum, NumberofEta))
    allocate( Chi_xyyy_tensor_L    (OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor_L    (OmegaNum, NumberofEta))
    allocate( Chi_xyyy_tensor_mpi_S(OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor_mpi_S(OmegaNum, NumberofEta))
    allocate( Chi_xyyy_tensor_mpi_L(OmegaNum, NumberofEta))
    allocate( Chi_yxxx_tensor_mpi_L(OmegaNum, NumberofEta))

    Chi_xyyy_tensor_S     = 0d0
    Chi_yxxx_tensor_S     = 0d0
    Chi_xyyy_tensor_L     = 0d0
    Chi_yxxx_tensor_L     = 0d0
    Chi_xyyy_tensor_mpi_S = 0d0
    Chi_yxxx_tensor_mpi_S = 0d0
    Chi_xyyy_tensor_mpi_L = 0d0
    Chi_yxxx_tensor_mpi_L = 0d0

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
        if (cpuid.eq.0 .and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif

        ikx= (ik-1)/(Nk2*Nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube  &
            + K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(Nk3)

        call sigma_INPHC_single_k(k, energy, NumberofEta, Eta_array, Chi_xyyy_k_S, Chi_xyyy_k_L, Chi_yxxx_k_S, Chi_yxxx_k_L)

        Chi_xyyy_tensor_mpi_S = Chi_xyyy_tensor_mpi_S + Chi_xyyy_k_S
        Chi_yxxx_tensor_mpi_S = Chi_yxxx_tensor_mpi_S + Chi_yxxx_k_S
        Chi_xyyy_tensor_mpi_L = Chi_xyyy_tensor_mpi_L + Chi_xyyy_k_L
        Chi_yxxx_tensor_mpi_L = Chi_yxxx_tensor_mpi_L + Chi_yxxx_k_L
    enddo ! ik


    call mpi_allreduce(Chi_xyyy_tensor_mpi_S, Chi_xyyy_tensor_S, size(Chi_xyyy_tensor_S), mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Chi_yxxx_tensor_mpi_S, Chi_yxxx_tensor_S, size(Chi_yxxx_tensor_S), mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Chi_xyyy_tensor_mpi_L, Chi_xyyy_tensor_L, size(Chi_xyyy_tensor_L), mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Chi_yxxx_tensor_mpi_L, Chi_yxxx_tensor_L, size(Chi_yxxx_tensor_L), mpi_dp,mpi_sum,mpi_cmw,ierr)

    Chi_xyyy_tensor_S = Chi_xyyy_tensor_S * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    Chi_yxxx_tensor_S = Chi_yxxx_tensor_S * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    Chi_xyyy_tensor_L = Chi_xyyy_tensor_L * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    Chi_yxxx_tensor_L = Chi_yxxx_tensor_L * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, NumberofEta
            write(etaname, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            if (include_m_spin) then
                write(ahcfilename, '(7a)')'sigma_INPHC_S_eta', trim(adjustl(etaname)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1'
                write(outfileindex, '("#",a)')' For 2D cases, you need to multiply the 3rd lattice vector in SI unit'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy', '\sigma_yxxx'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, Chi_xyyy_tensor_S(ie,ieta), &
                        Chi_yxxx_tensor_S(ie,ieta)
                enddo
                close(outfileindex)
            endif

            if (include_m_orb ) then
                write(ahcfilename, '(7a)')'sigma_INPHC_L_eta', trim(adjustl(etaname)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1'
                write(outfileindex, '("#",a)')' For 2D cases, you need to multiply the 3rd lattice vector in SI unit'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy', '\sigma_yxxx'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, Chi_xyyy_tensor_L(ie,ieta), &
                        Chi_yxxx_tensor_L(ie,ieta)
                enddo
                close(outfileindex)
            endif
        enddo
    endif

    deallocate(energy, Eta_array)
    deallocate(Chi_xyyy_k_S, Chi_yxxx_k_S, Chi_xyyy_tensor_S, Chi_xyyy_tensor_mpi_S, Chi_yxxx_tensor_S, Chi_yxxx_tensor_mpi_S)
    deallocate(Chi_xyyy_k_L, Chi_yxxx_k_L, Chi_xyyy_tensor_L, Chi_xyyy_tensor_mpi_L, Chi_yxxx_tensor_L, Chi_yxxx_tensor_mpi_L)

    return

end subroutine sigma_INPHC

#endif