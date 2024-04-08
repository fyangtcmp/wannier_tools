!> Do not parallel these codes over multiple nodes, cause it does not support the OPENMP.
!> the calculations on nonlinear transport are heavy, so the parallel version of wt.x is needed.

module nonlinear_transport
    use para, only: dp, eV2Hartree, Echarge, mu_B, Hartree2J, hbar, Num_wann, OmegaNum, zi, band_degeneracy_threshold, Eta_Arc
    implicit none

    !> magnetic moments in nonlinear planar Hall, see readinput.f90
    ! logical               :: include_m_spin = .false.
    ! logical               :: include_m_orb  = .true.

    !> Fermi-Dirac distribution
    real(dp) :: mu, Fermi, diffFermi

    !> eta = 8.617e-5*Temperature
    !> temperature lists:          20K       50K       70K      100K      200K      300K
    real(dp):: Eta_array_all(6) = [0.0017d0, 0.0043d0, 0.006d0, 0.0086d0, 0.0172d0, 0.0259d0]*eV2Hartree
    integer :: NumberofEta = 7
    
    character*40 :: ahcfilename, etaname

    real(dp) :: time_start, time_end
    
    real(dp), parameter :: SOAHC_unit_factor = Echarge**3/hbar/Hartree2J
    real(dp), parameter :: INPHC_unit_factor = Echarge**3/hbar/Hartree2J * mu_B

    !> loop 
    integer(8) :: n, m, l, ie, ieta, ikx, iky, ikz, ik, knv3
    integer    :: ierr

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

    
    subroutine Lambda_abc_df(k, energy, W, velocities, Eta_array, sigma_xyy_k, sigma_yxx_k)
        use para, only:  OmegaMin, OmegaMax
        implicit none

        real(dp), intent(in)  :: k(3)  !> not used, in order to keep the inputs the same as Lambda_abc_dG
        real(dp), intent(in)  :: energy(OmegaNum)
        real(dp), intent(in)  :: W(Num_wann)
        complex(dp),intent(in):: velocities(Num_wann, Num_wann, 3)
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: sigma_xyy_k(OmegaNum, NumberofEta)
        real(dp), intent(out) :: sigma_yxx_k(OmegaNum, NumberofEta)

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

                    !> very important! prevent NaN error
                    if ((W(m)-mu)/Eta_array(ieta)<50) then
                        diffFermi= Exp((W(m)-mu)/Eta_array(ieta))/(Exp((W(m)-mu)/Eta_array(ieta))+1d0)**2 /Eta_array(ieta)
                    else
                        diffFermi= 0.d0
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


    subroutine Lambda_abc_dG(k, energy, W, velocities, Eta_array, sigma_xyy_k, sigma_yxx_k)
        use para
        implicit none

        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
        real(dp), intent(in)  :: W(Num_wann)
        complex(dp),intent(inout):: velocities(Num_wann, Num_wann, 3)
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: sigma_xyy_k(OmegaNum, NumberofEta)
        real(dp), intent(out) :: sigma_yxx_k(OmegaNum, NumberofEta)

        real(dp) :: k_dx(3)
        real(dp) :: k_dy(3)
        
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: Amat(:, :)
        complex(dp), allocatable :: UU(:, :)

        real(dp), allocatable :: W_dx(:)
        real(dp), allocatable :: W_dy(:)

        complex(dp), allocatable :: vx(:, :),    vy(:, :)
        complex(dp), allocatable :: vx_dx(:, :), vy_dx(:, :)  
        complex(dp), allocatable :: vx_dy(:, :), vy_dy(:, :) 

        real(dp) :: G_xx, G_xy, G_yx, G_yy, G_yy_dx, G_xy_dy, G_xx_dy, G_yx_dx

        allocate( Hamk_bulk (Num_wann, Num_wann))
        allocate( Amat (Num_wann, Num_wann))
        allocate( UU (Num_wann, Num_wann))

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

        deallocate(vx, vy, Hamk_bulk, Amat, UU)
        deallocate(W_dx, W_dy, vx_dx, vx_dy, vy_dx, vy_dy)
        return
    end subroutine Lambda_abc_dG


    subroutine sigma_ISOAHC_single_k(k, energy, Eta_array, sigma_xyy_k, sigma_yxx_k)

        implicit none
    
        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
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
    
        call Lambda_abc_df(k, energy, W, velocities, Eta_array, sigma_xyy_k, sigma_yxx_k)
    
        deallocate( W, Hamk_bulk, UU, velocities)
        return
    end subroutine sigma_ISOAHC_single_k


    subroutine sigma_INPHC_single_k(k, energy, Eta_array, Chi_xyyy_k, Chi_yxxx_k)
        use magnetic_moments
        use para
        implicit none
    
        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: Chi_xyyy_k(OmegaNum, NumberofEta, 2, 2) !> the third index: 1=spin, 2=orbital
        real(dp), intent(out) :: Chi_yxxx_k(OmegaNum, NumberofEta, 2, 2)
    
        complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
        complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments
            
        real(dp) :: k_dx(3)
        real(dp) :: k_dy(3)
        
        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: Amat(:, :)
        complex(dp), allocatable :: UU(:, :)
        complex(dp), allocatable :: UU_dag(:, :)
    
        real(dp), allocatable :: W_dx(:)
        real(dp), allocatable :: W_dy(:)
    
        complex(dp), allocatable :: sx(:, :), sy(:, :) , sz(:, :)
        complex(dp), allocatable :: lx(:, :), ly(:, :) , lz(:, :)
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
            ! sz = -0.5d0 * Lande_g_S * M_S(:,:,3)
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W, velocities, M_L)
            lx = Lande_g_L * M_L(:,:,1)
            ly = Lande_g_L * M_L(:,:,2)
            ! lz = Lande_g_L * M_L(:,:,3) 
        endif    
    
        !> k + dk_x <===============================================================
        allocate( W_dx (Num_wann))   
        allocate( vx_dx(Num_wann, Num_wann), vy_dx(Num_wann, Num_wann))
    
        k_dx = k+(/Origin_cell%Rua(1)*dx , Origin_cell%Rub(1)*dx , Origin_cell%Ruc(1)*dx/)/twopi
    
        call ham_bulk_latticegauge(k_dx, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W_dx)
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
        call velocity_latticegauge_simple(k_dy, UU, velocities)
        vx_dy = velocities(:,:,1)
        vy_dy = velocities(:,:,2)
        !===========================================================================
    
        Chi_xyyy_k = 0d0
        Chi_yxxx_k = 0d0
    
        do n= 1, Num_wann
            if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle !> prevent NaN error
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

                    !> very important! prevent NaN error
                    if ((W(n)-mu)/Eta_array(ieta)<50) then
                        diffFermi= Exp((W(n)-mu)/Eta_array(ieta))/(Exp((W(n)-mu)/Eta_array(ieta))+1d0)**2 /Eta_array(ieta)
                    else
                        diffFermi= 0.d0
                    endif

                    if (include_m_spin) then
                        Chi_xyyy_k(ie,ieta, 1, 1) = Chi_xyyy_k(ie,ieta, 1, 1) + real( (vx(n,n)*Lambda_yyy_S - vy(n,n)*Lambda_xyy_S) ) * diffFermi
                        Chi_xyyy_k(ie,ieta, 1, 2) = Chi_xyyy_k(ie,ieta, 1, 2) + real( ((G_yy_dx - G_yy)/dx - (G_xy_dy - G_xy)/dy)*sy(n,n) ) * diffFermi
    
                        Chi_yxxx_k(ie,ieta, 1, 1) = Chi_yxxx_k(ie,ieta, 1, 1) + real( (vy(n,n)*Lambda_xxx_S - vx(n,n)*Lambda_yxx_S) ) * diffFermi
                        Chi_yxxx_k(ie,ieta, 1, 2) = Chi_yxxx_k(ie,ieta, 1, 2) + real( ((G_xx_dy - G_xx)/dy - (G_yx_dx - G_yx)/dx)*sx(n,n) ) * diffFermi
                    endif
                    if (include_m_orb) then
                        Chi_xyyy_k(ie,ieta, 2, 1) = Chi_xyyy_k(ie,ieta, 2, 1) + real( (vx(n,n)*Lambda_yyy_L - vy(n,n)*Lambda_xyy_L) ) * diffFermi
                        Chi_xyyy_k(ie,ieta, 2, 2) = Chi_xyyy_k(ie,ieta, 2, 2) + real( ((G_yy_dx - G_yy)/dx - (G_xy_dy - G_xy)/dy)*ly(n,n) ) * diffFermi
    
                        Chi_yxxx_k(ie,ieta, 2, 1) = Chi_yxxx_k(ie,ieta, 2, 1) + real( (vy(n,n)*Lambda_xxx_L - vx(n,n)*Lambda_yxx_L) ) * diffFermi
                        Chi_yxxx_k(ie,ieta, 2, 2) = Chi_yxxx_k(ie,ieta, 2, 2) + real( ((G_xx_dy - G_xx)/dy - (G_yx_dx - G_yx)/dx)*lx(n,n) ) * diffFermi
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


    subroutine drude_weight_single_k(k, energy, Eta_array, drude_k)
        use para
        implicit none
    
        real(dp), intent(in)  :: k(3)
        real(dp), intent(in)  :: energy(OmegaNum)
        real(dp), intent(in)  :: Eta_array(NumberofEta)
        real(dp), intent(out) :: drude_k(OmegaNum, NumberofEta, 2)
                
        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: UU(:, :)    
    
        complex(dp), allocatable :: velocities(:,:,:)
        complex(dp), allocatable :: vx(:, :), vy(:, :)

        allocate( W         (Num_wann))
        allocate( Hamk_bulk (Num_wann, Num_wann))
        allocate( UU        (Num_wann, Num_wann))
        allocate( velocities(Num_wann, Num_wann, 3))
        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
    
        call ham_bulk_latticegauge(k, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        call velocity_latticegauge_simple(k, UU, velocities)
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)
    
        drude_k = 0d0
    
        do n= 1, Num_wann
            if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle
            do ieta=1, NumberofEta
                do ie=1, OmegaNum
                    mu = energy(ie)

                    !> very important! prevent NaN error
                    if ((W(n)-mu)/Eta_array(ieta)<50) then
                        diffFermi= Exp((W(n)-mu)/Eta_array(ieta))/(Exp((W(n)-mu)/Eta_array(ieta))+1d0)**2 /Eta_array(ieta)
                    else
                        diffFermi= 0.d0
                    endif

                    drude_k(ie, ieta, 1) = drude_k(ie, ieta, 1) + real(vx(n,n))**2 *diffFermi
                    drude_k(ie, ieta, 2) = drude_k(ie, ieta, 2) + real(vy(n,n))**2 *diffFermi
                enddo ! ie
            enddo ! ieta
        enddo ! n
    
        deallocate(W, vx, vy, Hamk_bulk, UU, velocities)
        return
    end subroutine drude_weight_single_k

end module

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

    real(dp), allocatable :: Eta_array(:)

    real(dp) :: k(3)

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
    allocate( Eta_array(NumberofEta) )
    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:NumberofEta) = Eta_array_all(1:NumberofEta-1)

    allocate( sigma_xyy        (OmegaNum, NumberofEta))
    allocate( sigma_yxx        (OmegaNum, NumberofEta))
    allocate( sigma_xyy_k      (OmegaNum, NumberofEta))
    allocate( sigma_yxx_k      (OmegaNum, NumberofEta))
    allocate( sigma_xyy_mpi    (OmegaNum, NumberofEta))
    allocate( sigma_yxx_mpi    (OmegaNum, NumberofEta))

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

        call sigma_ISOAHC_single_k(k, energy, Eta_array, sigma_xyy_k, sigma_yxx_k)

        sigma_xyy_mpi = sigma_xyy_mpi + sigma_xyy_k
        sigma_yxx_mpi = sigma_yxx_mpi + sigma_yxx_k
    enddo ! ik

#if defined (MPI)
    call mpi_allreduce(sigma_xyy_mpi,sigma_xyy,size(sigma_xyy),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)

    call mpi_allreduce(sigma_yxx_mpi,sigma_yxx,size(sigma_yxx),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)
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

    real(dp) :: k(3)

    real(dp), allocatable :: Eta_array(:)

    real(dp), allocatable :: energy(:)  !> Fermi energy, dim= OmegaNum

    real(dp), allocatable :: Chi_xyyy_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_xyyy_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi(:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi(:,:,:,:)

    allocate( energy (OmegaNum))

    !> temperature, Eta = 1/k_B/T
    allocate( Eta_array(NumberofEta) )
    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:NumberofEta) = Eta_array_all(1:NumberofEta-1)

    allocate( Chi_xyyy_k         (OmegaNum, NumberofEta,2,2))
    allocate( Chi_yxxx_k         (OmegaNum, NumberofEta,2,2))
    allocate( Chi_xyyy_tensor    (OmegaNum, NumberofEta,2,2))
    allocate( Chi_yxxx_tensor    (OmegaNum, NumberofEta,2,2))
    allocate( Chi_xyyy_tensor_mpi(OmegaNum, NumberofEta,2,2))
    allocate( Chi_yxxx_tensor_mpi(OmegaNum, NumberofEta,2,2))

    Chi_xyyy_tensor     = 0d0
    Chi_yxxx_tensor     = 0d0
    Chi_xyyy_tensor_mpi = 0d0
    Chi_yxxx_tensor_mpi = 0d0

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

        call sigma_INPHC_single_k(k, energy, Eta_array, Chi_xyyy_k, Chi_yxxx_k)

        Chi_xyyy_tensor_mpi = Chi_xyyy_tensor_mpi + Chi_xyyy_k
        Chi_yxxx_tensor_mpi = Chi_yxxx_tensor_mpi + Chi_yxxx_k
    enddo ! ik

#if defined (MPI)
    call mpi_reduce(Chi_xyyy_tensor_mpi, Chi_xyyy_tensor, size(Chi_xyyy_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
    call mpi_reduce(Chi_yxxx_tensor_mpi, Chi_yxxx_tensor, size(Chi_yxxx_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
#endif

    Chi_xyyy_tensor = Chi_xyyy_tensor * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    Chi_yxxx_tensor = Chi_yxxx_tensor * INPHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, NumberofEta
            write(etaname, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            if (include_m_spin) then
                write(ahcfilename, '(7a)')'sigma_INPHC_S_eta', trim(adjustl(etaname)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1'
                write(outfileindex, '("#",a)')' Please refer to the Sec. III of the supplementary materials of 10.1103/PhysRevLett.130.126303, for the definition of term I and term II of the INPHE conductivities'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy_I', '\sigma_xyyy_II', '\sigma_xyyy_tot', '\sigma_yxxx_I', '\sigma_yxxx_II','\sigma_yxxx_tot'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,1,1), Chi_xyyy_tensor(ie,ieta,1,2), Chi_xyyy_tensor(ie,ieta,1,1) + Chi_xyyy_tensor(ie,ieta,1,2),&
                        Chi_yxxx_tensor(ie,ieta,1,1), Chi_yxxx_tensor(ie,ieta,1,2), Chi_yxxx_tensor(ie,ieta,1,1) + Chi_yxxx_tensor(ie,ieta,1,2)
                enddo
                close(outfileindex)
            endif

            if (include_m_orb ) then
                write(ahcfilename, '(7a)')'sigma_INPHC_L_eta', trim(adjustl(etaname)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1'
                write(outfileindex, '("#",a)')' For 2D cases, you need to multiply the 3rd lattice vector in SI unit'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy_I', '\sigma_xyyy_II', '\sigma_xyyy_tot', '\sigma_yxxx_I', '\sigma_yxxx_II','\sigma_yxxx_tot'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,2,1), Chi_xyyy_tensor(ie,ieta,2,2), Chi_xyyy_tensor(ie,ieta,2,1) + Chi_xyyy_tensor(ie,ieta,2,2),&
                        Chi_yxxx_tensor(ie,ieta,2,1), Chi_yxxx_tensor(ie,ieta,2,2), Chi_yxxx_tensor(ie,ieta,2,1) + Chi_yxxx_tensor(ie,ieta,2,2)
                enddo
                close(outfileindex)
            endif
        enddo
    endif

    deallocate(energy, Eta_array)
    deallocate(Chi_xyyy_k, Chi_yxxx_k, Chi_xyyy_tensor, Chi_xyyy_tensor_mpi, Chi_yxxx_tensor, Chi_yxxx_tensor_mpi)

    return

end subroutine sigma_INPHC

subroutine drude_weight

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp) :: k(3)
    
    real(dp), allocatable :: Eta_array(:)
    real(dp), allocatable :: energy(:) !> Fermi energy, dim= OmegaNum
    
    real(dp), allocatable :: drude      (:,:,:)
    real(dp), allocatable :: drude_mpi  (:,:,:)
    real(dp), allocatable :: drude_k    (:,:,:)

    allocate( energy(OmegaNum))

    allocate( drude    (OmegaNum, NumberofEta,2))
    allocate( drude_mpi(OmegaNum, NumberofEta,2))
    allocate( drude_k  (OmegaNum, NumberofEta,2))

    ! !> temperature, Eta = 1/k_B/T
    allocate( Eta_array (NumberofEta) )
    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:NumberofEta) = Eta_array_all(1:NumberofEta-1)

    ! !> Fermi energy in Hatree energy, not eV
    do ie=1, OmegaNum
        if (OmegaNum>1) then
            energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie- 1d0)/dble(OmegaNum- 1)
        else
            energy= OmegaMin
        endif
    enddo ! ie

    knv3= Nk1*Nk2*Nk3
    drude = 0d0
    drude_mpi = 0d0

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

        call drude_weight_single_k(k, energy, Eta_array, drude_k)

        drude_mpi = drude_mpi + drude_k
    enddo ! ik

#if defined (MPI)
    call mpi_reduce(drude_mpi, drude, size(drude_mpi), mpi_dp, mpi_sum, 0, mpi_cmw, ierr)
    call mpi_reduce(drude_mpi, drude, size(drude_mpi), mpi_dp, mpi_sum, 0, mpi_cmw, ierr)
#endif

    drude = drude * Echarge**2/hbar**2 *Hartree2J/Bohr_radius /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, NumberofEta
            write(etaname, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            write(ahcfilename, '(7a)')'drude_weight_eta', trim(adjustl(etaname)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Drude weight, in unit of S/m/(relaxation time)'
            write(outfileindex, '("#",a)')'  '

            write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', 'xx', 'yy' !, 'zz'
            do ie=1, OmegaNum
                write(outfileindex, '(20E16.5e3)')energy(ie)/eV2Hartree, drude(ie,ieta,1), drude(ie,ieta,2)
            enddo
            close(outfileindex)
        enddo
    endif

    deallocate(energy, Eta_array)
    deallocate(drude, drude_mpi, drude_k)
    return
end subroutine drude_weight

