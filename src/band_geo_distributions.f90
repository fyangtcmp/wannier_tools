!> Calculate the distributions of some band geometry properties which are used in the nonlinear transport
!> Do not parallel these codes over multiple nodes, cause it does not support the OPENMP.
!> the calculations on nonlinear transport are heavy, so the parallel version of wt.x is needed.


#if defined (MPI)
subroutine band_geo_props_kplane
    
    use wmpi
    use para
    implicit none
   
    integer  :: ik, i, j, ierr, nkmesh

    real(dp) :: k(3), kxy_plane(3)

    !> k points slice
    real(dp), allocatable :: kslice(:, :), kslice_xyz(:, :)
  
    real(dp) :: time_start, time_end

    real(dp), allocatable :: props(:,:), props_mpi(:,:)

    nkmesh= Nk1*Nk2
    allocate( props    (nkmesh, 6))
    allocate( props_mpi(nkmesh, 6))

    allocate( kslice(nkmesh, 3))
    allocate( kslice_xyz(nkmesh, 3))

    props= 0d0
    props_mpi= 0d0

    kslice=0d0
    kslice_xyz=0d0
   
    ik =0
    do i= 1, nk1
       do j= 1, nk2
          ik =ik +1
          kslice(ik, :)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                    + K3D_vec2*(j-1)/dble(nk2-1) ! - (K3D_vec1+ K3D_vec2)/2d0
          kslice_xyz(ik, :)= kslice(ik, 1)* Origin_cell%Kua+ kslice(ik, 2)* Origin_cell%Kub+ kslice(ik, 3)* Origin_cell%Kuc 
       enddo
    enddo

    call now(time_start)
    do ik= 1+ cpuid, nkmesh, num_cpu
        if (cpuid.eq.0 .and. mod(ik/num_cpu, 200).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, nkmesh, '  time left', (nkmesh-ik)*(time_end-time_start)/num_cpu/200d0/60d0
            time_start= time_end
        endif 

        k= kslice(ik, :)
        call ISOAHC_dist_single_k_Ef(k,  props_mpi(ik, 1:6))
        ! call INPHC_dist_single_k_Ef(k,  props_mpi(ik, 1:6))

    enddo ! ik

    call mpi_allreduce(props_mpi, props, size(props), mpi_dp,mpi_sum,mpi_cmw,ierr)

    !> write the Berry curvature to file
    outfileindex= outfileindex+ 1
    if (cpuid==0) then
        open(unit=outfileindex, file='band_geometry_kplane.dat')
        write(outfileindex, '(9a12)')'# k1', 'k2', 'k3', & ! "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", &
            'G_xx', 'G_xy', 'G_yx', 'G_yy', 'L_xyy', 'L_yxx' 

        ik= 0
        do i= 1, nk1
            do j= 1, nk2
                ik= ik+ 1
                ! call rotate_k3_to_kplane(kslice_xyz(ik, :), kxy_plane)
                write(outfileindex, '(3f9.3,6E12.3e3)') kslice(ik, :), props(ik, :) ! kxy_plane*Angstrom2atomic, &
            enddo
            ! write(outfileindex, *) ' '
        enddo

        close(outfileindex)

    endif

    deallocate( kslice, kslice_xyz )
    deallocate( props, props_mpi )

    return

end subroutine
#endif


subroutine ISOAHC_dist_single_k_Ef(k, props)
    
    use nonlinear_transport
    use para
    implicit none
   
    real(dp), intent(in)  :: k(3)
    real(dp), intent(out) :: props(6)

    ! eigen value of H
    real(dp),    allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

    complex(dp), allocatable :: vx(:, :), vy(:, :)
    complex(dp), allocatable :: velocities(:,:,:)

    real(dp) :: G_xy, G_yx, G_xx, G_yy

    allocate( W (Num_wann))
    allocate( Hamk_bulk (Num_wann, Num_wann))
    allocate( Amat (Num_wann, Num_wann))
    allocate( UU (Num_wann, Num_wann))
    allocate( UU_dag (Num_wann, Num_wann))

    allocate( velocities(Num_wann, Num_wann, 3))
    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))

    call ham_bulk_latticegauge(k, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    UU_dag= conjg(transpose(UU))
    call velocity_latticegauge_simple(k, UU, velocities)
    vx = velocities(:,:,1)
    vy = velocities(:,:,2)

    props = 0d0

    do m= 1, Num_wann
        !> calculate G for each band
        G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0

        do n= 1, Num_wann
            if (ABS(W(m)-W(n)) < band_degeneracy_threshold) cycle
            G_xx= G_xx+ 2.d0*real(vx(m, n)*vx(n, m)/((W(m)-W(n))**3))
            G_xy= G_xy+ 2.d0*real(vx(m, n)*vy(n, m)/((W(m)-W(n))**3))
            G_yx= G_yx+ 2.d0*real(vy(m, n)*vx(n, m)/((W(m)-W(n))**3))
            G_yy= G_yy+ 2.d0*real(vy(m, n)*vy(n, m)/((W(m)-W(n))**3))
        enddo ! n
        
        if (W(m)/Eta_Arc<50) then
            diffFermi= -Exp(W(m)/Eta_Arc)/(Exp(W(m)/Eta_Arc)+1d0)**2 /Eta_Arc
        else
            diffFermi=0d0
        endif

        props(1) = props(1) + G_xx * diffFermi
        props(2) = props(2) + G_xy * diffFermi
        props(3) = props(3) + G_yx * diffFermi
        props(4) = props(4) + G_yy * diffFermi
        props(5) = props(5) + (G_yy*real(vx(m,m))-G_xy*real(vy(m,m))) * diffFermi
        props(6) = props(6) + (G_xx*real(vy(m,m))-G_yx*real(vx(m,m))) * diffFermi
        
    enddo ! m

    deallocate(W, vx, vy, Hamk_bulk, Amat, UU, UU_dag, velocities)
    return
end subroutine


subroutine INPHC_dist_single_k_Ef(k, props)

    use nonlinear_transport
    use magnetic_moments
    use para
    implicit none

    real(dp), intent(in)  :: k(3)
    real(dp), intent(out) :: props(6)

    complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
    complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments
    
    ! eigen value of H
    real(dp),    allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

    complex(dp), allocatable :: sx(:, :), sy(:, :), sz(:, :)
    complex(dp), allocatable :: vx(:, :), vy(:, :)
    complex(dp), allocatable :: velocities(:,:,:)

    real(dp) :: G_xx, G_xy, G_yx, G_yy
    real(dp) :: Lambda_xyy, Lambda_yyy, Lambda_yxx, Lambda_xxx
    real(dp) :: dEnm, dEnm3, dEml, dEnl
    
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

    do n= 1, Num_wann
        G_xx= 0d0
        G_xy= 0d0
        G_yx= 0d0
        G_yy= 0d0
       
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

        if (W(n)/Eta_Arc<50) then
            diffFermi= -Exp(W(n)/Eta_Arc)/(Exp(W(n)/Eta_Arc)+1d0)**2 /Eta_Arc
        else
            diffFermi=0d0
        endif

        props(1) = props(1) + G_xx * diffFermi
        props(2) = props(2) + G_xy * diffFermi
        props(3) = props(3) + G_yx * diffFermi
        props(4) = props(4) + G_yy * diffFermi
        props(5) = props(5) + Lambda_xyy * diffFermi
        props(6) = props(6) + Lambda_yxx * diffFermi
    enddo ! n

    deallocate(W, vx, vy, Hamk_bulk, Amat, UU, UU_dag, velocities)
    deallocate(sx, sy, sz)
    deallocate(M_S, M_L)
    return

end subroutine