!> Calculate the distributions of some band geometry properties which are used in the nonlinear transport
!> Do not parallel these codes over multiple nodes, cause it does not support the OPENMP.
!> the calculations on nonlinear transport are heavy, so the parallel version of wt.x is needed.


subroutine band_geo_props_kplane
    
    use wmpi
    use para
    use nonlinear_transport
    implicit none
   
    integer  :: i, j, nkmesh

    real(dp) :: kxy_plane(3)

    !> k points slice
    real(dp), allocatable :: kslice(:, :), kslice_xyz(:, :)
  
    real(dp), allocatable :: props(:,:), props_mpi(:,:)

    nkmesh= Nk1*Nk2
    allocate( props    (nkmesh, 6))
    allocate( props_mpi(nkmesh, 6))

    allocate( kslice(nkmesh, 3))
    allocate( kslice_xyz(nkmesh, 3))

    allocate( energy(1) )
    energy = 0d0

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    props= 0d0
    props_mpi= 0d0

    kslice=0d0
    kslice_xyz=0d0
   
    ik =0
    do i= 1, nk1
        do j= 1, nk2
            ik =ik +1
            kslice(ik, :)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                        + K3D_vec2*(j-1)/dble(nk2-1) - (K3D_vec1+ K3D_vec2)/2d0
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

#if defined (MPI)
    call mpi_allreduce(props_mpi, props, size(props), mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif

    !> write the Berry curvature to file
    outfileindex= outfileindex+ 1
    if (cpuid==0) then
        open(unit=outfileindex, file='band_geometry_kplane.dat')
        write(outfileindex, '("#",1a11, 2a12, 12a15)') "kx (1/A)", "ky (1/A)", "kz (1/A)", &
            'G_xx', 'G_xy', 'G_yx', 'G_yy', 'L_xyy', 'L_yxx' 

        ik= 0
        do i= 1, nk1
            do j= 1, nk2
                ik= ik+ 1
                ! call rotate_k3_to_kplane(kslice_xyz(ik, :), kxy_plane)
                write(outfileindex, '(3f12.6, 12E15.4e3)') kslice_xyz(ik, :)/Ang2Bohr, props(ik, :) ! kxy_plane*Angstrom2atomic, &
            enddo
            ! write(outfileindex, *) ' '
        enddo

        close(outfileindex)

    endif

end subroutine


subroutine ISOAHC_dist_single_k_Ef(k_in, props)
    
    use nonlinear_transport
    use para
    implicit none
   
    real(dp), intent(in)  :: k_in(3)
    real(dp), intent(out) :: props(6)
    real(dp) :: diffFermi

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

    call ham_bulk_latticegauge(k_in, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    UU_dag= conjg(transpose(UU))
    call velocity_latticegauge_simple(k_in, UU, velocities)
    vx = velocities(:,:,1)
    vy = velocities(:,:,2)

    props = 0d0

    m = Numoccupied
    !do m= 1, Num_wann
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
        
    !enddo ! m

    deallocate(W, vx, vy, Hamk_bulk, Amat, UU, UU_dag, velocities)
    return
end subroutine


subroutine INPHC_dist_single_k_Ef(k_in, props)

    use nonlinear_transport
    use magnetic_moments
    use para
    implicit none

    real(dp), intent(in)  :: k_in(3)
    real(dp), intent(out) :: props(6)

    real(dp), allocatable :: Chi_xyyy_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_xyyx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxy_k         (:,:,:,:)

    allocate( Chi_xyyy_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyx_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxy_k         (OmegaNum, Eta_number,2,2))
    
    call sigma_INPHC_single_k(  k_in, Chi_xyyy_k, Chi_yxxx_k)
    call sigma_INPHC_single_k_2(k_in, Chi_xyyx_k, Chi_yxxy_k)

    props(1) = sum(Chi_xyyy_k(1,1,:,:))
    props(2) = sum(Chi_yxxx_k(1,1,:,:))
    props(3) = sum(Chi_xyyx_k(1,1,:,:))
    props(4) = sum(Chi_yxxy_k(1,1,:,:))
    props(5) = 0d0
    props(6) = 0d0


end subroutine