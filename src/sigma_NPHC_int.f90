subroutine sigma_NPHC_int
    !> Calculate the intrinsic nonlinear planar Hall conductivity, the xyyy and yxxx elements
    !
    !> usage: sigma_NPHC_int_calc = T
    !
    !> ref : 10.1103/PhysRevLett.130.126303
    !
    !> 2023/10/31 Fan Yang
    !
    !> The Lande g-factor is fixed to 2
    !> The input hr.dat must be from wannier90 v2.0 or newer versions

    use wmpi
    use para
    implicit none

    real(dp),parameter :: mu_B=9.274d-24 !> Bohr magneton in SI unit
    real(dp),parameter :: Lande_g=2d0

    integer :: ik, ikx, iky, ikz
    integer :: n, m, l, ie
    integer :: ierr, knv3
    real(dp) :: k(3)
    
    !> Fermi-Dirac distribution
    real(dp) :: mu, Beta_fake, diffFermi, diffFermi2, etmp

    real(dp) :: time_start, time_end
    
    ! eigen value of H
    real(dp), allocatable :: W(:)!> eig(:,:,:,:)
    real(dp), allocatable :: energy(:)!> Fermi energy, dim= OmegaNum
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

    complex(dp), allocatable :: Vmn_Ham(:, :, :)
    complex(dp), allocatable :: vx(:, :), vy(:, :), vz(:, :)!> velocities
    complex(dp), allocatable :: sx_oper(:, :), sy_oper(:, :), sz_oper(:, :)!> pauli matrix operator
    complex(dp), allocatable :: sx(:, :), sy(:, :), sz(:, :)! <psi| pauli matrix |psi>

    real(dp) :: Lambda_xyy, Lambda_yyy, Lambda_yxx, Lambda_xxx
    real(dp) :: G_xx, G_xy, G_yx, G_yy
    real(dp) :: Chi_xyyy, Chi_yxxx
    real(dp) :: dEnm, dEnm3, dEml, dEnl

    real(dp), allocatable :: Chi_xyyy_tensor(:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi(:)
    real(dp), allocatable :: Chi_yxxx_tensor(:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi(:)
        
    allocate( W (Num_wann))
    allocate( Vmn_Ham(Num_wann, Num_wann, 3))
    allocate( energy (OmegaNum))
    allocate( Hamk_bulk (Num_wann, Num_wann))
    allocate( Amat (Num_wann, Num_wann))
    allocate( UU (Num_wann, Num_wann))
    allocate( UU_dag (Num_wann, Num_wann))

    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann), vz(Num_wann, Num_wann))
    allocate( sx_oper(Num_wann, Num_wann), sy_oper(Num_wann, Num_wann), sz_oper(Num_wann, Num_wann))
    allocate( sx(Num_wann, Num_wann), sy(Num_wann, Num_wann), sz(Num_wann, Num_wann))

    allocate( Chi_xyyy_tensor    (OmegaNum))
    allocate( Chi_xyyy_tensor_mpi    (OmegaNum))
    allocate( Chi_yxxx_tensor    (OmegaNum))
    allocate( Chi_yxxx_tensor_mpi    (OmegaNum))

    Beta_fake=1d0/Eta_Arc
    diffFermi=0d0
    diffFermi2=0d0
    etmp=0d0

    Hamk_bulk= 0d0
    Amat= 0d0
    UU= 0d0
    UU_dag= 0d0   

    Chi_xyyy_tensor     = 0d0
    Chi_xyyy_tensor_mpi = 0d0
    Chi_yxxx_tensor     = 0d0
    Chi_yxxx_tensor_mpi = 0d0

    dEnm= 0d0
    dEnm3= 0d0
    dEml= 0d0
    dEnl= 0d0

    !> Fermi energy in Hatree energy, not eV
    do ie=1, OmegaNum
        if (OmegaNum>1) then
            energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie- 1d0)/dble(OmegaNum- 1)            
        else
            energy= OmegaMin
        endif
    enddo ! ie

    !> generate Pauli matrix
    sx_oper= 0d0; sy_oper= 0d0; sz_oper= 0d0
    do ie=1, Num_wann, 2
        sx_oper(ie  , ie+1)= (1d0, 0d0);    sx_oper(ie+1, ie  )= (1d0, 0d0)
        sy_oper(ie  , ie+1)= (0d0,-1d0);    sy_oper(ie+1, ie  )= (0d0, 1d0)
    enddo

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

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_latticegauge(k, Hamk_bulk)
  
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        UU_dag= conjg(transpose(UU))

        !> get velocity operator in Hamiltonian basis, without 1/hbar!!!
        call dHdk_latticegauge_Ham(k, W, UU, Vmn_Ham)
        vx = Vmn_Ham(:,:,1)
        vy = Vmn_Ham(:,:,2)
        vz = 0d0 ! Vmn_Ham(:,:,3)

        sx= 0d0; sy= 0d0; sz= 0d0
        call mat_mul(Num_wann, sx_oper, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, sx) 
        call mat_mul(Num_wann, sy_oper, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, sy) 
        ! call mat_mul(Num_wann, sz, UU, Amat) 
        ! call mat_mul(Num_wann, UU_dag, Amat, sz)

        do n= 1, Num_wann
            if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle 
            G_xx= 0d0
            G_xy= 0d0
            G_yx= 0d0
            G_yy= 0d0
            Lambda_xyy= 0d0
            Lambda_yyy= 0d0
            Lambda_yxx= 0d0
            Lambda_xxx= 0d0

            do m= 1, Num_wann
                if (m==n) cycle
                dEnm= W(n) - W(m)           
                if (ABS(dEnm)<1.d-6) cycle

                dEnm3= dEnm**3
                G_xx= G_xx+ 2.d0*real(vx(n, m)*vx(m, n)/dEnm3)
                G_xy= G_xy+ 2.d0*real(vx(n, m)*vy(m, n)/dEnm3)
                G_yx= G_yx+ 2.d0*real(vy(n, m)*vx(m, n)/dEnm3)
                G_yy= G_yy+ 2.d0*real(vy(n, m)*vy(m, n)/dEnm3) 

                Lambda_xyy= Lambda_xyy + 3.d0*real(vx(n, m)*vy(m, n)*(sy(n, n)-sy(m, m))/dEnm3/dEnm)
                Lambda_yyy= Lambda_yyy + 3.d0*real(vy(n, m)*vy(m, n)*(sy(n, n)-sy(m, m))/dEnm3/dEnm)
                Lambda_yxx= Lambda_yxx + 3.d0*real(vy(n, m)*vx(m, n)*(sx(n, n)-sx(m, m))/dEnm3/dEnm)
                Lambda_xxx= Lambda_xxx + 3.d0*real(vx(n, m)*vx(m, n)*(sx(n, n)-sx(m, m))/dEnm3/dEnm)
                
                do l= 1, Num_wann
                    dEml= W(m)-W(l)
                    dEnl= W(n)-W(l)
                    if (ABS(dEnl)>1.d-6) then
                        Lambda_xyy= Lambda_xyy - real((vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)*sy(n, l)) /dEnm3/dEnl)
                        Lambda_yyy= Lambda_yyy - real((vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)*sy(n, l)) /dEnm3/dEnl)
                        Lambda_yxx= Lambda_yxx - real((vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)*sx(n, l)) /dEnm3/dEnl)
                        Lambda_xxx= Lambda_xxx - real((vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)*sx(n, l)) /dEnm3/dEnl)
                    endif
                    if (ABS(dEml)>1.d-6) then
                        Lambda_xyy= Lambda_xyy - real((vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)*sy(m, l)) /dEnm3/dEml)
                        Lambda_yyy= Lambda_yyy - real((vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)*sy(m, l)) /dEnm3/dEml)
                        Lambda_yxx= Lambda_yxx - real((vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)*sx(m, l)) /dEnm3/dEml)
                        Lambda_xxx= Lambda_xxx - real((vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)*sx(m, l)) /dEnm3/dEml)
                    endif
                enddo ! l

            enddo ! m

            Lambda_xyy = Lambda_xyy * 2.d0
            Lambda_yyy = Lambda_yyy * 2.d0
            Lambda_yxx = Lambda_yxx * 2.d0
            Lambda_xxx = Lambda_xxx * 2.d0

            !> consider the Fermi-distribution according to the brodening Earc_eta
            do ie=1, OmegaNum
                mu = energy(ie)
                etmp=Exp(Beta_fake*(W(n)-mu))
                diffFermi = -Beta_fake   *etmp/(etmp+1d0)**2
                diffFermi2=  Beta_fake**2*2d0*Exp(2d0*Beta_fake*(W(n)-mu))/(etmp+1d0)**3 + Beta_fake*diffFermi
                
                Chi_xyyy = 0d0
                Chi_yxxx = 0d0

                Chi_xyyy = Chi_xyyy + real( vx(n,n)*Lambda_yyy - vy(n,n)*Lambda_xyy)          * diffFermi
                Chi_xyyy = Chi_xyyy + real((vx(n,n)*G_yy       - vy(n,n)*G_xy      )*sy(n,n)) * diffFermi2
                Chi_yxxx = Chi_yxxx + real( vy(n,n)*Lambda_xxx - vx(n,n)*Lambda_yxx)          * diffFermi
                Chi_yxxx = Chi_yxxx + real((vy(n,n)*G_xx       - vx(n,n)*G_yx      )*sx(n,n)) * diffFermi2

                
                Chi_xyyy_tensor_mpi(ie)= Chi_xyyy_tensor_mpi(ie) + Chi_xyyy
                Chi_yxxx_tensor_mpi(ie)= Chi_yxxx_tensor_mpi(ie) + Chi_yxxx
            enddo ! ie
        enddo ! n
          
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
    Chi_xyyy_tensor= Chi_xyyy_tensor *Echarge**3/hbar *(-Lande_g*mu_B) & !
        /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    Chi_yxxx_tensor= Chi_yxxx_tensor  & !*Echarge**2/hbar *(-Lande_g*mu_B)
        /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
       open(unit=outfileindex, file= 'sigma_NPHC_int.dat')
       write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of m*A*V^-2*T^-1'
       write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy', '\sigma_yxxx'
       do ie=1, OmegaNum
          write(outfileindex, '(200E17.8E3)') &
          energy(ie)/eV2Hartree, Chi_xyyy_tensor(ie), Chi_yxxx_tensor(ie)
       enddo
       close(outfileindex)
    endif

    deallocate(W, vx, vy, vz, Hamk_bulk, Amat, UU, UU_dag, energy)
    deallocate(sx_oper, sy_oper, sz_oper, sx, sy, sz)
    deallocate(Chi_xyyy_tensor, Chi_xyyy_tensor_mpi, Chi_yxxx_tensor, Chi_yxxx_tensor_mpi)

    return

end subroutine sigma_NPHC_int
