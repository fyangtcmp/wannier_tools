subroutine sigma_SOAHC_int
   !> Calculate the intrinsic second order anomalous hall conductivity, the xyy and yxx elements
   !
   !> usage: sigma_SOAHC_int_calc = T
   !
   !> ref1 : 10.1103/PhysRevLett.127.277201
   !> ref2 : 10.1103/PhysRevLett.127.277202
   !
   !> Original developed by Jianzhou Zhao
   !> 2022/07/15 Fan Yang, correct the units
   !> 2023/10/30 Fan Yang, update to wannier tools 2.7.0
   !
   
   use wmpi
   use para
   implicit none

   ! The energy threshold to judge the degenerate bands 
   ! The intrinsic SOAHC is usually calculated in PT-symmetric materials, with double degenerate bands.
   ! And for normal hr.dat, the numerical error between two physical degenerate bands is about 0.1~ meV.
   ! So we set it to 1meV = 3.6749d-5 Hartree, a relatively large value. It will lead to smoother curves of the conductivities.
   ! If you want larger peaks of the conductivities, or you have symmetrized your hr.dat to reduce the numerical error,
   ! you can set the threshold to a smaller value like 0.01meV = 3.6749d-7. 
   real(dp), parameter :: degeneracy_threshold = 3.6749d-5
   real(dp), parameter :: unit_factor = Echarge**3/hbar/Hartree2J
   real(dp), parameter :: adapted_threshold = 1d-7
   real(dp)            :: adapted_threshold_relative

   ! local k-points along each directions of the dense k-meshes
   integer :: Nk1_local = 5
   integer :: Nk2_local = 5
   integer :: Nk3_local
   integer :: knv3_local

   integer :: ik, ik_local, ikx, iky, ikz, knv3
   integer :: ik_index
   integer :: ie, icore
   integer :: ierr

   ! to decide the number of k-points which will be implemented with local dense k-meshes
   integer, allocatable :: displacement(:)
   integer ::              Nk_adapted_mpi
   integer, allocatable :: Nk_adapted(:) ! Nk_adapted on every cores
   integer ::              Nk_adapted_tol
   integer, allocatable :: klist_adapted_mpi(:)
   integer, allocatable :: klist_adapted    (:)
   
   real(dp) :: Beta_fake, k(3)
   
   real(dp) :: time_start, time_end

   real(dp), allocatable :: klist_local(:,:)
   !> conductivity  dim = OmegaNum
   real(dp), allocatable :: energy(:)

   real(dp), allocatable :: sigma_xyy(:)
   real(dp), allocatable :: sigma_yxx(:)

   real(dp), allocatable :: sigma_xyy_k(:)
   real(dp), allocatable :: sigma_yxx_k(:)

   real(dp), allocatable :: sigma_xyy_mpi(:)
   real(dp), allocatable :: sigma_yxx_mpi(:)

   allocate( Nk_adapted(num_cpu), displacement(num_cpu))
   allocate( klist_adapted_mpi(Nk1*Nk2), klist_adapted(Nk1*Nk2*num_cpu))
   
   allocate( energy(OmegaNum))

   allocate( sigma_xyy        (OmegaNum))
   allocate( sigma_yxx        (OmegaNum))
   allocate( sigma_xyy_k      (OmegaNum))
   allocate( sigma_yxx_k      (OmegaNum))
   allocate( sigma_xyy_mpi    (OmegaNum))
   allocate( sigma_yxx_mpi    (OmegaNum))

   knv3= Nk1*Nk2*Nk3

   Beta_fake= 1d0/Eta_Arc

   sigma_xyy_mpi    = 0d0
   sigma_yxx_mpi    = 0d0

   !=============================adapted k-mesh block======
   !> generate the local k-meshes
   if (Nk3 == 1) then !> 2D system
      Nk3_local = 1
   else
      Nk3_local = Nk1_local
   endif
   knv3_local = Nk1_local*Nk2_local*Nk3_local
   allocate( klist_local(knv3_local, 3) ) 

   do ik= 1, knv3_local
      ikx= (ik-1)/(Nk2_local*Nk3_local)+1
      iky= ((ik-1-(ikx-1)*Nk2_local*Nk3_local)/Nk3_local)+1
      ikz= (ik-(iky-1)*Nk3_local- (ikx-1)*Nk2_local*Nk3_local)
      klist_local(ik,:) = K3D_vec1_cube*(ikx-1)/dble(Nk1_local)  &
      + K3D_vec2_cube*(iky-1)/dble(Nk2_local)  &
      + K3D_vec3_cube*(ikz-1)/dble(Nk3_local)
   enddo

   !> shifted to centered k-meshes
   klist_local(:,1) = klist_local(:,1)/dble(Nk1)! - K3D_vec1_cube/dble(Nk1)/2.d0
   klist_local(:,2) = klist_local(:,2)/dble(Nk2)! - K3D_vec2_cube/dble(Nk2)/2.d0
   klist_local(:,3) = klist_local(:,3)/dble(Nk3)! - K3D_vec3_cube/dble(Nk3)/2.d0

   Nk_adapted_mpi      = 0
   klist_adapted_mpi   = 0

   adapted_threshold_relative = adapted_threshold/(unit_factor/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume)
   !=====================================================
   
   !> energy
   do ie=1, OmegaNum
      if (OmegaNum>1) then
         energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
      else
         energy= OmegaMin
      endif
   enddo ! ie

   write(stdout, '("(1/2): Testing the coarse k-mesh, the threshold for using a dense k-mesh is ", 14E10.2)') adapted_threshold
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

      call sigma_SOAHC_int_single_k(k, energy, degeneracy_threshold, sigma_xyy_k, sigma_yxx_k)
      
      ! sigma_xyy_mpi = sigma_xyy_mpi + sigma_xyy_k
      ! sigma_yxx_mpi = sigma_yxx_mpi + sigma_yxx_k

      if ((maxval(abs(sigma_xyy_k)) > adapted_threshold_relative) .or. (maxval(abs(sigma_yxx_k)) > adapted_threshold_relative)) then
         Nk_adapted_mpi  = Nk_adapted_mpi  + 1
         klist_adapted_mpi(Nk_adapted_mpi) = ik
      endif
   enddo ! ik

   Nk_adapted        = 0
   klist_adapted     = 0
   displacement      = 0

#if defined (MPI)
   call mpi_barrier(mpi_cmw,ierr)
   call mpi_allgather(Nk_adapted_mpi, 1, mpi_in, Nk_adapted, 1, mpi_in, mpi_cmw,ierr)
   do icore=2, size(Nk_adapted)
      displacement(icore)=sum(Nk_adapted(1:icore-1))
   enddo
   call mpi_allgatherv(klist_adapted_mpi, Nk_adapted_mpi, mpi_in, klist_adapted, Nk_adapted, &
      displacement, mpi_in, mpi_cmw,ierr)
#else
   Nk_adapted = Nk_adapted_mpi
   klist_adapted = klist_adapted_mpi
#endif

   Nk_adapted_tol = sum(Nk_adapted)
   write(stdout, '("There are ", i15, "/", i18, "  k-points hit the threshold")') Nk_adapted_tol, knv3
   write(stdout, '("(2/2): Scanning the local dense k-mesh")')

   call now(time_start) 
   do ik_index = 1+ cpuid, Nk_adapted_tol, num_cpu
      if (cpuid.eq.0.and. mod(ik_index/num_cpu, 100).eq.0) then
         call now(time_end) 
         write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
         ik_index, Nk_adapted_tol, '  time left', (Nk_adapted_tol - ik_index)*(time_end-time_start)/num_cpu/100d0/60d0
         time_start= time_end
      endif

      ik = klist_adapted(ik_index)
      ikx= (ik-1)/(Nk2*Nk3)+1
      iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
      ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
      k = K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
      + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
      + K3D_vec3_cube*(ikz-1)/dble(Nk3)

      do ik_local = 1, knv3_local
         call sigma_SOAHC_int_single_k(k+klist_local(ik_local,:), energy, degeneracy_threshold, sigma_xyy_k, sigma_yxx_k)
         
         sigma_xyy_mpi = sigma_xyy_mpi + sigma_xyy_k/dble(knv3_local)
         sigma_yxx_mpi = sigma_yxx_mpi + sigma_yxx_k/dble(knv3_local)
      enddo ! ik_local 
   enddo ! ik


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
   sigma_xyy= sigma_xyy * unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
          
   sigma_yxx= sigma_yxx * unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
          
   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='sigma_SOAHC_int_3D.dat')
      write(outfileindex, '("#",a)')' Intrinsic 2nd anomalous hall conductivity, in unit of A.V^-2 for 3D cases.'
      write(outfileindex, '("#",a)')' For 2D cases, you need to multiply the 3rd lattice vector in SI unit'
      write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyy', '\sigma_yxx'
      do ie=1, OmegaNum
         write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, sigma_xyy(ie), &
                                                      sigma_yxx(ie)
      enddo
      close(outfileindex)
   endif

   deallocate( Nk_adapted, displacement, klist_adapted_mpi, klist_adapted)
   deallocate( energy, klist_local)
   deallocate( sigma_xyy, sigma_yxx, sigma_xyy_mpi, sigma_yxx_mpi )

   return
end subroutine sigma_SOAHC_int

subroutine sigma_SOAHC_int_single_k(k, energy, degeneracy_threshold, sigma_xyy_k, sigma_yxx_k)
   
   use wmpi
   use para
   implicit none

   real(dp), intent(in)  :: k(3)
   real(dp), intent(in)  :: energy(OmegaNum)
   real(dp), intent(in)  :: degeneracy_threshold
   real(dp), intent(out) :: sigma_xyy_k(OmegaNum)
   real(dp), intent(out) :: sigma_yxx_k(OmegaNum)

   integer :: m, n, ie
   integer :: ierr
   
   real(dp) :: Beta_fake, mu, diffFermi
   real(dp) :: G_xy, G_yx, G_xx, G_yy
   
   ! eigen value of H
   real(dp),    allocatable :: W(:)
   complex(dp), allocatable :: Hamk_bulk(:, :)
   complex(dp), allocatable :: UU(:, :)
   
   !> velocities
   complex(dp), allocatable :: Vmn_Ham(:, :, :)
   complex(dp), allocatable :: vx(:, :), vy(:, :)     
       
   allocate( W(Num_wann))
   allocate( Hamk_bulk(Num_wann, Num_wann))
   allocate( UU(Num_wann, Num_wann))

   allocate( Vmn_Ham(Num_wann, Num_wann, 3))
   allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
   
   Hamk_bulk=0d0
   UU= 0d0

   sigma_xyy_k        = 0d0
   sigma_yxx_k        = 0d0

   Beta_fake= 1d0/Eta_Arc

   ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
   call ham_bulk_latticegauge(k, Hamk_bulk)

   !> diagonalization by call zheev in lapack
   UU=Hamk_bulk
   call eigensystem_c( 'V', 'U', Num_wann, UU, W)

   !> get velocity operator in Hamiltonian basis, without 1/hbar!!!
   call dHdk_latticegauge_Ham(k, W, UU, Vmn_Ham)
   vx = Vmn_Ham(:,:,1)
   vy = Vmn_Ham(:,:,2)

   do m= 1, Num_wann
      !> At room temperature, the derivatives of the Fermi-Dirac distribution
      !> at (E-Ef)=0.5 is seven orders smaller than (E-Ef)=0.1.
      !> So we choose the energy truncation as [-0.5eV 0.5eV],
      !> it will not affect the precsion and will accelerate the calculations
      if (W(m)<OmegaMin- 2.d-2 .or. W(m)>OmegaMax+ 2.d-2) cycle   

      !> calculate G for each band
      G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0

      do n= 1, Num_wann
         if (ABS(W(m)-W(n)) < degeneracy_threshold) cycle 
         G_xx= G_xx+ 2.d0*real(vx(m, n)*vx(n, m)/((W(m)-W(n))**3))
         G_xy= G_xy+ 2.d0*real(vx(m, n)*vy(n, m)/((W(m)-W(n))**3))
         G_yx= G_yx+ 2.d0*real(vy(m, n)*vx(n, m)/((W(m)-W(n))**3))
         G_yy= G_yy+ 2.d0*real(vy(m, n)*vy(n, m)/((W(m)-W(n))**3))                                 
      enddo ! n

      !> consider the Fermi-distribution according to the brodening Earc_eta
      do ie=1, OmegaNum
         mu = energy(ie)

         !> the if...else... statement here is to avoid infinite values,
         !> generally, the 'Beta_fake*(W(m)-mu)<50' condition will be satisfied.
         if (Beta_fake*(W(m)-mu)<50) then
            diffFermi=-Beta_fake*Exp(Beta_fake*(W(m)-mu))/(Exp(Beta_fake*(W(m)-mu))+1d0)**2
         else
            diffFermi=0.d0
         endif
         
         sigma_xyy_k(ie)= sigma_xyy_k(ie) &
            + (G_yy*real(vx(m,m))-G_xy*real(vy(m,m)))*diffFermi
         sigma_yxx_k(ie)= sigma_yxx_k(ie) &
            + (G_xx*real(vy(m,m))-G_yx*real(vx(m,m)))*diffFermi
      enddo ! ie
      
   enddo ! m

   deallocate( W, Vmn_Ham, vx, vy, Hamk_bulk, UU)
   return
end subroutine sigma_SOAHC_int_single_k
