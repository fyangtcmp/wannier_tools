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

   integer :: ik, ikx, iky, ikz
   integer :: m, n, ie
   integer :: ierr, knv3
   
   real(dp) :: Beta_fake, mu, diffFermi
   real(dp) :: k(3)
   
   real(dp) :: time_start, time_end

   ! The energy threshold to judge the degenerate bands 
   ! The intrinsic SOAHC is usually calculated in PT-symmetric materials, with double degenerate bands.
   ! And for normal hr.dat, the numerical error between two physical degenerate bands is about 0.1~ meV.
   ! So we set it to 1meV = 3.6749d-5 Hartree, a relatively large value. It will lead to smoother curves of the conductivities.
   ! If you want larger peaks of the conductivities, or you have symmetrized your hr.dat to reduce the numerical error,
   ! you can set the threshold to a smaller value like 0.01meV = 3.6749d-7. 
   real(dp), parameter :: degeneracy_threshold = 3.6749d-5
   ! real(dp), parameter :: degeneracy_threshold = 3.6749d-7

   real(dp), parameter :: Hartree2J = 4.359748d-18 
   
   ! eigen value of H
   real(dp), allocatable    :: W(:)
   complex(dp), allocatable :: Hamk_bulk(:, :)
   complex(dp), allocatable :: Amat(:, :)
   complex(dp), allocatable :: UU(:, :)
   
   !> velocities
   complex(dp), allocatable :: Vmn_Ham(:, :, :)
   complex(dp), allocatable :: vx(:, :), vy(:, :)     
   !> Gij 
   real(dp):: G_xy,G_yx,G_xx,G_yy
   
   !> conductivity  dim = OmegaNum
   real(dp), allocatable :: energy(:)

   real(dp), allocatable :: sigma_xyytensor_ahc(:)
   real(dp), allocatable :: sigma_yxxtensor_ahc(:)

   real(dp), allocatable :: sigma_xyytensor_ahc_mpi(:)
   real(dp), allocatable :: sigma_yxxtensor_ahc_mpi(:)
         
   allocate( W(Num_wann))
   allocate( Hamk_bulk(Num_wann, Num_wann))
   allocate( Amat(Num_wann, Num_wann))
   allocate( UU(Num_wann, Num_wann))

   allocate( Vmn_Ham(Num_wann, Num_wann, 3))
   allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
   
   allocate( energy(OmegaNum))

   allocate( sigma_xyytensor_ahc        (OmegaNum))
   allocate( sigma_yxxtensor_ahc        (OmegaNum))
   allocate( sigma_xyytensor_ahc_mpi    (OmegaNum))
   allocate( sigma_yxxtensor_ahc_mpi    (OmegaNum))

   Hamk_bulk=0d0
   Amat= 0d0
   UU= 0d0

   sigma_xyytensor_ahc        = 0d0
   sigma_yxxtensor_ahc        = 0d0
   sigma_xyytensor_ahc_mpi    = 0d0
   sigma_yxxtensor_ahc_mpi    = 0d0

   Beta_fake= 1d0/Eta_Arc
   
   !> energy
   do ie=1, OmegaNum
      if (OmegaNum>1) then
         energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
      else
         energy= OmegaMin
      endif
   enddo ! ie

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
      k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1)  &
      + K3D_vec2_cube*(iky-1)/dble(Nk2)  &
      + K3D_vec3_cube*(ikz-1)/dble(Nk3)

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

            !> the if...else... statements here is for robustness,
            !> generally, the 'Beta_fake*(W(m)-mu)<50' condition will be satisfied.
            if (Beta_fake*(W(m)-mu)<50) then
               diffFermi=-Beta_fake*Exp(Beta_fake*(W(m)-mu))/(Exp(Beta_fake*(W(m)-mu))+1d0)**2
            else
               diffFermi=0.d0
            endif
            
            sigma_xyytensor_ahc_mpi(ie)= sigma_xyytensor_ahc_mpi(ie) &
               + (G_yy*real(vx(m,m))-G_xy*real(vy(m,m)))*diffFermi
            sigma_yxxtensor_ahc_mpi(ie)= sigma_yxxtensor_ahc_mpi(ie) &
               + (G_xx*real(vy(m,m))-G_yx*real(vx(m,m)))*diffFermi
         enddo ! ie
         
      enddo ! m
         
   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(sigma_xyytensor_ahc_mpi,sigma_xyytensor_ahc,size(sigma_xyytensor_ahc),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)

   call mpi_allreduce(sigma_yxxtensor_ahc_mpi,sigma_yxxtensor_ahc,size(sigma_yxxtensor_ahc),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   sigma_xyytensor_ahc= sigma_xyytensor_ahc_mpi
   sigma_yxxtensor_ahc= sigma_yxxtensor_ahc_mpi
#endif
   !> the sigma_xyytensor_ahc contains an additional [energy]^-1 dimension, so besides e^3/hbar, we need to convert hartree to joule
   sigma_xyytensor_ahc= sigma_xyytensor_ahc/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume &
      *Echarge**3/hbar/Hartree2J 
          
   sigma_yxxtensor_ahc= sigma_yxxtensor_ahc/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume &
      *Echarge**3/hbar/Hartree2J
          
   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='sigma_SOAHC_int_3D.dat')
      write(outfileindex, '("#",a)')' Intrinsic 2nd anomalous hall conductivity, in unit of A.V^-2 for 3D cases.'
      write(outfileindex, '("#",a)')' For 2D cases, you need to multiply the 3rd lattice vector in SI unit'
      write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyy', '\sigma_yxx'
      do ie=1, OmegaNum
         write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, sigma_xyytensor_ahc(ie), &
                                                      sigma_yxxtensor_ahc(ie)
      enddo
      close(outfileindex)
   endif

   deallocate( W, Vmn_Ham, vx, vy, Hamk_bulk, Amat, UU, energy)
   deallocate( sigma_xyytensor_ahc, sigma_yxxtensor_ahc, sigma_xyytensor_ahc_mpi, sigma_yxxtensor_ahc_mpi )

   return
end subroutine sigma_SOAHC_int
