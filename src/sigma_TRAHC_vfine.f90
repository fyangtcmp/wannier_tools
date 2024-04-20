subroutine sigma_TRAHC
    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), parameter :: TRAHC_factor = Echarge**4/hbar**2 *Bohr_radius/Hartree2J

    real(dp), allocatable :: sigma_xxyy(:,:)
    real(dp), allocatable :: sigma_xxxx(:,:)
    real(dp), allocatable :: sigma_yyxx(:,:)
    real(dp), allocatable :: sigma_yyyy(:,:)
    real(dp), allocatable :: sigma_xxyytensor(:,:)
    real(dp), allocatable :: sigma_xxxxtensor(:,:)
    real(dp), allocatable :: sigma_yyxxtensor(:,:)
    real(dp), allocatable :: sigma_yyyytensor(:,:)
    real(dp), allocatable :: sigma_xxyytensor_mpi(:,:)
    real(dp), allocatable :: sigma_xxxxtensor_mpi(:,:)
    real(dp), allocatable :: sigma_yyxxtensor_mpi(:,:)
    real(dp), allocatable :: sigma_yyyytensor_mpi(:,:)

    real(dp), allocatable :: sigma_yyxxfine(:,:), sigma_yyyyfine(:,:)
    real(dp), allocatable :: sigma_xxyyfine(:,:), sigma_xxxxfine(:,:)

    allocate( energy(OmegaNum))

    allocate( sigma_xxyy             (OmegaNum, Eta_number))
    allocate( sigma_xxxx             (OmegaNum, Eta_number))
    allocate( sigma_yyxx             (OmegaNum, Eta_number))
    allocate( sigma_yyyy             (OmegaNum, Eta_number))
    allocate( sigma_yyxxtensor       (OmegaNum, Eta_number))
    allocate( sigma_yyyytensor       (OmegaNum, Eta_number))
    allocate( sigma_xxyytensor       (OmegaNum, Eta_number))
    allocate( sigma_xxxxtensor       (OmegaNum, Eta_number))
    allocate( sigma_xxyytensor_mpi   (OmegaNum, Eta_number))
    allocate( sigma_xxxxtensor_mpi   (OmegaNum, Eta_number))
    allocate( sigma_yyxxtensor_mpi   (OmegaNum, Eta_number))
    allocate( sigma_yyyytensor_mpi   (OmegaNum, Eta_number))

    allocate( sigma_xxyyfine         (OmegaNum, Eta_number))
    allocate( sigma_xxxxfine         (OmegaNum, Eta_number))
    allocate( sigma_yyxxfine         (OmegaNum, Eta_number))
    allocate( sigma_yyyyfine         (OmegaNum, Eta_number))

    allocate( Nk_adaptive (num_cpu))
    allocate( displacement(num_cpu))
    allocate( ik_adapt_list_mpi(Nk1*Nk2*Nk3))
    allocate( ik_adapt_list    (Nk1*Nk2*Nk3*num_cpu))
    
    sigma_yyxxtensor    = 0d0
    sigma_yyyytensor    = 0d0
    sigma_yyxxtensor_mpi    = 0d0
    sigma_yyyytensor_mpi    = 0d0

    sigma_xxyytensor    = 0d0
    sigma_xxxxtensor    = 0d0
    sigma_xxyytensor_mpi    = 0d0
    sigma_xxxxtensor_mpi    = 0d0

    call Fermi_energy_list(energy)

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    !> each k point
    knv3= int8(Nk1)*Nk2*Nk3

    Nk_adaptive         = 0
    Nk_adaptive_mpi     = 0
    displacement        = 0 
    ik_adapt_list_mpi  = 0
    ik_adapt_list      = 0

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif
        
        call ik_to_kpoint(ik,k)

        call sigma_TRAHC_k(k, sigma_xxxx, sigma_xxyy, sigma_yyxx, sigma_yyyy)

        if (maxval(abs(sigma_xxxx))>lim .or. maxval(abs(sigma_yyyy))>lim .or. maxval(abs(sigma_yyxx))>lim .or. maxval(abs(sigma_xxyy))>lim) then
            Nk_adaptive_mpi  = Nk_adaptive_mpi  + 1
            ik_adapt_list_mpi(Nk_adaptive_mpi) = ik
        else
            sigma_yyxxtensor_mpi = sigma_yyxxtensor_mpi + sigma_yyxx
            sigma_yyyytensor_mpi = sigma_yyyytensor_mpi + sigma_yyyy    
            sigma_xxyytensor_mpi = sigma_xxyytensor_mpi + sigma_xxyy
            sigma_xxxxtensor_mpi = sigma_xxxxtensor_mpi + sigma_xxxx         
        endif

    enddo ! ik

#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allgather(Nk_adaptive_mpi, 1, mpi_in, Nk_adaptive, 1, mpi_in, mpi_cmw,ierr)
    do icore=2, size(Nk_adaptive)
        displacement(icore)=sum(Nk_adaptive(1:icore-1))
    enddo
    call mpi_allgatherv(ik_adapt_list_mpi, Nk_adaptive_mpi, mpi_in, ik_adapt_list, Nk_adaptive, &
        displacement, mpi_in, mpi_cmw, ierr)
#endif

    Nk_adaptive_tol = sum(Nk_adaptive)
    if (cpuid .eq. 0) then
        write(stdout, '(" ")')
        write(stdout, '("There are ", i15, "/", i18, "  k-points hit the threshold")') Nk_adaptive_tol, knv3
        write(stdout, '(" ")')
        write(stdout, '("Start to scan the local fine k-grids")')
    endif

    if (Nk3<2) then
        knv3_fine = Nk_fine**2
    else 
        knv3_fine = Nk_fine**3
    endif
    allocate(k_fine_list(knv3_fine,3))

    do ikfine=1, knv3_fine
        if (Nk3<2) then
            ikfinex= (ikfine-1)/(Nk_fine)+1
            ikfiney= (ikfine-1-(ikfinex-1)*Nk_fine)+1
            ikfinez= 1
        else 
            ikfinex= (ikfine-1)/(Nk_fine*Nk_fine)+1
            ikfiney= ((ikfine-1-(ikfinex-1)*Nk_fine*Nk_fine)/Nk_fine)+1
            ikfinez= (ikfine-(ikfiney-1)*Nk_fine- (ikfinex-1)*Nk_fine*Nk_fine)
        endif
        k_fine_list(ikfine,:) = K3D_vec1_cube*(ikfinex-1)/dble(Nk1*Nk_fine)  &
            + K3D_vec2_cube*(ikfiney-1)/dble(Nk2*Nk_fine)  &
            + K3D_vec3_cube*(ikfinez-1)/dble(Nk3*Nk_fine)
    enddo
    
    call now(time_start)
    do ik2= 1+ cpuid, Nk_adaptive_tol, num_cpu
        if (cpuid.eq.0.and. mod(ik2/num_cpu, 200).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/Nk_adaptive_tol', &
            ik2, Nk_adaptive_tol, '  time left', (Nk_adaptive_tol-ik2)*(time_end-time_start)/num_cpu/200d0/60d0
            time_start= time_end
        endif
        ik = ik_adapt_list(ik2)
        call ik_to_kpoint(ik,k)

        sigma_yyxxfine=0d0
        sigma_yyyyfine=0d0
        sigma_xxyyfine=0d0
        sigma_xxxxfine=0d0
 
        do ikfine=1, knv3_fine
            call sigma_TRAHC_k( k + k_fine_list(ikfine,:), sigma_xxxx, sigma_xxyy, sigma_yyxx, sigma_yyyy )
            sigma_yyxxfine = sigma_yyxxfine + sigma_yyxx
            sigma_yyyyfine = sigma_yyyyfine + sigma_yyyy
            sigma_xxyyfine = sigma_xxyyfine + sigma_xxyy
            sigma_xxxxfine = sigma_xxxxfine + sigma_xxxx
        enddo
        sigma_yyxxtensor_mpi = sigma_yyxxtensor_mpi + sigma_yyxxfine/dble(knv3_fine)
        sigma_yyyytensor_mpi = sigma_yyyytensor_mpi + sigma_yyyyfine/dble(knv3_fine)
        sigma_xxyytensor_mpi = sigma_xxyytensor_mpi + sigma_xxyyfine/dble(knv3_fine)
        sigma_xxxxtensor_mpi = sigma_xxxxtensor_mpi + sigma_xxxxfine/dble(knv3_fine)
    enddo ! ik


#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allreduce(sigma_yyxxtensor_mpi,sigma_yyxxtensor,size(sigma_yyxxtensor),&
                    mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(sigma_yyyytensor_mpi,sigma_yyyytensor,size(sigma_yyyytensor),& 
                    mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(sigma_xxyytensor_mpi,sigma_xxyytensor,size(sigma_xxyytensor),&
                    mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(sigma_xxxxtensor_mpi,sigma_xxxxtensor,size(sigma_xxxxtensor),& 
                    mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif

    sigma_yyxxtensor= sigma_yyxxtensor/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor 
    sigma_yyyytensor= sigma_yyyytensor/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor 
    sigma_xxyytensor= sigma_xxyytensor/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor 
    sigma_xxxxtensor= sigma_xxxxtensor/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor  

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree
            write(ahcfilename, '(7a)')'sigma_TRAHC_eta', trim(adjustl(Eta_name)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Third-order anomalous hall conductivity without \tau, in unit of m.A.V^-3.s^-1 '
            write(outfileindex, '("#",a13, 20a16)') 'Energy (eV)', ' xxxx', ' xxyy', ' yyxx',  ' yyyy'
            do ie=1, OmegaNum
                write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                    sigma_xxxxtensor(ie,ieta), sigma_xxyytensor(ie,ieta), sigma_yyxxtensor(ie,ieta), sigma_yyyytensor(ie,ieta)
            enddo
            close(outfileindex)
        enddo
    endif

end subroutine sigma_TRAHC


subroutine sigma_TRAHC_k(kin, sigma_xxxxtensor, sigma_xxyytensor, sigma_yyxxtensor, sigma_yyyytensor)

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), intent(in)  :: kin(3)
    real(dp), intent(out) :: sigma_xxxxtensor(OmegaNum, Eta_number)
    real(dp), intent(out) :: sigma_xxyytensor(OmegaNum, Eta_number)
    real(dp), intent(out) :: sigma_yyxxtensor(OmegaNum, Eta_number)
    real(dp), intent(out) :: sigma_yyyytensor(OmegaNum, Eta_number)

    integer  :: iR, ikm
    real(dp) :: kdotr
    !real(dp) :: sqkvec1(3),sqkvec2(3),sqkvec3(3)
    !real(dp) :: kyp(3), kym(3)
    real(dp) :: kmat(3,3)
    real(dp) :: pyG_yy, pyG_yx, pyG_xx
    real(dp) :: pxG_xx, pxG_xy, pxG_yy

    !> Gij
    real(dp) :: G_xy(3),G_yx(3),G_xx(3),G_yy(3)

    !> eigen value of H
    real(dp), allocatable :: W(:),eig(:),eig2(:,:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

    !> velocities
    complex(dp), allocatable :: vx(:, :), vy(:, :)
    complex(dp), allocatable :: vxmat(:, :, :), vymat(:, :, :)
    complex(dp), allocatable :: vxx(:, :), vxy(:, :), vyy(:, :)
    ! real(dp) :: vxx_2, vxy_2, vyy_2

    real(dp), allocatable :: Fshort(:) !> short notation of the denominator of the Fermi distribution
    real(dp), allocatable :: diffFermi(:), diff2Fermi(:)

    allocate( Fshort(OmegaNum), diffFermi(OmegaNum), diff2Fermi(OmegaNum))

    allocate( W (Num_wann))
    allocate( eig(Num_wann))
    allocate( eig2(Num_wann,3))
    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))!, vz(Num_wann, Num_wann)
    allocate( vxx(Num_wann, Num_wann), vxy(Num_wann, Num_wann),vyy(Num_wann, Num_wann))
    allocate( vxmat(Num_wann, Num_wann,3), vymat(Num_wann, Num_wann,3))!, vzmat(Num_wann, Num_wann,Nk1,Nk2,Nk3)
    allocate( Hamk_bulk(Num_wann, Num_wann))
    allocate( Amat(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))
    allocate( UU_dag(Num_wann, Num_wann))

    Hamk_bulk=0d0
    Amat= 0d0
    UU_dag=0d0
    UU= 0d0
    eig=0d0
    eig2=0d0
    W=0d0

    kmat(1,:)=kin - (/Origin_cell%Rua(1)*dkx , Origin_cell%Rub(1)*dkx , Origin_cell%Ruc(1)*dkx/)/twopi
    kmat(2,:)=kin - (/Origin_cell%Rua(2)*dky , Origin_cell%Rub(2)*dky , Origin_cell%Ruc(2)*dky/)/twopi
    kmat(3,:)=kin

    do ikm= 1,3
        k=kmat(ikm,:)
        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_latticegauge(k, Hamk_bulk)
        
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        UU_dag= conjg(transpose(UU))
        eig2(:,ikm)=W
        vx= 0d0; vy= 0d0!; vz= 0d0
        vxx= 0d0; vxy =0d0; vyy=0d0
        do iR= 1, Nrpts
            kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
            vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            ! vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            ! if (ikm == 3) then
            !     vxx= vxx - crvec(1, iR)*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            !     vxy= vxy - crvec(1, iR)*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            !     vyy= vyy - crvec(2, iR)*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            ! endif
        enddo ! iR
    
        !> unitility rotate velocity
        UU_dag= conjg(transpose(UU))
        call mat_mul(Num_wann, vx, UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, vx)
        call mat_mul(Num_wann, vy, UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, vy)
        vxmat(:,:,ikm)=vx
        vymat(:,:,ikm)=vy
        ! if (ikm == 3) then
        !     call mat_mul(Num_wann, vxx, UU, Amat)
        !     call mat_mul(Num_wann, UU_dag, Amat, vxx)
        !     call mat_mul(Num_wann, vxy, UU, Amat)
        !     call mat_mul(Num_wann, UU_dag, Amat, vxy)
        !     call mat_mul(Num_wann, vyy, UU, Amat)
        !     call mat_mul(Num_wann, UU_dag, Amat, vyy)
        ! endif
    enddo !ikm

    eig=eig2(:,3)
    vx=vxmat(:,:,3)
    vy=vymat(:,:,3)
    sigma_xxxxtensor=0d0
    sigma_xxyytensor=0d0
    sigma_yyxxtensor=0d0
    sigma_yyyytensor=0d0

    do m= 1, Num_wann
        if (eig(m)<OmegaMin- 2.d-2 .or. eig(m)>OmegaMax+ 2.d-2) cycle !> prevent NaN error
        !> calculate G for each band
        !> G_xy == G_yx, we calculate both of them to check it
        G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0
        ! pyyG_xx=0d0; pyyG_yy=0d0; pyyG_yy=0d0
        ! vxx_2=0d0; vxy_2=0d0; vyy_2=0d0
        
        !> sum of all other bands n
        do n= 1, Num_wann
            do ikm = 1,3
                if (ABS(eig(m)-eig(n)) < band_degeneracy_threshold) cycle

                G_xx(ikm)= G_xx(ikm)+ 2.d0*real(vxmat(m, n, ikm)*vxmat(n, m, ikm)*((eig2(m,ikm)-eig2(n,ikm))/((eig2(m,ikm)-eig2(n,ikm))**2))**3)
                G_xy(ikm)= G_xy(ikm)+ 2.d0*real(vxmat(m, n, ikm)*vymat(n, m, ikm)*((eig2(m,ikm)-eig2(n,ikm))/((eig2(m,ikm)-eig2(n,ikm))**2))**3)
                G_yx(ikm)= G_yx(ikm)+ 2.d0*real(vymat(m, n, ikm)*vxmat(n, m, ikm)*((eig2(m,ikm)-eig2(n,ikm))/((eig2(m,ikm)-eig2(n,ikm))**2))**3)
                G_yy(ikm)= G_yy(ikm)+ 2.d0*real(vymat(m, n, ikm)*vymat(n, m, ikm)*((eig2(m,ikm)-eig2(n,ikm))/((eig2(m,ikm)-eig2(n,ikm))**2))**3)
                
                ! if (ikm == 3) then
                !     vxx_2= vxx_2+2d0*real(vx(m, n)*vx(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2+Eta_Arc**2))
                !     vxy_2= vxy_2+2d0*real(vx(m, n)*vy(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2+Eta_Arc**2))
                !     vyy_2= vyy_2+2d0*real(vy(m, n)*vy(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2+Eta_Arc**2))
                ! endif
            enddo
        enddo ! n
    
        pyG_yy=(G_yy(3)-G_yy(2))/dky
        pyG_xx=(G_xx(3)-G_xx(2))/dky
        pyG_yx=(G_xy(3)-G_xy(2))/dky  

        pxG_xx=(G_xx(3)-G_xx(1))/dkx
        pxG_xy=(G_yx(3)-G_yx(1))/dkx
        pxG_yy=(G_yy(3)-G_yy(1))/dkx
        
        do ieta=1, Eta_number
            Fshort = Exp((eig(m)-energy)/Eta_array(ieta))
            !> this format is very important! prevent NaN error
            diffFermi  = -1d0 / (Fshort+1d0) / (1d0/Fshort+1d0) / Eta_array(ieta)
            diff2Fermi = diffFermi * ( 1d0 - 2d0 / (1d0/Fshort+1d0)) / Eta_array(ieta)

            sigma_xxxxtensor(:,ieta) = sigma_xxxxtensor(:,ieta) + pxG_xx*real(vx(m,m))*diffFermi + real(vx(m,m)*vx(m,m))*G_xx(3)*diff2Fermi/2d0

            sigma_yyyytensor(:,ieta) = sigma_yyyytensor(:,ieta) + pyG_yy*real(vy(m,m))*diffFermi + real(vy(m,m)*vy(m,m))*G_yy(3)*diff2Fermi/2d0

            sigma_xxyytensor(:,ieta) = sigma_xxyytensor(:,ieta) +  real(vy(m,m))* (2d0*pxG_xy + pyG_xx) *diffFermi/3d0 &
                                                        + (real(vx(m,m)*vx(m,m))*G_yy(3) + 2d0*real(vx(m,m)*vy(m,m))*G_xy(3) )*diff2Fermi/6d0

            sigma_yyxxtensor(:,ieta) = sigma_yyxxtensor(:,ieta) +  real(vx(m,m))* (2d0*pyG_yx + pxG_yy) *diffFermi/3d0 &
                                                        + (real(vy(m,m)*vy(m,m))*G_xx(3)+2d0*real(vy(m,m)*vx(m,m))*G_yx(3)  )*diff2Fermi/6d0

        enddo ! ieta
    enddo !m

end subroutine sigma_TRAHC_k
