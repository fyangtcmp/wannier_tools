subroutine sigma_TRAHC
    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), parameter :: TRAHC_factor_tau1 = Echarge**4/hbar**2 *Bohr_radius/Hartree2J
    real(dp), parameter :: TRAHC_factor_tau3 = Echarge**4/hbar**4 *Bohr_radius/Hartree2J
    real(dp) :: lim=1d7

    real(dp), allocatable :: sigma           (:,:,:)   !> the second index = xxxx xxyy yyxx yyyy
    real(dp), allocatable :: sigma_tensor    (:,:,:)
    real(dp), allocatable :: sigma_tensor_mpi(:,:,:)
    real(dp), allocatable :: sigma_fine      (:,:,:)

    allocate( energy(OmegaNum))

    allocate( sigma              (OmegaNum, 8, Eta_number))
    allocate( sigma_tensor       (OmegaNum, 8, Eta_number))
    allocate( sigma_tensor_mpi   (OmegaNum, 8, Eta_number))
    allocate( sigma_fine         (OmegaNum, 8, Eta_number))
    
    sigma_tensor     = 0d0
    sigma_tensor_mpi = 0d0

    call Fermi_energy_list(energy)

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    !> each k point
    knv3= int8(Nk1)*Nk2*Nk3

    allocate( Nk_adapt_icore   (num_cpu))
    allocate( displacement     (num_cpu))
    allocate( ik_adapt_list_mpi(ceiling(real(knv3)/num_cpu)))         !> reduce memory
    allocate( ik_adapt_list    (ceiling(real(knv3)/num_cpu)*num_cpu)) !> reduce memory

    Nk_adapt_icore      = 0
    Nk_adapt_icore_mpi  = 0
    ik_adapt_list       = 0
    ik_adapt_list_mpi   = 0

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif
        
        call ik_to_kpoint(ik,k)

        call sigma_TRAHC_k(k, sigma)

        if (maxval(abs(sigma(:,1:4,:)))>lim) then
            Nk_adapt_icore_mpi  = Nk_adapt_icore_mpi  + 1
            ik_adapt_list_mpi(Nk_adapt_icore_mpi) = ik
        else
            sigma_tensor_mpi = sigma_tensor_mpi + sigma      
        endif

    enddo ! ik

    displacement = 0 
#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allgather(Nk_adapt_icore_mpi, 1, mpi_in, Nk_adapt_icore, 1, mpi_in, mpi_cmw, ierr) !> Use integer8 here will cause error, but the reason is unknown
    do icore=2, size(Nk_adapt_icore)
        displacement(icore)=sum(Nk_adapt_icore(1:icore-1))
    enddo
    call mpi_allgatherv(ik_adapt_list_mpi, Nk_adapt_icore_mpi, mpi_integer8, ik_adapt_list, Nk_adapt_icore, &
        displacement, mpi_integer8, mpi_cmw, ierr)
#endif

    Nk_adapt = sum(Nk_adapt_icore)
    if (cpuid .eq. 0) then
        write(stdout, '(" ")')
        write(stdout, '("There are ", i15, "/", i18, "  k-points hit the threshold")') Nk_adapt, knv3
        write(stdout, '(" ")')
        write(stdout, '("Start to scan the local fine k-grids")')
    endif

    if (Nk3<2) then
        knv3_fine = Nk_fine**2
    else 
        knv3_fine = Nk_fine**3
    endif
    allocate(k_fine_list(knv3_fine,3))

    k_fine_list = 0d0
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
    do ik2= 1+ cpuid, Nk_adapt, num_cpu
        if (cpuid.eq.0.and. mod(ik2/num_cpu, 200).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/Nk_adapt', &
            ik2, Nk_adapt, '  time left', (Nk_adapt-ik2)*(time_end-time_start)/num_cpu/200d0/60d0
            time_start= time_end
        endif
        ik = ik_adapt_list(ik2)
        call ik_to_kpoint(ik,k)

        sigma_fine=0d0
 
        do ikfine=1, knv3_fine
            call sigma_TRAHC_k( k + k_fine_list(ikfine,:), sigma )
            sigma_fine = sigma_fine + sigma
        enddo
        sigma_tensor_mpi = sigma_tensor_mpi + sigma_fine/dble(knv3_fine)
    enddo ! ik


#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allreduce(sigma_tensor_mpi,sigma_tensor,size(sigma_tensor),&
                    mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif

    sigma_tensor(:,1:4,:)= sigma_tensor(:,1:4,:)/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor_tau1 
    sigma_tensor(:,5:8,:)= sigma_tensor(:,5:8,:)/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor_tau3 

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree
            write(ahcfilename, '(7a)')'sigma_TRAHC_eta', trim(adjustl(Eta_name)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Third-order anomalous hall conductivity without tau, in the international System of Units '
            write(outfileindex, '("#",a13, 20a16)') 'Energy (eV) ', 'xxxx/tau', ' xxyy/tau', ' yyxx/tau',  ' yyyy/tau', &
                'xxxx/tau**3', ' xxyy/tau**3', ' yyxx/tau**3',  ' yyyy/tau**3'
            do ie=1, OmegaNum
                write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                    sigma_tensor(ie,1,ieta), sigma_tensor(ie,2,ieta), sigma_tensor(ie,3,ieta), sigma_tensor(ie,4,ieta),&
                    sigma_tensor(ie,5,ieta), sigma_tensor(ie,6,ieta), sigma_tensor(ie,7,ieta), sigma_tensor(ie,8,ieta)
            enddo
            close(outfileindex)
        enddo
    endif

end subroutine sigma_TRAHC


subroutine sigma_TRAHC_k(kin, sigma_tensor)

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), intent(in)  :: kin(3)
    real(dp), intent(out) :: sigma_tensor(OmegaNum, 8, Eta_number)

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
    real(dp), allocatable :: W(:),eig(:),eigmat(:,:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

    !> velocities
    complex(dp), allocatable :: vx(:, :), vy(:, :)
    complex(dp), allocatable :: vxmat(:, :, :), vymat(:, :, :)
    complex(dp), allocatable :: vxx(:, :), vxy(:, :), vyy(:, :)

    real(dp) :: vxx_2, vxy_2, vyy_2, exx, exy, eyy
    real(dp), allocatable :: Fshort(:) !> short notation of the denominator of the Fermi distribution
    real(dp), allocatable :: diffFermi(:), diff2Fermi(:)

    allocate( Fshort(OmegaNum), diffFermi(OmegaNum), diff2Fermi(OmegaNum))

    allocate( W (Num_wann))
    allocate( eig(Num_wann))
    allocate( eigmat(Num_wann,3))
    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
    allocate( vxmat(Num_wann, Num_wann,3), vymat(Num_wann, Num_wann,3))
    allocate( vxx(Num_wann, Num_wann), vxy(Num_wann, Num_wann), vyy(Num_wann, Num_wann))
    allocate( Hamk_bulk(Num_wann, Num_wann))
    allocate( Amat(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))
    allocate( UU_dag(Num_wann, Num_wann))

    Hamk_bulk=0d0
    Amat= 0d0
    UU_dag=0d0
    UU= 0d0
    eig=0d0
    eigmat=0d0
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
        eigmat(:,ikm)=W
        vx= 0d0; vy= 0d0!; vz= 0d0
        vxx= 0d0; vxy =0d0; vyy=0d0
        do iR= 1, Nrpts
            kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
            vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            if (ikm == 3) then
                vxx= vxx - crvec(1, iR)*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
                vxy= vxy - crvec(1, iR)*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
                vyy= vyy - crvec(2, iR)*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
            endif
        enddo ! iR
    
        !> unitility rotate velocity
        UU_dag= conjg(transpose(UU))
        call mat_mul(Num_wann, vx, UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, vx)
        call mat_mul(Num_wann, vy, UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, vy)
        vxmat(:,:,ikm)=vx
        vymat(:,:,ikm)=vy
        if (ikm == 3) then
            call mat_mul(Num_wann, vxx, UU, Amat)
            call mat_mul(Num_wann, UU_dag, Amat, vxx)
            call mat_mul(Num_wann, vxy, UU, Amat)
            call mat_mul(Num_wann, UU_dag, Amat, vxy)
            call mat_mul(Num_wann, vyy, UU, Amat)
            call mat_mul(Num_wann, UU_dag, Amat, vyy)
        endif
    enddo !ikm

    eig= eigmat(:,3)
    vx = vxmat(:,:,3)
    vy = vymat(:,:,3)
    sigma_tensor=0d0

    do m= 1, Num_wann
        if (eig(m)<OmegaMin- 2.d-2 .or. eig(m)>OmegaMax+ 2.d-2) cycle !> prevent NaN error
        !> calculate G for each band
        !> G_xy == G_yx, we calculate both of them to check it
        G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0
        vxx_2=0d0; vxy_2=0d0; vyy_2=0d0
        
        !> sum of all other bands n
        do n= 1, Num_wann

            do ikm = 1,3
                if (ABS(eig(m)-eig(n)) < band_degeneracy_threshold) cycle

                G_xx(ikm)= G_xx(ikm)+ 2.d0*real(vxmat(m, n, ikm)*vxmat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)
                G_xy(ikm)= G_xy(ikm)+ 2.d0*real(vxmat(m, n, ikm)*vymat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)
                G_yx(ikm)= G_yx(ikm)+ 2.d0*real(vymat(m, n, ikm)*vxmat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)
                G_yy(ikm)= G_yy(ikm)+ 2.d0*real(vymat(m, n, ikm)*vymat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)         
            enddo

            vxx_2= vxx_2+2d0*real(vx(m, n)*vx(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2))
            vxy_2= vxy_2+2d0*real(vx(m, n)*vy(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2))
            vyy_2= vyy_2+2d0*real(vy(m, n)*vy(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2))
        enddo ! n
    
        pyG_yy=(G_yy(3)-G_yy(2))/dky
        pyG_xx=(G_xx(3)-G_xx(2))/dky
        pyG_yx=(G_yx(3)-G_yx(2))/dky  

        pxG_xx=(G_xx(3)-G_xx(1))/dkx
        pxG_xy=(G_xy(3)-G_xy(1))/dkx
        pxG_yy=(G_yy(3)-G_yy(1))/dkx

        exx=real(vxx(m,m))+vxx_2
        exy=real(vxy(m,m))+vxy_2
        eyy=real(vyy(m,m))+vyy_2
        
        do ieta=1, Eta_number
            Fshort = Exp((eig(m)-energy)/Eta_array(ieta))
            !> this format is very important! prevent NaN error
            diffFermi  = -1d0 / (Fshort+1d0) / (1d0/Fshort+1d0) / Eta_array(ieta)
            diff2Fermi = diffFermi * ( 1d0 - 2d0 / (1d0/Fshort+1d0)) / Eta_array(ieta)

            ! xxxx/tau
            sigma_tensor(:,1,ieta) = sigma_tensor(:,1,ieta) + pxG_xx*real(vx(m,m))*diffFermi + real(vx(m,m)*vx(m,m))*G_xx(3)*diff2Fermi/2d0

            ! xxyy/tau
            sigma_tensor(:,2,ieta) = sigma_tensor(:,2,ieta) +  real(vy(m,m))* (2d0*pxG_xy + pyG_xx) *diffFermi/3d0 &
                + ( real(vx(m,m)*vx(m,m))*G_yy(3) + 2d0*real(vx(m,m)*vy(m,m))*G_xy(3) )*diff2Fermi/6d0

            ! yyxx/tau
            sigma_tensor(:,3,ieta) = sigma_tensor(:,3,ieta) +  real(vx(m,m))* (2d0*pyG_yx + pxG_yy) *diffFermi/3d0 &
                + ( real(vy(m,m)*vy(m,m))*G_xx(3) + 2d0*real(vy(m,m)*vx(m,m))*G_yx(3) )*diff2Fermi/6d0

            ! yyyy/tau
            sigma_tensor(:,4,ieta) = sigma_tensor(:,4,ieta) + pyG_yy*real(vy(m,m))*diffFermi + real(vy(m,m)*vy(m,m))*G_yy(3)*diff2Fermi/2d0

            ! xxxx/tau**3
            sigma_tensor(:,5,ieta) = sigma_tensor(:,5,ieta) + exx*exx*diffFermi + exx*real(vx(m,m)*vx(m,m))*diff2Fermi

            ! xxyy/tau**3
            sigma_tensor(:,6,ieta) = sigma_tensor(:,6,ieta) + (exx*eyy*diffFermi + exx*real(vy(m,m)*vy(m,m))*diff2Fermi &
                +2d0*exy*exy*diffFermi + 2d0*exy*real(vx(m,m)*vy(m,m))*diff2Fermi)/3d0

            ! yyxx/tau**3
            sigma_tensor(:,7,ieta) = sigma_tensor(:,7,ieta) + (eyy*exx*diffFermi + eyy*real(vx(m,m)*vx(m,m))*diff2Fermi &
                +2d0*exy*exy*diffFermi + 2d0*exy*real(vx(m,m)*vy(m,m))*diff2Fermi)/3d0
            
            ! yyyy/tau**3
            sigma_tensor(:,8,ieta) = sigma_tensor(:,8,ieta) + eyy*eyy*diffFermi + eyy*real(vy(m,m)*vy(m,m))*diff2Fermi
            
        enddo ! ieta
    enddo !m

end subroutine sigma_TRAHC_k
