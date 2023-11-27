module magnetic_moments
    use para, only: dp, zi, mu_B, eV2Hartree, Num_wann, Echarge, hbar, Bohr_radius
    use nonlinear_transport, only: band_degeneracy_threshold
    implicit none

    !> Lande g-factor
    real(dp), parameter :: Lande_g_S = 2d0
    real(dp), parameter :: Lande_g_L = 1d0

contains
    subroutine spin_magnetic_moments(M_S)
        !> extend the 2x2 pauli matrices to the wannier basis, without any units
        !> It is operators, so you need to do <psi| M_S |psi> yourself

        use para, only: Package
        implicit none

        integer :: j, nwann
        complex(dp), intent(out) :: M_S(Num_wann, Num_wann, 3)

        M_S= 0d0
        nwann = Num_wann/2

        !> generate Pauli matrix
        if (index( Package, 'VASP6'  )/=0 .or. &
            index( Package, 'QE')/=0) then

            do j=1, nwann
                M_S(2*j-1, 2*j,   1)=  1d0
                M_S(2*j,   2*j-1, 1)=  1d0
                M_S(2*j-1, 2*j,   2)=  -zi
                M_S(2*j,   2*j-1, 2)=   zi
                M_S(2*j-1, 2*j-1, 3)=  1d0
                M_S(2*j,   2*j,   3)= -1d0
            enddo
        else if (index( Package, 'VASP'  )/=0 .or. &
            index( Package, 'Wien2k')/=0 .or. &
            index( Package, 'Abinit')/=0 .or. &
            index( Package, 'openmx')/=0) then

            do j=1, nwann
                M_S(j,       nwann+j, 1)=  1d0
                M_S(j+nwann, j,       1)=  1d0
                M_S(j,       nwann+j, 2)=  -zi
                M_S(j+nwann, j,       2)=   zi
                M_S(j,       j,       3)=  1d0
                M_S(j+nwann, j+nwann, 3)= -1d0
            enddo
        endif
        return
    end subroutine spin_magnetic_moments

    subroutine orbital_magnetic_moments(W, velocities, M_L) !> without units, so we divide the M_L(Hartree/T) by mu_B(Hartree/T)
        !> SciPost Phys. 14, 118 (2023), Eq 24b
        !> It is NOT operators, which means M_tot = <psi| M_S |psi> + M_L

        implicit none

        integer :: l, n, p
        real(dp),    intent(in)  :: W(Num_wann)
        complex(dp), intent(in)  :: velocities(Num_wann, Num_wann,3)
        complex(dp), intent(out) :: M_L(Num_wann, Num_wann, 3)

        real(dp) :: dEpl, dEpn, inv_omega_plpn
        M_L = 0d0

        do l= 1, Num_wann
            do n= 1, Num_wann
                do p= 1, Num_wann
                    dEpl = W(p) - W(l)
                    dEpn = W(p) - W(n)
                    if ((ABS(dEpl) < band_degeneracy_threshold) .or. (ABS(dEpn) < band_degeneracy_threshold)) cycle
                    
                    inv_omega_plpn = (1/dEpl + 1/dEpn)
                    M_L(l,n,1) = M_L(l,n,1) + inv_omega_plpn * velocities(l,p,2) * velocities(p,n,3)
                    M_L(l,n,2) = M_L(l,n,2) + inv_omega_plpn * velocities(l,p,3) * velocities(p,n,1)
                    M_L(l,n,3) = M_L(l,n,3) + inv_omega_plpn * velocities(l,p,1) * velocities(p,n,2)
                enddo !p
            enddo !n
        enddo !l
        
        M_L = M_L * 4/zi * Echarge / hbar * Bohr_radius**2 /mu_B
        return
    end subroutine orbital_magnetic_moments
end module
