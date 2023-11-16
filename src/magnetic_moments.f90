module magnetic_moments
    use para, only: dp, eV2Hartree, Echarge, hbar, Num_wann
    implicit none

    !> Lande g-factor
    real(dp), parameter :: Lande_g_S = 2d0
    real(dp), parameter :: Lande_g_L = 1d0

    real(dp), parameter :: mu_B = 5.788d-5*eV2Hartree !> Bohr magneton, Hartree/Tesla

contains
    subroutine spin_magnetic_moments(M_S)
        !> extend the pauli matrices to the wannier basis, without any units

        use para, only: Package, Num_wann, zi
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
    end subroutine spin_magnetic_moments

    subroutine orbital_magnetic_moments(M_L)
        use para, only: Package, Num_wann, zi
        implicit none

        integer :: j, nwann
        complex(dp), intent(out) :: M_L(Num_wann, Num_wann, 3)

        M_L = 0d0
        
    end subroutine orbital_magnetic_moments
end module
