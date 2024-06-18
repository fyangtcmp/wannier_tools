# WannierTools NL

A modified edition of **WannierTools**. It supports the calculations of some nonlinear transport properties. 

## Citation
Please cite the original edition of **WannierTools**

```
@article{WU2018,
title = "WannierTools : An open-source software package for novel topological materials",
journal = "Computer Physics Communications",
volume = "224",
pages = "405 - 416",
year = "2018",
doi = "https://doi.org/10.1016/j.cpc.2017.09.033",
url = "http://www.sciencedirect.com/science/article/pii/S0010465517303442",
issn = "0010-4655",
preprint = "arXiv:1703.07789",
author = "QuanSheng Wu and ShengNan Zhang and Hai-Feng Song and Matthias Troyer and Alexey A. Soluyanov",
keywords = "Novel topological materials, Topological number, Surface state, Tight-binding model"
}
```

## Usage

### Drude weight ($\sigma_{xx}/\tau$)
```
&CONTROL
drude_weight_calc = .TRUE.
/

&SYSTEM
SOC = 0
E_FERMI =  7.1370
/

&PARAMETERS
OmegaNum =  601
OmegaMin = -0.5
OmegaMax =  0.5
Nk1 =  201
Nk2 =  201
Nk3 =  201
Eta_Arc = 0.001
/

KCUBE_BULK
  0.00  0.00  0.00   ! Original point for 3D k plane
  1.00  0.00  0.00   ! The first vector to define 3d k space plane
  0.00  1.00  0.00   ! The second vector to define 3d k space plane
  0.00  0.00  1.00   ! The third vector to define 3d k cube
```

### Intrinsic second order Hall effect  ($\chi_{abc}$)

Ref: 10.1103/PhysRevLett.127.277201, 10.1103/PhysRevLett.127.277202

```
&CONTROL
sigma_SOAHC_int_calc = .TRUE. ! static mpi, fixed k-mesh
/

&SYSTEM
SOC = 0
E_FERMI =  7.1370
/

&PARAMETERS
OmegaNum =  601
OmegaMin = -0.5
OmegaMax =  0.5
Nk1 =  2001          !  very high value for fixed k-mesh method
Nk2 =  2001
Nk3 =  2001
Eta_Arc = 0.001
/

KCUBE_BULK
  0.00  0.00  0.00   ! Original point for 3D k plane
  1.00  0.00  0.00   ! The first vector to define 3d k space plane
  0.00  1.00  0.00   ! The second vector to define 3d k space plane
  0.00  0.00  1.00   ! The third vector to define 3d k cube
```

### Intrinsic nonlinear planar Hall effect ($\chi_{abcd}$)

Ref: 10.1103/PhysRevLett.130.126303

### Drude-like nonlinear planar Hall effect ($\chi_{abcd}/\tau^{2}$)

Ref: 10.1103/PhysRevB.108.075155

### Third order Hall effect ($\chi_{abcd}/\tau^{1}$ and $\chi_{abcd}/\tau^{3}$)

Ref: 10.1103/PhysRevB.105.045118

```
&CONTROL
sigma_NPHC_int_calc  = .TRUE.  ! dynamical mpi, auto adapted k-mesh
! sigma_NPHC_tau2_calc = .TRUE.  ! dynamical mpi, auto adapted k-mesh
! sigma_TRAHC_calc     = .TRUE.  ! dynamical mpi, auto adapted k-mesh
/

&SYSTEM
SOC = 0
E_FERMI =  7.1370
/

&PARAMETERS
OmegaNum =  601
OmegaMin = -0.5
OmegaMax =  0.5
Nk1 =  201          !  it will be auto adapted with 5x5x5 local k-mesh
Nk2 =  201
Nk3 =  201
Eta_Arc = 0.001
/

KCUBE_BULK
  0.00  0.00  0.00   ! Original point for 3D k plane
  1.00  0.00  0.00   ! The first vector to define 3d k space plane
  0.00  1.00  0.00   ! The second vector to define 3d k space plane
  0.00  0.00  1.00   ! The third vector to define 3d k cube
```

### Distributions of intrinsic nonlinear planar Hall effect
```
&CONTROL
band_geo_props_kplane_calc = .TRUE.
/

&SYSTEM
SOC = 0             ! without
E_FERMI = 6.3906    ! e-fermi in the hr.dat
/

&PARAMETERS
Eta_arc = 0.02      ! infinite small value, like brodening
E_arc = 0.0         ! energy for calculate Fermi Arc
Nk1 = 301           ! number k points  odd number would be better
Nk2 = 301           ! number k points  odd number would be better
/

KPLANE_BULK         ! unit is the reciprocal lattice vectors
 0.00  0.00  0.00   ! center point for 3D k plane
 1.00  0.00  0.00   ! The first vector to define 3d k space plane
 0.00  1.00  0.00   ! The second vector to define 3d k space plane
```
