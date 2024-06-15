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
drude_weight_calc = .TRUE.  ! static mpi, fixed k-mesh
```

### Intrinsic second order Hall effect  ($\chi_{abc}$)
```
sigma_SOAHC_int_calc = .TRUE. ! static mpi, fixed k-mesh
```

Ref: 10.1103/PhysRevLett.127.277201, 10.1103/PhysRevLett.127.277202

### Intrinsic nonlinear planar Hall effect ($\chi_{abcd}$)
```
sigma_NPHC_int_calc = .TRUE. ! dynamical mpi, auto adapted k-mesh
```

Ref: 10.1103/PhysRevLett.130.126303

### Drude-like nonlinear planar Hall effect ($\chi_{abcd}/\tau^{2}$)
```
sigma_NPHC_tau2_calc = .TRUE. ! dynamical mpi, auto adapted k-mesh
```

Ref: 10.1103/PhysRevB.108.075155

### Third order Hall effect ($\chi_{abcd}/\tau^{1}$ and $\chi_{abcd}/\tau^{3}$)
```
sigma_TRAHC_calc = .TRUE. ! dynamical mpi, auto adapted k-mesh
```

Ref: 10.1103/PhysRevB.105.045118

### Distributions of intrinsic nonlinear planar Hall effect
```
band_geo_props_kplane_calc = .TRUE.
```
