elASPECT - elasticity-based ASPECT
==================================

About
-----

elASPECT is mainly designed for numerical simulation of lithospheric scale deformation.
It is born out of the well-known mantle convection simulator ASPECT 
(https://github.com/geodynamics/aspect), but replaces the Stokes equation by the 
elasticity-based equilibrium equation in order to capture the elastic response under
tectonic loads more precisely.


Installation instructions
-------------------------

elASPECT is configured using CMake and has the following requirements:
- CMake 3.1.0 or newer
- GCC with C++14 support
- [deal.II](https://github.com/dealii/dealii) 9.4.1 configured with:
  - MPI, Trilinos, p4est (required)
  - BLAS/LAPACK, zlib (strongly recommended)
  - HDF5 (optional)
