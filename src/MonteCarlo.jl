"""
!-------------------------------------------------------------------------
!               MODULE FOR MONTE CARLO RADIATIVE TRANSFER
!
! This module does radiative transfer in dust continuum. The module
! can do the following things:
!
!  1) Thermal Monte Carlo according to the Bjorkman & Wood method 
!     for computing the dust temperatures
!
!  2) Monochromatic Monte Carlo for computing the scattering source 
!     function or for computing the mean intensity.
!
! The dust scattering can be done to various degrees of realism. This is
! regulated by the scattering_mode integer (see dust_module.f90). The
! various possible values of scattering_mode are:
!
!     =0 Dust scattering not included
!
!     =1 Isotropic scattering approximation
!
!     =2 Include anisotropic scattering with Henyey-Greenstein function
!
!     =3 Include anisotropic scattering with full phase function, based
!        on the full scattering matrix information
!
!     =4 Include scattering with full phase function and polarization 
!        (for randomly oriented particles), but only in the scattering 
!        source function (i.e. the Monte Carlo photon packages remain
!        unpolarized, and the polarization of the scattering source
!        function only takes into account the last scattering before
!        observation)
!
!     =5 Include scattering with full phase function and polarization 
!        (for randomly oriented particles), full treatment.
!
!-------------------------------------------------------------------------
"""