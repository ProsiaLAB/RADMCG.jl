#=
=======================================================================
!                       LAMBDA OPERATOR MODULE
!
! For the Lambda Iteration scheme (and its sister, the Accelerated
! Lambda Iteration scheme) a Lambda Operator is required.
!
! A Lambda operator is simply an operator (subroutine) that computes
! the mean intensity at some wavelength using the integrals over a large
! number of rays. In contrast to the mean intensity subroutine in the
! montecarlo_module, for the Lambda operator we do not "follow a photon
! along its way" but we instead "integrate the formal transfer equation
! along a ray". If no an-isotropic dust scattering is included in the
! problem, and the scattering source function is known, this gives a 
! more accurate answer. Also, depending on which method of Lambda
! operator you use, you get less noise. Finally, one can use the
! Accelerated Lambda Iteration scheme to speed up convergence of the
! Lambda Iteration.
======================================================================
=#