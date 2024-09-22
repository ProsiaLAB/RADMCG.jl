"""
!========================================================================
!                          POLARIZATION MODULE
!
!                Based on subroutines kindly provided by 
!                Michiel Min, Aug 2009 (Thanks Michiel!!)
!
!------------------------------------------------------------------------
!
! DEFINITIONS OF THE STOKES VECTOR COMPONENTS
!
! The Stokes vector is (I,Q,U,V) or, dependent on how you look at it,
! (E,Q,U,V). This vector has meaning in a plane perpendicular to the
! direction of motion of the photon or light ray. There is a degree of
! freedom of chosing the angle of orientation of the (x',y') coordinate
! system in that plane and its right/left-handedness. It is important
! to define exactly how the coordinate system is set up, and how the
! Q, U and V Stokes components are defined with respect to these.
!
! CONVENTIONS: 
! There are different conventions for how to set up the coordinate system
! and define the Stokes vectors. Our definition follows the IAU 1974 
! definition as described in:
!
!   Hamaker & Bregman (1996) A&AS 117, pp.161
!
! In this convention the x' axis points to the north on the sky, while
! the y' axis points to the east on the sky. The z' axis points to the
! observer. This coordinate system is positively right-handed. The 
! radiation moves toward positive z'. Angles in the (x',y') plane
! are measured counter-clockwise (angle=0 means positive x' direction,
! angle=pi/2 means positive y' direction).
!
! In the following we will (completely consistent still with the IAU
! definitions above) define "up" to be positive y' (=east) and "right" to be
! positive x' (=north). So, the (x',y') coordinates are in a plane
! perpendicular to the photon propagation, and oriented as seen by the
! observer of that photon. So the direction of propagation is toward you,
! while y' points up and x' points to the right, just as one would normally
! orient it. Compared to the IAU definition we will point our head in east
! direction, so to speak. 
!
! Note that this is fully equivalent to adjusting the IAU definition to have
! x' pointing west and y' pointing north, which is perhaps a bit more
! intuitive, since most images in the literature have this orientation.  So
! for convenience of communication, let us simply adjust the IAU 1974
! definition to have positive $x'$ ("right") pointing west and positive $y'$
! ("up") pointing north. It will have no further consequences for the
! definitions and internal workings of RADMC-3D because RADMC-3D does not
! know what ``north'' and ``east'' are.
!
! The (Q,U) definition (linear polarization) is such that a linearly
! polarized ray with Q=+I,U=V=0 has the electric field in the (x',y')=(1,0)
! direction, while Q=-I,U=V=0 has the electric field in the (x',y')=(0,1)
! direction. If we have Q=0, U=+I, V=0 then the E-field points in the x'=y'
! direction, while Q=0, U=-I, V=0 the E-field points in the x'=-y'
! direction.
!
! The (V) definition (circular polarization) is such that (quoting directly
! from the Hamaker & Bregman paper): "For *right*-handed circularly
! polarized radiation, the position angle of the electric vector at any
! point increases with time; this implies that the y' component of the field
! lags the x' component. Also the electric vectors along the line of sight
! form a *left*-handed screw. The Stokes V is positive for right-handed
! circular polarization."
!
! We can put these definitions into the standard formulae:
!
!   Q = I*cos(2*beta)*cos(2*chi)
!   U = I*cos(2*beta)*sin(2*chi)
!   V = I*sin(2*beta)
!
! The angle chi is the angle of the E-field in the (x',y') coordinates,
! measured counter-clockwise from x' (consistent with our definition
! of angles). Example: chi = 45 deg = pi/4, then cos(2*chi)=0 and
! sin(2*chi)=1, meaning that Q=0 and U/I=+1. Indeed this is consistent 
! with the above definition that U/I=+1 is E_x'=E_y'. 
!
! The angle 2*beta is the phase difference between the y'-component of the
! E-field and the x'-component of the E-field such that for 0<beta<pi/2 the
! E-field rotates in a counter-clockwise sense. In other words: the y'-wave
! lags 2*beta behind the x' wave. Example: if we have beta=pi/4, i.e.
! 2*beta=pi/2, then cos(2*beta)=0 and sin(2*beta)=1, so we have Q=U=0 and
! V/I=+1. This corresponds to the y' wave being lagged pi/2 behind the x'
! wave, meaning that we have a counter-clockwise rotation. If we use the
! right-hand-rule and point the thumb into the direction of propagation
! (toward us) then the fingers indeed point in counter-rotating direction,
! meaning that V/I=+1 is righthanded polarized radiation.
!
! In terms of the electric fields:
!
!   E_x'(t) = E_h cos(wt-Delta_h)
!   E_y'(t) = E_v cos(wt-Delta_v)
!
! (with E_h>0 and E_v>0 and Delta_h,v are the phase lags of the components
! with respect to some arbitrary phase) we can write the Stokes components
! as:
!
!   I = E_h^2 + E_v^2
!   Q = E_h^2 - E_v^2
!   U = 2 E_h E_v cos(Delta)
!   V = 2 E_h E_v sin(Delta)
!
! with 
!
!   Delta = Delta_v - Delta_h = 2*beta
!
! NOTE: The IAU 1974 definition is different from the definitions used
!       in the Planck mission, for instance. So be careful. There is
!       something said about this on the website of the healpix software:
!       http://healpix.jpl.nasa.gov/html/intronode12.htm
!
!       The IAU 1974 definition is also different from the Mishchenko
!       book and papers (see below).
!
! If we want to use these definitions for a photon moving somewhere in
! space (where "north" and "east" are not defined), all these definition
! are therefore only meaningful if we define the orientation of the (x',y')
! basis. To do this we introduce a unit vector (Sx,Sy,Sz) that is
! perpendicular to the direction of motion and is, by definition, pointing
! in the (x',y')=(0,1) direction. In other words: the S-vector points into
! the negative-Q direction.
!
! We can transform from a '-basis to a ''-basis by rotating the S-vector
! counter-clockwise (as seen by the observer watching the radiation) by an
! angle ang. Any vector (x',y') in the '-basis will become a vector
! (x'',y'') in a ''-basis, given by the transformation:
!
!  (x'') = ( cos ang   sin ang ) (x')
!  (y'') = ( -sin ang  cos ang ) (y')
!
! NOTE: We choose (x',y') to be the usual counter-clockwise basis for the
! observer seeing the radiation. Rotating the basis in counter-clockwise
! direction means rotating the vector in that basis in clockwise direction,
! hence the sign convention in the matrix. 
!
! If we have (I,Q,U,V) in the '-basis (which we might have written as
! (I',Q',U',V') but by convention we drop the '), the (I'',Q'',U'',V'') in
! the ''-basis becomes
!
!  (I'')   (   1        0            0         0   ) (I)
!  (Q'') = (   0    cos(2*ang)  sin(2*ang)     0   ) (Q)
!  (U'')   (   0   -sin(2*ang)  cos(2*ang)     0   ) (U)
!  (V'')   (   0        0            0         1   ) (V)
!
! This form of matrix and its sign convention can be understood in the
! following way: The Q and U can be seen as I*cos(2*beta)*UnitVector, where
! UnitVector can be expressed as a complex number exp(2*i*chi). However, the
! actual vector pointing in the major-axis direction is exp(i*chi) in the
! complex plane. The transformed complex number (rotating clockwise! by ang)
! means multiplying that with exp(-i*ang).  So exp(i*chi'') =
! exp(-i*ang)*exp(i*chi) = exp(i*(chi-ang)). So chi''=chi-ang (not
! surprisingly!). For UnitVector'' ==exp(i*chi'') this means UnitVector'' =
! exp(-2*i*ang)*UnitVector.  In other words: it is like rotating with angle
! 2*ang.
!
!------------------------------------------------------------------------
"""





mutable struct Photon
    E::Float64
    Q::Float64
    U::Float64
    V::Float64
    n::Array{Float64,1}{3}
    k::Array{Float64,1}{3}
end