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
!
! SCATTERING OFF RANDOMLY ORIENTED NON-SPHERICAL PARTICLES
!
! Now let's look at the scattering. For this we need the scattering Mueller
! matrix Z(theta,phi). Here theta is the degree of directional change in
! radian, with theta=0 meaning forward scattering, theta=pi/2 is
! perpendicular scattering and theta=pi is back-scattering. The phi angle is
! the scattering direction in the (x',y') plane. It defines the "scattering
! plane" spanned by the incoming and outgoing photon direction. For phi=0
! this scattering plane is in the (x',y')=(1,0) direction. In that case,
! the S-vector acts as the rotation axis of the scattering. 
!
! NOTE: Our convention is that the local 2-D (x',y') basis is such that it
! defines a 2-D plane perpendicular to the INCOMING photon. This photon is
! pointing toward the observer, i.e. if we make 3-D vectors of the x'- and
! y'-basis: e_x' and e_y', then e_z'=cross(e_x',e_y') points toward the
! observer and gives the direction of propagation of the incoming photon.
! The S-vector is by definition e_y'. If phi=0 and 0<theta<pi, the
! scattering takes place in the x'-direction, i.e. in 3-D this is the plane
! spanned by (e_x',e_z'), i.e. the plane perpendicular to S==e_y'. For phi=0
! and 0<theta<pi the new photon propagation direction after scattering will
! be in positive x' direction, i.e. n_new = a_z*e_z' + a_x*e_x', with
! a_x>=0.
!
! The Stokes parameters of the scattered flux off a single grain can be
! written in terms of the Stokes parameters of the incoming flux using the
! scattering matrix Z:
! 
!   ( I_scat )            ( Z_11   Z_12   Z_13   Z_14 ) ( I_inc )
!   ( Q_scat )    m_grain ( Z_21   Z_22   Z_23   Z_24 ) ( Q_inc )
!   ( U_scat )  = ------- ( Z_31   Z_32   Z_33   Z_34 ) ( U_inc ) 
!   ( V_scat )      r^2   ( Z_41   Z_42   Z_43   Z_44 ) ( V_inc )
!
! where (I,Q,U,V)_inc = flux_inc and (I,Q,U,V)_scat = flux_scat. The
! scattering matrix is Z is the matrix-valued differential cross section for
! the scattering per unit mass of dust (hence the factor m_grain in the
! equation above). Note that the book of van de Hulst calls this matrix F
! instead of Z. The book by Mishchenko calls it Z for aligned particles, but
! F for randomly oriented particles, because in these two cases they
! define the (x',y')_inc and (x',y')_scat basis vectors differently
! (in the fixed alignment case the particle orientation is fixed in
! the lab frame while in the randomly oriented particle case the 
! x',y' for incoming radiation is oriented w.r.t. the outgoing scattering
! angle). For exact normalizations of Z, see below in the section
! comparing our conventions to standard books.
!
! For randomly oriented particles we need to specify Z only as a function of
! theta. This is, however, only true for randomly oriented particles for
! which each helicity has equal amount of counter-helicity. If we have
! magnetic or gravitational alignment of particles then the full
! Z(theta,phi) must be specified.
!
! The scattering matrix Z(theta) is thus given for a '-basis such that the
! S-vector is perpendicular to the scattering plane. Typically your photon
! initially has some arbitrary perpendicular S-vector. Scattering may then
! take place in any arbitrary scattering plane. If we want to use the
! Z(theta), for that scattering direction, you first have to rotate the
! '-basis to the ''-basis for which the NEW S-vector (i.e. the
! x''-direction) is again perpendicular to the scattering plane. Then we can
! apply the Z(theta), giving us the new (I,Q,U,V). This also automatically
! defines the new orientation of the (x',y') plane which, after the
! scattering event, is of course also rotated in 3-D space. The 3-D version
! of the y' axis stays unaltered during scattering and the 3-D version of
! the x'-axis is rotated around the y'-axis (=S-vector) by an angle theta
! in right-handed direction (if you put the S-vector pointing toward you,
! the x' axis rotates counter-clockwise).
!
! If we want to calculate the scattering source function in a given
! direction of the observer, for a given incident angle of radiation, we
! know precisely the theta and phi angles in advance. We can then indeed
! simply rotate the '-basis, multiply the scattering matrix, and rotate to
! the user-specified S-vector. This is relatively simple, and is done
! by the subroutine polarization_randomorient_scatsource().
! 
! However, if we want to choose a new photon direction and polarization for
! a photon package in a Monte Carlo code, then we do not know the theta and
! phi in advance. In fact, we run into the problem that we cannot a-priori
! rotate the '-basis into the scattering direction.  On the contrary, we
! want to know the scattering probability as a function of theta and phi, so
! we NEED the full Z(theta,phi) rather than just Z(theta). Fortunately we
! can still do this with rotation, because multiplying with Z(theta,phi) is
! the same as (1) first rotating the basis by angle ang=phi (using the
! formulae above), (2) multiplying with the Z(theta) matrix, and (3)
! rotating back to the original basis by choosing ang=-phi. 
!
! So let us assume that we have the full scattering matrix Z(theta,phi). Then
! if we integrate the outgoing I_out(theta,phi) over 4*pi for a given
! unpolarized input intensity I, we should get kappa_scat*I:
!
!      /+1             /2pi                             
!      |   dcos(theta) |   dphi  Z(theta,phi)[1,1] * I   
!      /-1             /0                                
! 
!                               = kappa_scat * I
!
! That is how we define the normalization of the scattering Mueller
! matrix. The Z(theta,phi)[1,1] is thus the angular differential cross
! section of scattering per unit mass of dust. For isotropic scattering it
! would thus be Z11=kappa_scat/(4*pi).
!
! For polarized input intensity this should hold as well, as long as the
! grains are randomly oriented and have no netto helicity (i.e. equal
! left- and right-helicity). 
!
! We can verify this by doing the integrals for the q, u and v components:
!
!                      /+1     /2pi      4                      
! kappa_scat(q,u,v) =  |   dmu |   dphi Sum Z(theta,phi)[1,i] * (1,q,u,v)[i]
!                      /-1     /0       i=1                         
!
! where mu=cos(theta0), q=Q/I, u=U/I, v=V/I (i.e. all values between -1 and 1).
! We can work out this integral a bit deeper, because the integration over
! dphi can be done as follows:
! 
!                      /+1      4                  /2pi         
! kappa_scat(q,u,v) =  |   dmu Sum Z(theta)[1,i] * |   dphi Rotate_phi{(1,q,u,v)}[i]
!                      /-1     i=1                 /0               
!
! i.e. we now have only Z(mu) instead of Z(mu,phi). We have
!
!  Rotate_phi{(1,q,u,v)} = 
!     (1,cos(2*phi)*q+sin(2*phi)*u,-sin(2*phi)*q+cos(2*phi)*u,v)
!
! The integral over 2*pi of the middle two components is 0, so we obtain
!
!                          /+1      4 
! kappa_scat(q,u,v) = 2 pi |   dmu Sum Z(theta)[1,i] * (1,0,0,v)[i]
!                          /-1     i=1
!
!                          /+1      
!                   = 2 pi |   dmu ( Z[1,1] + Z[1,4]*v )
!                          /-1     
!
! As one can see, linearly polarized light has the same scattering cross
! section as unpolarized light. This is because of the random orientation
! of the dust particles, and will NOT be the case when the particles are
! somehow aligned. But note that circularly polarized light might have
! a different scattering cross-section if Z[1,4].ne.0. This means that if
! we have netto helicity of our particles, then we also have dependence
! on the input polarization state. 
!
! For scattering off randomly oriented particles with a plane of symmetry
! (i.e. non-helical) the scattering matrix takes the form (e.g. Mishchenko et
! al. 1996 JQRST 55, 535):
!
!           ( Z_11   Z_12    0      0   )
!           ( Z_12   Z_22    0      0   )
!  Z(mu) =  (  0      0     Z_33   Z_34 )
!           (  0      0    -Z_34   Z_44 )
!
! This has 6 independent matrix elements.
!
! For the somewhat more special case of Mie scattering off spherical
! particles the scattering matrix is even somewhat simpler:
!
!           ( Z_11   Z_12    0      0   )
!           ( Z_12   Z_11    0      0   )
!  Z(mu) =  (  0      0     Z_33   Z_34 )
!           (  0      0    -Z_34   Z_33 )
!
! This has only 4 independent matrix elements! A famous code that can
! compute these, and is freely available on the internet, is the code
! by Bohren & Huffman from the appendix of their book. A version that
! has been improved by Bruce Draine in the early 90s is the bhmie code
! from the scatterlib library:
!
!   http://code.google.com/p/scatterlib/
!
! The callbhmie() program from that bhmie library computes not only
! the scattering and absorption cross sections, but also the Z_11,
! Z_12, Z_33 and Z_34 (which they call S_11, S_12, S_33 and S_34,
! but please see below for the different normalization). 
!
! In this code we allow for the 6-element matrix form, i.e. more general
! than the Mie scattering case, i.e. allowing for prolate or oblate
! particles. 
!
! In any case, for these randomly-oriented or spherical grains we have
! Z_14=0, so that
!
!                           /+1      
!  kappa_scat(q,u,v) = 2 pi |   dmu Z(mu)[1,1]
!                           /-1     
!
! in other words:
!
!  kappa_scat(q,u,v) = kappa_scat
!
! This means that to find out the scattering cross section, we do not need
! any information about the polarization state. This is good, because
! otherwise even the unpolarized Monte Carlo scattering model would be
! wrong, because the polarization state would then affect the further
! scattering behavior. This is, in the case of randomly oriented and/or
! spherical grains not the case. For now we shall use this as a working
! hypothesis. NOTE: It is still not perfectly correct, because while the
! total cross section for scattering is the same, the direction in which the
! photons scatter may be different at the second and later scatterings. But
! it is usually not that problematic.
!
! For Monte Carlo we are concerned with photons (or photon packages). We use
! tau_scat = rho*kappa_scat*pathlength to find out IF a photon scatters.
! 
! ONCE we know that it scatters, the angular probability distribution
! function for the direction of scattering is:
!
!               Z(mu,phi)[1,1] + Z(mu,phi)[1,2] * q + Z(mu,phi)[1,3] * u + Z(mu,phi)[1,4] * v
!  P(mu,phi) = -------------------------------------------------------------------------------
!                                            kappa_scat
!
! where the q, u and v are those of the incoming radiation.
!
! The P(mu,phi) is normalized in (mu,phi) space to unity:
! 
!   /+1     /2pi
!   |   dmu |   dphi P(mu,phi) = 1
!   /-1     /0
!
! so P(mu,phi) is the probability per dmu and per dphi. 
!
! So to find the new direction of scattering we use P(mu,phi) as the
! random-number distribution function and pick a direction. Once we know the
! direction, we still need to figure out what the q,u,v are. For that we
! apply the Z(mu,phi) in the 3-stage way we discussed above.
!
! We can work out P(mu,phi) in terms of Z(mu) matrix elements by performing,
! at a given value of phi, a basis rotation of (x',y') to (x'',y'') such
! that the x''-axis points into the direction of scattering, i.e.
! the e_x'' basis vector equals the e_x' vector rotated by phi in 
! counter-clockwise direction. According to the above rules for
! rotation (with ang=phi) the new q'', u'' and v'' then become:
!
!   q'' =  cos(2*phi)*q + sin(2*phi)*u
!   u'' = -sin(2*phi)*q + cos(2*phi)*u
!
! Since in the (x'',y'') coordinates the scattering is now in the positive
! x''-direction, we have that the probability for scattering is
! P=(Z(mu)[1,1]+Z(mu)[1,2]*q'')/kappa_scat because the other matrix elements
! are zero in this coordinate system.  Inserting the above expression for
! q'' yields:
!
!               Z(mu)[1,1] + Z(mu)[1,2] * ( cos(2*phi)*q + sin(2*phi)*u )
!  P(mu,phi) = -----------------------------------------------------------
!                                    kappa_scat
!
! We can integrate this over phi to get the P(mu):
!
!                   Z(mu)[1,1]
!   P(mu) =  2 pi  ------------
!                   kappa_scat
!
! This shows that, for the probability of scattering with a certain
! angle theta=acos(mu), the input polarization state is irrelevant.
! In other words: the mu-phase function for scattering is entirely
! given by Z(mu)[1,1]. Only the phi-phase function depends on the
! input polarization state.
!
! We can also integrate P(mu,phi) over mu:
!
!
!             Zint[1,1] + Zint[1,2] * ( cos(2*phi)*q + sin(2*phi)*u )
!   P(phi) = ---------------------------------------------------------
!                           2*pi * kappa_scat
!
! with 
!
!                    /+1
!   Zint[1,i] = 2 pi |  dmu Z(mu)[1,i]
!                    /-1
!
! Note that Zint[1,1] == kappa_scat.
!
! Now let's return to the question how to find the new direction of the
! photon using P(mu,phi) as probability distribution function. We do this in
! two stages.
!
!  Stage 1: Find phi using P(phi)
!
!  Stage 2: Compute P(mu|phi), which is the probability function for
!           scattering at an angle mu, given angle phi. P(mu|phi) is:
!
!                                    /+1
!            P(mu|phi) = P(mu,phi) / |  dmu' P(mu',phi)
!                                    /-1 
!
!                        P(mu,phi)
!                      = ---------
!                         P(phi)
!
!            Find mu from P(mu|phi).
!
! Once we found both mu and phi, we can rotate (1,q,u,v) to the new phi, and
! multiply with the scattering matrix to get (i'',q'',u'',v''). However, since
! by definition i'' should become 1, because we have dealt with the
! amplitude of the scattering by using kappa_scat and the P(mu,phi), we
! obtain the new q,u,v as:
!
!   q_new = q''/i'',   u_new = u''/i''  and v_new = v''/i''
!   
! That should be it for the Monte Carlo scattering event. 
!
! NOTE: For randomly-oriented and mirror-symmetric particle ensembles
!       the extinction matrix is scalar, i.e. the extinction coefficient
!       is the same for I,Q,U,V and there is no mixing. In other words:
!       K = diag(1,1,1,1)*kappa. This will NOT be the case for oriented
!       particle ensembles.
! 
!------------------------------------------------------------------------
!
! SCATTERING OFF ALIGNED NON-SPHERICAL PARTICLES
!
! Suppose we have particles which have a plane of symmetry (i.e. no
! helicity) and can be aligned. Take, for instance, prolate or oblate
! ellipsoids. If we have gravity or a magnetic field, then we can align
! these particles with their axis of symmetry along the gravity or 
! B-field direction. Then the scattering Mueller matrix depends on
! three angles instead of just one. 
!
! We rotate the (x',y') plane such that the symmetry axis of the dust
! particle is in the y'-z'-plane. The angle of the particle can be 
! written as xi, such that for xi=0 the symmetry axis of the
! particle is pointing in the y'-direction and for xi=pi/2 it is
! pointing in the z'-direction, and for xi=pi/4 it is in the y'=z'
! direction. The scattering angles remain theta and phi. In this 
! case we can, unfortunately, not reduce angles through symmetries.
! We will have to account for the full Z(xi,theta,phi)-dependence.
! Also: We cannot assume that the matrix has the upper-right and
! lower-left matrix elements zero.
!
! If we wish to include this at some point, we must embed the 
! cross-section routine (e.g. the T-matrix code of Mishchenko)
! into the code, because it would make no sense to precalculate
! and store a 4-dimensional array (3 angles and 1 freq), I think.
!
! It is going to be difficult to implement this into the Monte
! Carlo code, because that involves integrals over the cross
! section. For a "last scattering method" this would, however,
! be no problem.
!
! FOR NOW WE DO NOT INCLUDE THIS POSSIBILITY IN RADMC-3D
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