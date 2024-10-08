#=
!=======================================================================
!                     MOLECULAR/ATOMIC LINE MODULE
!
! NOTE:
!   The molecular line data file format is that of the Leiden Atomic and
!   Molecular Database (LAMBDA), see
!   http://www.strw.leidenuniv.nl/~moldata/.  The paper belonging to that
!   data is: Schöier, F.L., van der Tak, F.F.S., van Dishoeck E.F., Black,
!   J.H. 2005, A&A 432, 369-379.
!
! The lines_mode tells how the level populations are computed. The 
! sign of the lines_mode has the following meaning:
!
!   lines_mode < 0    means on-the-fly computation (during the 
!                     ray-tracing, i.e. local; can be slow)
!   lines_mode > 0    means computing the populations beforehand, 
!                     storing them in a global array, and only then
!                     doing ray-tracing. This requires substantial
!                     memory.
! 
! The precise value of lines_mode has the following meaning: 
!
!   lines_mode = -1,+1    LTE populations
!   lines_mode = -2,+2    User-defined populations. RADMC-3D call the
!                         subroutine userdef_compute_levelpop() for that.
!                         Rest of line transfer is done by RADMC-3D.
!   lines_mode = -10      User-defined mode II: Line transfer is much
!                         more in the hands of the userdef_module.
!                         Using userdef_general_compute_levelpop().
!   lines_mode = 50       Level populations are read from file
!   lines_mode = -3,3     Large Velocity Gradient method (local non-LTE)
!   lines_mode = -4,4     Optically thin populations (local non-LTE)
!   lines_mode = 20..49   [FUTURE] These will be the non-local non-LTE
!                         modes such as lambda iteration, accelerated
!                         lambda iteration etc.
!
!=======================================================================
=#