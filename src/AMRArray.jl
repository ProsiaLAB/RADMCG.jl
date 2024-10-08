#=
==========================================================================
           Module for tracking straight line through AMR grid
              and for other formal radiative transfer stuff
==========================================================================

 Before using this AMRRay module you must have set up a grid with the AMR
 module. This is done by calling amr_initialize() with the appropriate
 arguments. You can also call (one of the?) reading routines of the AMR
 library so that the initialization is done for you. But it is important to
 set nrbranchesmax to the maximum number of branches you think you need.
 This is important because the AMRRay module will allocate some data 
 arrays for internal use, and it has to know how many branches there will
 be at maximum. If you set up a grid with, say, 1000 branches, but you
 wish to be able to add branches/leafs later by adaptive grid refinement,
 then you must specify nrbranchesmax to the max number you think you need.
 But do not take it TOO large, because this will waste computer memory.

 Now make sure that the amr_coordsystem value is set to any value between:
     0 <= amr_coordsystem < 100    = Cartesian
   100 <= amr_coordsystem < 200    = Spherical
   200 <= amr_coordsystem < 300    = Cylindrical [not yet implemented]
 The reason why there is this funny range of values all meaning the same
 coordinate system is because the author of this module also made a 
 hydrodynamic package which has more closely specified kinds of coordinate
 systems, and we may in the future have compatibility. So for now you can
 simply take 0, 100 or 200 for cartesian, spherical and cylindrical. 
 Internally in the amrray module this will be translated into 0, 1 or 2
 in the amrray_icoord variable. 

 After setting up the AMR grid, the next thing to do is to call the routine
 amrray_initialize(). This will set up all precalculated stuff required for
 the ray-tracing (especially in spherical coordinates this is necessary). 
 And it will perform a series of checks on the grid setup to make sure
 that it is conform what is expected for the AMRRay module. Also this is
 mainly for the spherical coordinates.

 Now you can start using the module. Once you want to end the program,
 you should call amrray_cleanup() to deallocate all the allocated arrays.


 NOTE: The xc values of each branch MUST be exactly in the middle of
       the cell: a%xc(idir) = 0.5d0 * ( a%xi(1,idir) + a%xi(2,idir) ).
       This is default in the AMR module, but not required in the AMR
       module. But it is strictly required here!
--------------------------------------------------------------------------
=#