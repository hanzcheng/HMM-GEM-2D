ELLAM or MMOC scheme for an advection-reaction problem:

Main file: 

ELLAM_MMOC_advection_reaction - solves an advection-reaction problem (with given velocity) via an ELLAM or MMOC scheme

Initialisation: 
1. gravity_centers - computes the cell center of gravity
2. edge_points- gives n equally spaced points along the edges of each cell
3. nPtsPerEdge - determines how many points to track along the edges of each cell
4. initCond - specifies the initial condition
5. compute_qPoints - computes "quadrature points", to be used when computing a reference solution

Files related to characteristic tracking: 

1. FO_Euler_exactVel - performs the characteristic tracking by using a forward Euler scheme and using micro timesteps
2. complete_tracking_exact_vel - performs the full tracking (i.e. loops and perform FO_Euler_exactVel until the time step is exhausted)
3. adjust_volumes_for_local_mass_divFree - local volume correction algorithm for divergence free velocity fields

Processing:
1. compArea_cell - computes area of intersection bet. 2 polygons
2. compute_refSol - computes a reference solution (only works when there are no source terms)
3. compute_errors - measures L1 and L2 errors
4. PolygonClip - clips out the polygonal region of intersection bet. 2 polygons (source: http://www.cs.man.ac.uk/~toby/gpc/)
5. write_solution_vtk - writes a vtk file for visualising the solution profile
6. write_solution_ensight - writes an ensight file for visualising the solution profile (note: Here, the results are written in a different folder, in this case, specified to be ensight/Un. Users may opt to change where the results are written by changing lines 121-134.) 

A sample test result is given in the folder ensight. This is obtained by running ELLAM_MMOC_advection_reaction as is. Of course, the mesh, time step, source term, velocity field, and initial conditions may be modified to run different test cases on different types of meshes. Convergence can also be tested by refining the mesh and reducing the time step.