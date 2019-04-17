HMM--GEM scheme for the miscible flow problem:

Main files: 
1.ELLAM_MMOC_efficient_tracking_mixed_local_mass_cons - for GEM (with the option of switching to ELLAM)
2.ELLAM_pure_MMOC - for MMOC

Initialisation: 
1. create_submesh - creates a triangulation of each cell so that RT velocity fields can be reconstructed on each sub-triangle
2. gravity_centers - computes the cell center of gravity
3. classifyPtsandRegions - determines which points/ regions correspond to ELLAM/MMOC
4. edge_centers and edge_points- gives the midpoints of the edges of each cell;gives n equally spaced points along the edges of each cell
5. getMbSub - computes the slope and intercept of the lines that determine the edges of the cells (to be used to determine possible exit points after tracking)
6. locateWells - determines which cells the injection and production wells are on
7. nPtsPerEdge - determines how many points to track along the edges of each cell

HMM:
1. local_flux_matrix - generates the HMM matrix

Diffusion tensor and Darcy velocity:
1. difftens - creates the standard diffusion tensor
2. difftens_Mod - creates a modified diffusion tensor (corresponding to vanishing diffusion test case)
3. kReconstruction - reconstructs darcy velocity
4. interiorFluxes - interior fluxes (KR method)
5. interiorFluxes_consistent - interior fluxes (C method)
6. interiorFluxes_aux - interior fluxes (A method)
7. seeDarcyU - plots streamlines of the darcy velocity

 
Files related to characteristic tracking: 

1. charTracking1 - performs the characteristic tracking on a given triangle for a fixed period of time (until the point exits)
2. tExit - determines the time it takes for a point to exit a given triangle
3. detTri - determines the (next) triangle for a point to track into
4. complete_tracking_efficient - performs the full tracking (i.e. loops and perform charTracking1 until the time step is exhausted)
5. adjust_volumes_for_local_mass - local volume correction algorithm 

Processing:
1. compute_area - computes area of a polygonal region
2. compAreaF - computes area of intersection bet. 2 polygons
3. PolygonClip - clips out the polygonal region of intersection bet. 2 polygons (source: http://www.cs.man.ac.uk/~toby/gpc/)
4. write_solution_vtk - writes a vtk file for visualising the solution profile
