HMM for the pure diffusion problem

Main files: 

1. convergence_HMM_diffusion - solves the diffusion problem on a sequence of meshes, and measures the rates of convergence in the L2 and Linf norms

2. HMM_diffusion_with_post_processing - solves the diffusion problem on a given mesh, and post-processes the solution to obtain an RT0 velocity on a submesh

Initialisation: 

1. gravity_centers - computes the cell center of gravity
2. setBC - determines (gives a tag) whether to impose Neumann or Dirichlet BC on a given boundary edge
3. ue and derivu - exact solution u and its derivatives
4. lambda - diffusion matrix
5. source_term - gives the source term

HMM: 

1. local_flux_matrix - HMM matrix  
2. Neumann_BC - boundary condition to be imposed on a boundary edge tagged as Neumman

post-processing:

A. velocity reconstruction:

	1. kReconstruction - reconstructs an RT0 velocity (especially useful if an H-div velocity field would be needed for later computations)
	2. interiorFluxes - interior fluxes (KR method)
	3. interiorFluxes_consistent - interior fluxes (C method)
	4. interiorFluxes_aux - interior fluxes (A method)

B. visualisation and error measurement

	1. compute_errors - computes the error in the gradient
	2. create_submesh - creates a triangulation of each cell (the RT velocity fields generated in kReconstruction corresponds to the cells (triangles) in this submesh)
	3. write_solution_vtk - writes a vtk file for visualising the solution profile

