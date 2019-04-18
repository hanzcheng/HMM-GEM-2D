# HMM-GEM-2D

HMM-GEM scheme for generic meshes in 2D

meshes: Different mesh types are given in meshes_HMM_GEM

Problems: 

1. Diffusion model (found in HMM_diffusion_with_postprocessing)

2. Advection-reaction model for solenoidal fields (found in ELLAM_MMOC_advection_reaction)

3. Miscible flow model (found in HMM-GEM_miscible_flow)

Some references: 

[1] J. Droniou, R. Eymard, T. Gallouet, and R. Herbin. A unified approach to mimetic finite
difference, hybrid finite volume and mixed finite volume methods. Math. Models Methods
Appl. Sci., 20(2):265-295, 2010. https://arxiv.org/abs/0812.2097

[2] H. M. Cheng and J. Droniou. An HMM-ELLAM scheme on generic
polygonal meshes for miscible incompressible flows in porous media.
Journal of Petroleum Science and Engineering, 172:707-723, 2019. https://arxiv.org/abs/1706.02452

[3] H. M. Cheng, J. Droniou, and K.-N. Le. A combined GDM-ELLAM-MMOC scheme for
advection dominated PDEs. ArXiv e-prints, May 2018. https://arxiv.org/abs/1805.05585

Notes: 

1. For a detailed discussion of the HMM, see [1].

2. The full HMM--ELLAM scheme for the miscible flow model can be found in ([2] Chapter 2).

3. The post-processing of the HMM solutions used for velocity reconstructions can be found in ([2] Section 2.3)

4. Discussion on how to perform characteristic tracking using the reconstructed RT0 velocities is found in ([2] Section 2.4). 

5. The HMM--GEM scheme for the miscible flow model is presented in ([3] Chapter 5). 

6. The combined ELLAM-MMOC scheme for advection PDEs is presented in ([3] Chapter 4).

7. The local volume adjustments can be found in ([3] Section 2.5). 

