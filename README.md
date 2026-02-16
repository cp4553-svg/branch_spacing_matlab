# branch_spacing_matlab
Readme file for matlab codes (migration.m; proliferation.m; with_branches.m)
proliferation.m and migration.m files simulate the behaviors of mesenchymal cells near a source that either enhances proliferation (proliferation.m) or acts as a chemotactic source (migration.m).
with_branches.m simulates the behavior of epithelial and mesenchymal cells in the presence of a chemoattractant source.

Common:
Domain
The system evolves in a square domain with side length Lx = Ly = 20.
Reflective boundary conditions are applied at all walls.

Symbol	Meaning
N	initial number of cells
v0	motility speed
D_r	rotational diffusion
alpha_m	chemotactic strength
k	cell-cell repulsion
k_w	cell-wall repulsion
k_bend	wall stiffness
k_press	wall pressure
k_spring	wall elasticity
d_cutoff	interaction range
gamma1	cell mobility
gamma2	wall mobility

Central Source
A circular region centered at (10,10) with radius R=1 acts as a chemoattractant target and a geometric exclusion zone.
Depending on the simulation, the source either enhances proliferation or acts as a chemoattractant.

File: proliferation.m
The number of cells inside the source-adjacent region N_region controls population growth.
Every cell within N_region proliferates periodically.

File: migration.m
The propulsion direction of cells located in the source-adjacent region aligns toward the source.
The strength of alignment is controlled by alpha_m.

File: with_branches.m
Mesenchymal cells located in the source-adjacent region migrate toward the source.
Interactions between mesenchymal and epithelial cells deform the epithelial branch, which is governed by interaction forces, bending stiffness, pressure, and spring force.
For non-interacting scenario, alpha_m is set to 0.



