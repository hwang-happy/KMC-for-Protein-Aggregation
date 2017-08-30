# KMC-for-Protein-Aggregation

Kinetic Monte Carlo <br />





![alt text](https://img.memecdn.com/kiss-my-ass_o_452063.webp)

Protein Aggregation <br />

Consider in a tube or in a cell, there are whole bunch of proteins floating in it. These proteins has tendency to bind with 
each other and thus they could form protein aggregation (PA).  


** The diffustion equation reads <br />
d phi / dt = D^2 d^2 phi / dx^2. 

==== Old ==================
Consider a particle is placed on a lattice with N grids. The particle can either go right or left with half of a grid spacing h with the same rate k. 
The relation between the transition rate and the diffusion coefficient is<br />
k = D / h^2.

Now let each sample runs for a while and then plot the number of particles on the grid. This could compare to the analytical solution, which read as probability distribution in this context.<br />
==== Old ====================



** Compartment-based approach to diffusion

To simulate diffusion, processing of information from different locations is needed, which can be realized by either setting up grids or making compartments. The basic idea of the compartment-based approach relies on the "hopping" of a single particle from one compartment to another. If we view the hopping process as a "chemical reaction" between two compartments, then the Gillespie SSA is applicable to diffusion, with the rate constant being D/h^2, where h is the length of the compartment.


** Gillespie SSA

The idea of Gillespie algorithm is similar to KMC. Taking diffusion for example, the "chemical raction" describing the diffusion process is A_i -> A_(i+1) with a rate constant d, where the subscript i indicates different compartments. Now, we initialize the system to be in A_i, then through the same process as KMC, a random number would "chose" a particular i, and thus chose where a diffusion process happens.


** Note for the diffusion1d_interaction.py


