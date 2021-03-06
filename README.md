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


** Note for the diffusion1d_aatoa2.py

It's a script solving the reaction-diffusion equation using Gillespie SSA. The chemical reactions, besides diffusion, are<br />
a + a -> a2 with rate constant kp<br />
a2 -> a + a with rate constant km<br />.
Diffusion of a2 particles are ignored for now. 

** This is only for dimer, a2. If we need a3, a4, ... to an where an could be an aggregation of a thousand particles, it seems this algorithm is not applicable. Maybe we could invent a new algorithm to attack this problem.

To be more specific, consider the following aggregation process:<br />
a + a -> a2<br />
a2 + a -> a3<br />
a3 + a -> a4<br />
...<br />
a(n) + a -> a(n+1)<br />

and these linear aggregates could dissolve into shorter aggregates:<br />
a(n) -> a(m) + a(n-m), where 1 < m <n.
For the algorightm I used in diffusion1d_aatoa2.py, it requires 3 layers of if and while loop for 4 kinds of interactions, diffusion to the left, diffusion to the right, a+a->a2 and a2->a+a. If we want to deal with a(n), we would need at least n*(n-1)/2 layers, which is absurd. 




