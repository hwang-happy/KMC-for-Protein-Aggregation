#@ Ace Shen, 8/22/2017
# This is a python script solving diffusion PDE, d phi / d t = D d^2 phi / d x^2, in 1-d using kinetic Monte Carlo method. The data is compared to the analytical solution in the plot. Certainly, increasing number of sample would improve the accuracy of the KMC simulation.



import numpy as np
import matplotlib.pyplot as plt
import random
import scipy
from scipy import special


# Total number of samples
nsample = 50000;

# 1-D lattice

#N = int(input('Enter number of grid points: '))
N = 1000; # Number of grids
L = N-1;   # System size, scaling
h = L/float(N-1) /2.  # Grid spacing, divided by 2 since for even number of moves and starting from origin, you won't stop at, say, the first grid to the right. So we only record those grids that will have particles at the end of simulation, or we allow the particle to move half of the grid spacing for a single step.


#nStep = int(input('Enter number of steps: '))
nStep = 200; # Number of steps

# Constants
# The timescale for a particle to diffuse to the nearest neighbor is tau ~ h^2/D and thus
# the rate for the diffusive process is k = 1/tau

D = 1; # Diffusion coefficient, also used in the analytical solution

kplus = D/h**2; # Rate that the particle moves to the right
kminus = D/h**2; # Rate that the particle moves to the left

# Constructing the cumulative rate vetor
k1 = kplus;
ktot = kplus + kminus;

# Time step

#tau = float(input('Enter time step: '))
tau = 1/(ktot);


# Grid vector
x = np.arange(0,N+1) - N/2
# Vector that records the number of particles at each grid point
xposlist = np.zeros(len(x));


# Initial time
t = 0;

# Initial position
x0 = 0;



## MAIN LOOP ##

for isample in range(1, nsample+1):

    # reset the initial position for each sample
    xposition = x0;
    

    for iStep in range(1,nStep+1):


        # Generating random number for each process
        r = random.uniform(0,1);
    
        # Determine which state the system go
        if( r*ktot <= k1 ):
            xposition = xposition + h;
        
        else:
            xposition = xposition - h;
        


        # Record the position in the position vector for each sample
        if( iStep == nStep ):
            xposlist = xposlist + ( x == xposition)*1;
    
    






#t = t - np.log(random.uniform(0,1)+0.0000001)/ktot;
    
# Total elapsed time
t = nStep *tau;

# Normalize
xposlist = xposlist / nsample;


# Analytical solution

y =  np.exp(-x**2/(4 * D * t))/np.sqrt(4*np.pi*D*t);



# Plotting
plt.figure(1)   # Clear figure 1 window and bring forward
plt.plot(x,xposlist,'o',x,y,'-')
plt.legend(['t='+str(t)])
plt.xlim([-50,50])
plt.xlabel('x')
plt.ylabel('$\phi (x)$',fontsize = 18)
plt.grid(True)
plt.show()



















