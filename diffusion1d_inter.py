#@ Ace Shen, 8/29/2017
# 
# Gillespie Algorithm, solving reaction-diffusion equation.
# Ver. 2
# Including nucleation, a + a -> a2

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy
import mpl_toolkits
from matplotlib import cm           
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import special
from mpl_toolkits.mplot3d import Axes3D

# Total number of samples
nsample = 1;

# Total number of particles
N = 50000;

# Total length of the system
L = 1; # mm

# Number of compartments
compartment = 51;

# Compartment length
h = L/compartment;


# Rate constant
kp = 0.01; # a + a -> a2
km = 0.001; # a2 -> a + a

# Diffusion rate constant
d = 1/(h*h);


# Initialize the number of molecules in each compartment
a = np.zeros(compartment+1);
a2 = np.zeros(compartment+1);


# Initial condition, the initial number of molecule in each compartment
#a[16] = N/2;
#a[17] = N/2;


a[20:30+1] = N/10;


# Initial time
t = 0;

# Number of steps
nStep = 10000000;

for iStep in range(1,nStep+1):
    
    # Generate two random numbers uniformly distributed in (0,1)
    r1 = random.uniform(0,1)+0.000001;
    r2 = random.uniform(0,1)+0.000001;
    
    
    # The propensity function
    alpha = 2 * d * sum(a) - d * a[1] - d * a[-1] + kp * sum( a*(a-1) ) + km * sum(a2);
    
    # The time interval for which the next reaction (here diffusion process) happens
    t = t + np.log(1/r1)/alpha;
    
    # Determine where does the precess happen
    
    c = 0;
    partialsum = 0;
    
    while( partialsum <= alpha * r2 and c < compartment - 1 ):
        c = c + 1;
        partialsum = partialsum + a[c] * d;
    
    
    # a particle in c compartment moves to the next right compartment
    if( r2 * alpha < partialsum ):
        a[c] = a[c] - 1;
        a[c+1] = a[c+1] + 1;
                
    else:
        c = 1;
        
        while( partialsum <= alpha * r2 and c < compartment ):
            c = c + 1;
            partialsum = partialsum + a[c] * d;
        
        if( r2 * alpha < partialsum ):
            a[c] = a[c] - 1;
            a[c-1] = a[c-1] + 1;
        
        else:
            c = 0;
            
            while( partialsum <= alpha * r2 and c < compartment ):
                c = c + 1;
                partialsum = partialsum + a[c] * (a[c]-1) * kp;
                
            if( r2 * alpha < partialsum ):
                a[c] = a[c] - 2;
                a2[c] = a2[c] + 1;
            
            else:
                c = 0;
                
                while( partialsum <= alpha * r2 and c < compartment ):
                    c = c + 1;
                    partialsum = partialsum + a2[c] * km;
                
                a[c] = a[c] + 2;
                a2[c] = a2[c] - 1;
                    
                    
                    
        
        
        
        
        
        
        
        
        
        
        
print('t = ',str(t))



xa = np.arange(1,compartment + 1);


#y =  np.exp(-((xa-26)*h)**2/(4 * 1 * t))/np.sqrt(4*np.pi*1*t)/compartment;

plt.figure(1);
plt.plot(xa,a[1:compartment+1],'o',xa, a2[1:compartment+1],'o')
plt.grid(True)
plt.xlabel('x', fontsize = 16)
plt.ylabel('Number of molecules', fontsize = 16);
plt.tick_params(axis='both', labelsize = 13)
plt.legend(['A','AA'], fontsize = 12)
plt.title('Diffusion in 1d, t='+str(round(t,3)),fontsize = 18)
plt.savefig("1d.pdf",format="pdf",dpi=1200)

#plt.axis([0, 40, 0, 50])

plt.show()














