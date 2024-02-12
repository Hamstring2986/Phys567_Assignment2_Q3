"""
This code simmulates the cross-sections of two vortex rings that leapfrogs.
The vortices coss-section are modeled using line vortex where the velocity field id given by
u = k/r in the azimuthal direction.

@author: Marc-Antoine Leclerc
@collab: List the names of your collaborators if any
February 12th, 2024
"""
#######################################################
#            Note to the grader
########################################################
# I plotted the vortex position without deleting them every iteration so that the overall motion of the vortices is well vizualized.
# If you run this code, you will see the two top vortex oscillate up and down which is what we are supposed to see for leapfrogging.
# Same thing for the two bottom vortex.
# I have included the math in my pdf submission 



####################################################
# Libraries, Constants and Functions
####################################################


import numpy as np
import matplotlib.pyplot as pl

def compute_velocity_field(l):
    """
    Compute and update the velocity field for the l'th vortex

    Parameters:
        l (int): Vortex index.
    Return
        None
    """
   
    # Iterating through the grid
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):

            # These conditional statements are to avoid any singularities
            if X[i, j] == x_v[l] and Y[i, j] == y_v[l]:
                continue

            elif X[i, j] == x_v[l]:
                prefactor = k_v[l] / np.sqrt((Y[i,j] - y_v[l])**2)

            elif Y[i, j] == y_v[l]:
                prefactor = k_v[l] / np.sqrt((X[i,j] - x_v[l])**2)

            else:
                prefactor = k_v[l] / np.sqrt((X[i,j] - x_v[l])**2 + (Y[i,j] - y_v[l])**2)

            # Computing the velocity field using the formula provided in the submitted pdf
            vel_x[i,j] += prefactor*(-(Y[i,j] - y_v[l]))
            vel_y[i,j] += prefactor*(X[i,j] - x_v[l])

# Defining the timestep and the number of timesteps
dt = 0.05
Nsteps = 30
r_mask = 0.1

###############################################################
#Setting up initial conditions and initial plot
###############################################################

# Vortex rings initial position
y_v = np.array([2.0, -2.0,2.0,-2.0], dtype=float)
x_v = np.array([-2.0,-2.0,-2.5,-2.5],dtype=float)

# Circulation = 2*pi*k
# Vortex constant k
k_v = np.array([2.0,-2.0,2.0,-2.0],dtype=float) 

# Setting up the plot
pl.ion()
fig, ax = pl.subplots(1,1)

# mark the positions of vortices
color = ["ro", "bo", "go","ko"]
for i in range(len(color)):
    ax.plot(x_v[i], y_v[i], color[i], markersize=5) 

# Setting up the grid
ngrid = 4
Y, X = np.mgrid[-ngrid:ngrid:100j, -ngrid:ngrid:100j] 

# Initiliazing the x and y velocities to 0
vel_x = np.zeros(np.shape(Y),dtype=float) 
vel_y = np.zeros(np.shape(Y),dtype=float) 

#looping over each vortex to initialize the velocity field
for i in range(len(x_v)): 
    compute_velocity_field(i)

# Creating the masked fields
mask = np.zeros_like(X, dtype=bool)
for i in range(len(x_v)):
    mask |= np.sqrt((X - x_v[i])**2 + (Y - y_v[i])**2) <= r_mask
    colors = ["ro", "bo", "go","ko"]
    ax.plot(x_v[i], y_v[i], color[i], markersize=5)

vel_x_masked = np.where(mask, float("NaN"), vel_x)
vel_y_masked = np.where(mask, float("NaN"), vel_y)

# set up the boundaries of the simulation box
ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

# initial plot of the streamlines
ax.streamplot(X, Y, vel_x_masked, vel_y_masked, density=[1, 1]) 
fig.canvas.draw()
pl.pause(10)



#####################################################
#                   Evolution
#####################################################

count = 0
while count < Nsteps:
    
    ### Compute and update advection velocity

    # Re-initializing the total velocity field for computation of total velocity field
    
    vel_x = np.zeros(np.shape(Y),dtype=float) 
    vel_y = np.zeros(np.shape(Y),dtype=float) 

    # Computing total velocity field
    for i in range(len(x_v)):
        compute_velocity_field(i)
    
    # Updating the positions of vortices
    for i in range(len(x_v)):

        # Get the indices of the grid position corresponding to the closest position of the vortices
        distances = np.sqrt((X - x_v[i])**2 + (Y - y_v[i])**2)
        min_index = np.unravel_index(np.argmin(distances), distances.shape)
        closest_position = (Y[min_index], X[min_index])

        # New position
        x_v[i] += dt * vel_x[min_index]
        y_v[i] += dt * vel_y[min_index]

    ### update plot

    # Re-initializing the total velocity field 
    vel_x = np.zeros(np.shape(Y),dtype=float) 
    vel_y = np.zeros(np.shape(Y),dtype=float)

    # Re-calculating the total velocity field for plotting
    mask = np.zeros_like(X, dtype=bool)
    for i in range(len(x_v)):
    
        compute_velocity_field(i)

        # Calculating the mask
        mask |= np.sqrt((X - x_v[i])**2 + (Y - y_v[i])**2) <= r_mask
        colors = ["ro", "bo", "go","ko"]
        ax.plot(x_v[i], y_v[i], color[i], markersize=5)

    # Getting the masked field
    vel_x_masked = np.where(mask, float("NaN"), vel_x)
    vel_y_masked = np.where(mask, float("NaN"), vel_y)
    
    # the following lines clear out the previous streamlines
    for collection in ax.collections:
        collection.remove()
    for patch in ax.patches:
        patch.remove()

    # to plot
    ax.streamplot(X, Y, vel_x_masked, vel_y_masked, density=[1, 1])
    fig.canvas.draw()
    pl.pause(0.1)

    count += 1
