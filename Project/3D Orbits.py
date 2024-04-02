#!/usr/bin/env python
# coding: utf-8

# # Clohessy-Wiltshire-Hill's Equations

# ### For a circular orbit without an external force acting in the Hill ($\mathcal{H}$) frame, the linearized equations are:
# 
# 
# $$ \ddot{x} - 2n\dot{y} - 3n^2x = 0$$ 
# $$ \ddot{y} + 2n\dot{x} = 0$$ 
# $$ \ddot{z} + n^2z = 0$$
# 
# * $n$ is the mean motion of the reference orbit;
# * In the $\mathcal{H}$ frame, the $(x,y,z)$ described by the differential equations correspond to the directions $(\hat e_{r}, \hat e_{\theta}, \hat e_{z})$

# In[1]:


#importing necessary libraries for calculations and simulations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D


# ## Solving the equations numerically

# Defining $\dot{x}=v_x$, $\dot{y}= v_y$ and $\dot{z}=v_z$, we get six first order ODE's:
# 
# Then
# 
# $$\vec{S} = \begin{bmatrix} x \\ v_x \\ y \\ v_y \\ z \\ v_z \end{bmatrix} \hspace{10mm} \implies \hspace{10mm} \frac{d\vec{S}}{dt} = \begin{bmatrix} \dot{x}\\ \dot{v_x} \\ \dot{y} \\ \dot{v_y} \\ \dot{z} \\ \dot{v_z} \end{bmatrix} =  \begin{bmatrix} v_x\\ 2nv_y+3n^2x \\ v_y \\ -2nv_x \\ v_z \\ -n^2z \end{bmatrix}$$

# ## 3D visualization

# In[2]:


def dSdt_linear(t, S):
    x, v_x, y, v_y, z, v_z = S
    dSdt = [v_x, 
            2*n*v_y + 3*n**2*x, 
            v_y, 
            -2*n*v_x, 
            v_z, 
            -n**2*z]
    
    return dSdt

def dSdt_nonlinear(t, S): #with j2 term already #retirar j2
    x, v_x, y, v_y, z, v_z = S
    J2 = 1*10**-3
    mu = 1
    r0 = 1
    dSdt = [v_x, 
           2*n*v_y+n**2*x-mu*(r0+x)/((r0+x)**2+y**2+z**2)**(3/2),#+mu/r0**2+J2*x/(x**2+y**2+z**2)**(7/2)*(6*z**2-3/2*(x**2+y**2)),
            v_y,
            - 2 * n * v_x + n**2 * y - mu*y/((r0+x)**2+y**2+z**2)**(3/2),# + J2*y/(x**2+y**2+z**2)**(7/2)*(6*z**2- 3/2*(x**2+y**2)),
            v_z,
            -mu*z/((r0+x)**2+y**2+z**2)**(3/2)]# + J2*z/(x**2+y**2+z**2)**(7/2)*(3*z**2- 9/2*(x**2+y**2))]
            
    return dSdt

# def dSdt_J2(t, S): #ignore this
#     x, v_x, y, v_y, z, v_z = S
    
#     a=0.2
    
#     r = np.sqrt(x**2+y**2+z**2)
#     v_r = (v_x*x+v_y*y+v_z*z)/np.sqrt(x**2+y**2+z**2)
#     theta = np.arccos(z/np.sqrt(x**2+y**2+z**2))
#     v_theta = (z*(v_x*x+v_y*y)-v_z*(x**2+y**2))/(np.sqrt(x**2+y**2)*(x**2+y**2+z**2))
#     phi = np.arctan(y/x)
#     v_phi = (v_y*x-v_x*y)/(x**2+y**2)
    
#     dSdt = [v_r,
#            r*v_theta**2 + r*(np.sin(theta))**2*v_phi**2-1/r**2+3/2*a/r**4*(3*(np.cos(theta))**2-1),
#            v_theta,
#            np.sin(theta)*np.cos(theta)*v_phi**2-2*v_r/r*v_theta+3*a/r**5*np.sin(theta)*np.cos(theta),
#            v_phi,
#            -2*v_r/r*v_phi-2*np.cos(theta)/np.sin(theta)*v_phi*v_theta]
    
    return dSdt

def hill_reference(dSdt, t):
    S0 = (x0, v_x0, y0, v_y0, z0, v_z0)
    sol = odeint(dSdt, y0=S0, t=t, tfirst=True)
    
#     if dSdt == dSdt_J2:
#         x = sol.T[0]*np.sin(sol.T[2])*np.cos(sol.T[4])
#         y = sol.T[2]*np.sin(sol.T[2])*np.sin(sol.T[4])
#         z = sol.T[4]*np.cos(sol.T[2])
#         r2 = np.array([x,y,z]) #orbit of the follower satellite
    
#     else:
    r2 = np.array([sol.T[0], sol.T[2], sol.T[4]])
    
    return r2

def inertial_reference(dSdt, t):
    S0 = (x0, v_x0, y0, v_y0, z0, v_z0)
    sol = odeint(dSdt, y0=S0, t=t, tfirst=True)
    
    #defining x,y,z in inertial frame
    x = sol.T[0]
    y = sol.T[2]
    z = sol.T[4]
    
    x_hill = np.sqrt(x**2+y**2+z**2)*np.cos(np.arctan(y/x))*np.sin(np.arccos(z/np.sqrt(x**2+y**2+z**2)))
    y_hill = np.sqrt(x**2+y**2+z**2)*np.sin(np.arctan(y/x))*np.sin(np.arccos(z/np.sqrt(x**2+y**2+z**2)))
    z_hill = z
    
    # Defiining r1 and r2
    r1 = np.array([np.cos(n*t), np.sin(n*t), [0] * len(t)]) #motion of the leader satellite
    
    x_inertial = r1[0] + x_hill
    y_inertial = r1[1] + y_hill
    z_inertial = r1[2] + z_hill
    r2 = np.array([x_inertial, y_inertial, z_inertial]) #orbit of the follower satellite
    
    return r1, r2


# In[3]:


fig = None #Initialize fig, follower as global variables
follower = None
leader = None

def plot(dSdt, reference, t): # Creating a 3D plot
    global fig, follower, leader # Use the global variables
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    if reference == hill_reference:
        r2 = hill_reference(dSdt, t)
        leader = ax.plot(0, 0, 0, 'bo', markersize=7, label = 'Leader')
        follower = ax.scatter([], [], [], c='r', marker='.', linewidths=0.5, label = 'Follower') # Scatter plot for the followersatellite  
        ax.set_title('Orbit of follower satellite in the Hill frame')
        
    elif reference == inertial_reference:
        r2 = inertial_reference(dSdt, t)[1]
        earth = ax.plot(0, 0, 0, 'go', markersize=7, label = 'Earth')
        leader = ax.scatter([], [], [], c='b', marker='.', linewidths=0.5, label = 'Leader')
        follower = ax.scatter([], [], [], c='r', marker='.', linewidths=0.5, label = 'Follower') # Scatter plot for the follower satellite
        ax.set_title('Orbit of follower satellite in the Inertial frame')    
        
    ax.legend()
    ax.set_xlabel('$ê_{x}$')  
    ax.set_ylabel('$ê_{y}$')
    ax.set_zlabel('$ê_{z}$')
    
    # Calculate the range for each axis
    x_range = np.max(r2[0]) - np.min(r2[0])
    y_range = np.max(r2[1]) - np.min(r2[1])
    z_range = np.max(r2[2]) - np.min(r2[2])

    # Set a percentage margin (adjust as needed)
    margin_percent = 0.1

    # Calculate the margin for each axis
    x_margin = x_range * margin_percent
    y_margin = y_range * margin_percent
    z_margin = z_range * margin_percent

    # Set axis bounds for the 3D plot with dynamic margin
    ax.set_xlim(np.min(r2[0]) - x_margin, np.max(r2[0]) + x_margin)
    ax.set_ylim(np.min(r2[1]) - y_margin, np.max(r2[1]) + y_margin)
    ax.set_zlim(np.min(r2[2]) - z_margin, np.max(r2[2]) + z_margin)


# In[4]:


def anim_orbits(dSdt, reference, t): # to plot and create the animation
    
    plot(dSdt, reference, t)
    S0 = (x0, v_x0, y0, v_y0, z0, v_z0)
    sol = odeint(dSdt, y0=S0, t=t, tfirst=True)
    
    def update(frame, reference):
        if reference == hill_reference:
            r2 = hill_reference(dSdt, t)
            follower._offsets3d = (r2[0, :frame], r2[1, :frame], r2[2, :frame])
            return [follower]
    
        elif reference == inertial_reference:
            r1 = inertial_reference(dSdt, t)[0]
            r2 = inertial_reference(dSdt, t)[1]
            leader._offsets3d = (r1[0, :frame], r1[1, :frame], r1[2, :frame])
            follower._offsets3d = (r2[0, :frame], r2[1, :frame], r2[2, :frame])
            return [leader, follower]
    
    ani = FuncAnimation(fig, update, frames=len(sol), interval=len(sol), fargs=(reference,), blit=True)
    
    # display the animation in jupyter notebook
    return HTML(ani.to_jshtml())


# In[5]:


t = np.linspace(0, 50, 1000) #sets the number of laps taken


# In[6]:


n=1
x0 = 0.01
v_x0 = 0
y0 = 0
v_y0 = -2*n*x0
z0 = 0
v_z0 = 0


# In[7]:


anim_orbits(dSdt_linear, hill_reference, t)


# In[8]:


anim_orbits(dSdt_linear, inertial_reference, t)


# In[9]:


v_y0 = -3/2*n*x0


# In[10]:


anim_orbits(dSdt_linear, hill_reference, t)


# In[11]:


anim_orbits(dSdt_linear, inertial_reference, t)


# In[12]:


omega = n = 1
mu = 0.01215
r0 = 1
x0 = 0.5
v_x0 = 0
y0 = 0
v_y0 = -2*n*x0
z0 = 0
v_z0 = 0


# In[13]:


anim_orbits(dSdt_nonlinear, hill_reference, t)


# In[14]:


anim_orbits(dSdt_nonlinear, inertial_reference, t)


# In[15]:


v_y0 = -3/2*n*x0


# In[16]:


anim_orbits(dSdt_nonlinear, hill_reference, t)


# In[17]:


anim_orbits(dSdt_nonlinear, inertial_reference, t)

