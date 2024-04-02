#!/usr/bin/env python
# coding: utf-8

# ## 2D visualization ($x$ and $y$ plane)

# In[1]:


#importing necessary libraries for calculations and simulations

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D


# In[2]:


def dSdt(t, S):
    x, v_x, y, v_y, z, v_z = S
    dSdt = [v_x, 
           2 * n * v_y + n**2 * x - mu*(r0+x)/((r0+x)**2+y**2+z**2)**(3/2) + mu/r0**2,
            v_y,
            -2 * n * v_x + n**2 * y - mu*y/((r0+x)**2+y**2+z**2)**(3/2),
            v_z,
            -mu*z/((r0+x)**2+y**2+z**2)**(3/2)]
    return dSdt


# In[3]:


def anim_orbits_hill(t, save_path=None): # to plot and create the animation
    
    S_0 = (x_0, v_x0, y_0, v_y0, z_0, v_z0)
    sol = odeint(dSdt, y0=S_0, t=t, tfirst=True)
    
    r2 = np.array([sol.T[0], sol.T[2]]) #orbit of the follower satellite
    
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    line2, = ax.plot([], [], 'r', ls ="--", linewidth=1.5)
    ball, = ax.plot(0, 0, 'bo', markersize=5, label = 'Leader') 
    ball2, = ax.plot([], [], 'ro', markersize=3, label = 'Follower')  
    ax.legend()
    ax.set_xlabel('$ê_{x}$')  
    ax.set_ylabel('$ê_{y}$')
    ax.set_title('Orbit of follower satellite in the Hill frame')
    
    ax.set_xlim(-.3, .3)
    ax.set_ylim(-.3, .3)
    
    def update(frame):
        line2.set_data(r2[0, :frame], r2[1, :frame])
        
        ball2.set_data(r2[0, frame], r2[1, frame]) # update the position of the balls for each frame
        
        ax.relim()
        ax.autoscale_view()
                
        return line2, ball, ball2

    ani = FuncAnimation(fig, update, frames=len(n*t), interval=50, blit=True)
    
    if save_path:
        # Save animation as a GIF file
        ani.save(save_path, writer='imagemagick')
    
    return HTML(ani.to_jshtml()) # display the animation

def anim_orbits_inertial(t, save_path=None): # to plot and create the animation
    
    S_0 = (x_0, v_x0, y_0, v_y0, z_0, v_z0)
    sol = odeint(dSdt, y0=S_0, t=t, tfirst=True)
    
    # Defiining r1 and r2
    r1 = np.array([np.cos(n*t),np.sin(n*t)]) #motion of the leader satellite
    x = r1[0] + np.multiply(sol.T[0], np.cos(n*t)) + np.multiply((sol.T[2]), np.cos(n*t+np.pi/2))
    y = r1[1] + np.multiply(sol.T[0], np.sin(n*t)) + np.multiply((sol.T[2]), np.sin(n*t+np.pi/2))
    r2 = np.array([x, y]) #orbit of the follower satellite
    
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    line1, = ax.plot([], [], 'b', linewidth=1)
    line2, = ax.plot([], [], 'r', ls ="--", linewidth=1.5)
    ball, = ax.plot(0, 0, 'go', markersize=5, label = 'Earth')
    ball1, = ax.plot([], [], 'bo', markersize=5, label = 'Leader') 
    ball2, = ax.plot([], [], 'ro', markersize=3, label = 'Follower')  
    ax.set_xlabel('$ê_{x}$')  
    ax.set_ylabel('$ê_{y}$')
    ax.legend()
    ax.set_title('Orbit of follower satellite in the Inertial frame')    
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
        
    def update(frame):
        line1.set_data(r1[0, :frame], r1[1, :frame])
        line2.set_data(r2[0, :frame], r2[1, :frame])
        
        # Update the position of the 'balls' for the most recent frame
        ball1.set_data(r1[0, frame], r1[1, frame])
        ball2.set_data(r2[0, frame], r2[1, frame])
        
        return line1, line2, ball1, ball2

    ani = FuncAnimation(fig, update, frames=len(n*t), interval=50, blit=True)
    
    if save_path:
        # Save animation as a GIF file
        ani.save(save_path, writer='imagemagick')
    
    # display the animation in jupyter notebook
    return HTML(ani.to_jshtml())


# ### First mode

# In[10]:


n = 1
mu = 1
r0 = 1
x_0 = 0.01 #must be small for the approximation be true
v_x0 = 0
y_0 = 0
v_y0 = -3/2*n*x_0 #initial conditions for the first mode
z_0 = 0
v_z0 = 0

t = np.linspace(0, 500, 1000) #sets the number of laps taken


# In[5]:


anim_orbits_hill(t, save_path='nonlinear_hill_animation_firstmode.gif')


# In[6]:


anim_orbits_inertial(t, save_path='nonlinear_eci_animation_firstmode.gif')


# ### Second mode

# In[11]:


x_0 = 0.05
v_y0 = -2*n*x_0 #second mode initial conditions


# In[12]:


anim_orbits_hill(t, save_path='nonlinear_hill_animation_secondmode.gif')


# In[13]:


anim_orbits_inertial(t, save_path='nonlinear_eci_animation_secondmode.gif')


# ### Third mode

# In[10]:


x_0 = 0.0 #must be small for the approximation be true
v_x0 = 2/n*y_0 #third mode initial conditions
y_0 = 0.1
v_y0 = 0 #initial conditions for the first mode
z_0 = 0
v_z0 = 0


# In[11]:


anim_orbits_hill(t, save_path='nonlinear_hill_animation_secondmode.gif')


# In[12]:


anim_orbits_inertial(t, save_path='nonlinear_eci_animation_secondmode.gif')

