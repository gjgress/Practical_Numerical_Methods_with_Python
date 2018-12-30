# This module is essentially a repeat of the first steps of the CFD notebooks I've been working through
# However, this module introduces animations-- and so I will use this opportunity as an exercise to learn animations with python, 
# as well as utilizing array operations

import numpy
from matplotlib import pyplot
import math
from matplotlib import animation
from IPython.core.display import display, HTML

nx = 41
dx = 2 / (nx - 1)
nu = 0.3   #the value of viscosity
sigma = .2 #sigma is our CFL constant (Courant number)
dt = sigma * dx**2 / nu #dt is defined using sigma ... more later!
nt = int(.25/dt)    #the number of timesteps we want to calculate

# Get the grid point coordinates.
x = numpy.linspace(0.0, 2, num=nx)

u0 = numpy.ones(nx)
u0[int(.5/dx):int(1/dx+1)] = 2 # creating our hat functon again

def diffusion(u0, sigma=0.5, nt=20):
    """
    Computes the numerical solution of the 1D diffusion equation
    over the time steps.
    
    Parameters
    ----------
    u0 : numpy.ndarray
        The initial conditions as a 1D array of floats.
    sigma : float, optional
        The value of nu * dt / dx^2;
        default: 0.5.
    nt : integer, optional
        The number of time steps to compute;
        default: 20.
    
    Returns
    -------
    u_hist : list of numpy.ndarray objects
        The history of the numerical solution.
    """

    u = numpy.copy(u0) # our temp array for iteration
    u_hist = [u0.copy()]
    for i in range(nt):
        u[1:-1] = u[1:-1] + nu*dt/dx**2 * (u[2:] - 2*u[1:-1] + u[:-2])
        u_hist.append(u.copy())
    return u_hist

u_hist = diffusion(u0,sigma=sigma,nt=nt)

fig = pyplot.figure(figsize=(6.0,4.0))
pyplot.xlabel('x')
pyplot.ylabel('u')
pyplot.grid()
line = pyplot.plot(x,u0, color = 'C0', linestyle='-', linewidth=2)[0]
pyplot.xlim(0.0,2)
pyplot.ylim(0.5,2.5)
fig.tight_layout()

def update_plot(n,u_hist):
    """
    Update the line y-data of the Matplotlib figure.
    
    Parameters
    ----------
    n : integer
        The time-step index.
    u_hist : list of numpy.ndarray objects
        The history of the numerical solution.
    """
    fig.suptitle('Time step {:0>2}'.format(n))
    line.set_ydata(u_hist[n])



# Create an animation.
anim = animation.FuncAnimation(fig, update_plot, frames=nt, fargs=(u_hist,), interval=100)

anim.save('1D_Diffusion.mp4', fps=30, extra_args=['-vcodec', 'libx264'])


# # Display the video.
# HTML(anim.to_html5_video())


