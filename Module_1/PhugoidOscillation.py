import numpy
from matplotlib import pyplot

T = 100.0 # Length of time interval
dt = .02 # delta t (time step)
N = int(T/dt) + 1 # number of time steps
t = numpy.arange(0,T+dt,dt) #time grid

# Initial conditions

z0 = 100.0 # altitude
b0 = 10.0 # z' ; upward velocity resulting from gust
zt = 100.0 # trim altitude -- trim altitude is the steady state altitude where the aerodone travels horizontally
g = 9.81 # you know what this is

# Initial value for u (our solution for 2nd order ODE)
u = numpy.array([z0, b0])

# Array storing z at each time step
z = numpy.zeros(N)
z[0] = z0

for n in range(1,N):
    uprime = numpy.array([u[1], g - g*u[0]/zt])
    u = u + dt * uprime
    z[n] = u[0]

# Set the font family and size to use for Matplotlib figures.
pyplot.rcParams['font.family'] = 'serif'
pyplot.rcParams['font.size'] = 16

# Plot the solution of the elevation.
pyplot.figure(figsize=(9.0, 4.0))  # set the size of the figure
pyplot.title('Elevation of the phugoid over the time')  # set the title
pyplot.xlabel('Time [s]')  # set the x-axis label
pyplot.ylabel('Elevation [m]')  # set the y-axis label
pyplot.xlim(t[0], t[-1])  # set the x-axis limits
pyplot.ylim(40.0, 160.0)  # set the y-axis limits
pyplot.grid()  # set a background grid to improve readability
pyplot.plot(t, z, color='C0', linestyle='-', linewidth=2)


# The exact solution
z_exact = (b0 * (zt / g)**0.5 * numpy.sin((g / zt)**0.5 * t) +
           (z0 - zt) * numpy.cos((g / zt)**0.5 * t) + zt)

# Plot versus the exact values to examine error

pyplot.figure(figsize=(9.0, 4.0))  # set the size of the figure
pyplot.title('Elevation of the phugoid over the time')  # set the title
pyplot.xlabel('Time [s]')  # set the x-axis label
pyplot.ylabel('Elevation [m]')  # set the y-axis label
pyplot.xlim(t[0], t[-1])  # set the x-axis limits
pyplot.ylim(40.0, 160.0)  # set the y-axis limits
pyplot.grid()  # set a background grid to improve readability
pyplot.plot(t, z, label='Numerical',
            color='C0', linestyle='-', linewidth=2)
pyplot.plot(t, z_exact, label='Analytical',
            color='C1', linestyle='-', linewidth=2)
pyplot.legend()  # set the legend

# Convergence check

# dt>0
dt_values = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]

# Initializing a list that will contain our solutions for each dt
z_values = []

for dt in dt_values:
    N = int(T/dt) + 1 # number of time steps
    t = numpy.arange(0,T+dt,dt) #time grid
    # Initial conditions
    u = numpy.array([z0,b0])
    z = numpy.empty_like(t)
    z[0] = z0

    # If I had known this was coming, I'd have made the stuff above a function. Maybe later.
    for n in range(1,N):
        uprime = numpy.array([u[1], g-g*u[0]/zt])
        u = u +dt*uprime
        z[n] = u[0] #Store each z at the timestep of n
    z_values.append(z) # Adds the z values of the timestep to our list

    # Now we're making the error function

def L1_error(z,z_exact,dt):
    """
    Shamelessly taken from the notebook I'm studying from-- I can write my own headers later...

    Computes and returns the error
    (between the numerical and exact solutions)
    in the L1 norm.
    
    Parameters
    ----------
    z : numpy.ndarray
        The numerical solution as an array of floats.
    z_exact : numpy.ndarray
        The analytical solution as an array of floats.
    dt : float
        The time-step size.
        
    Returns
    -------
    error: float
        L1-norm of the error with respect to the exact solution.
        """
    error = dt*numpy.sum(numpy.abs(z-z_exact))
    return error

# Back to checking convergence, we're gonna make a list to store our errors

error_values = []
for z, dt in zip(z_values,dt_values):
    N = int(T/dt) + 1 # number of time steps
    t = numpy.arange(0,T+dt,dt) #time grid
    # exact solution
    z_exact = (b0 * (zt / g)**0.5 * numpy.sin((g / zt)**0.5 * t) + (z0 - zt) * numpy.cos((g / zt)**0.5 * t) + zt)
    # Add the error into the list
    error_values.append(L1_error(z, z_exact, dt))

# Plot the errors on a double log scale!

# Plot the error versus the time-step size.
pyplot.figure(figsize=(6.0, 6.0))
pyplot.title('L1-norm error vs. time-step size')  # set the title
pyplot.xlabel('$\Delta t$')  # set the x-axis label
pyplot.ylabel('Error')  # set the y-axis label
pyplot.grid()
pyplot.loglog(dt_values, error_values,
              color='C0', linestyle='--', marker='o')  # log-log plot
pyplot.axis('equal')  # make axes scale equally
pyplot.show()

# The challenge is to implement the euler method as a function. You are now obligated to do it.