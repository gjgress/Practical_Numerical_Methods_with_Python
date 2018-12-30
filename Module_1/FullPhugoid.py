import math
import numpy
from matplotlib import pyplot

# Set the font family and size to use for Matplotlib figures.
pyplot.rcParams['font.family'] = 'serif'
pyplot.rcParams['font.size'] = 16

# Parameter and initial values

g = 9.81 # you know what this is
vt = 4.9 # trim velocity (in m/s) -- the velocity the aerodone remains horizontal
CD = 1.0/5 # drag coefficient
CL = 1.0 # lift coefficient (normal valus for CD/CL are ~40-60 for a sailplane)

v0 = 6.5 # initial velocity - we will start at trim velocity
theta0 = -0.1 # trajectory angle
x0 = 0.0 # horizontal initial position
y0 = 10.0 # vertical initial position (altitude)

def uprime(u, CL, CD, g, vt):
    """
    Returns the derivative of our system for the phugoid system of equations

    Parameters
    ----------
    u: list / numpy.ndarray
       this is an array containing our solution from the previous timestep, in the form of a vector
    CL: float
        the coefficient for lift
    CD: float
        the coefficient for drag
    g: float
        gravity, duh
    vt: float
        trim velocity

    Returns
    -------
    uprime: numpy.ndarray
            This is the iterated derivative of u, so we can multiply by our timestep and use euler's method
    """
    v, theta, x, y = u
    uprimevector = numpy.array([-g * math.sin(theta) - CD/CL * g / vt**2 * v**2, -g * math.cos(theta)/v+g/vt**2 * v, v*math.cos(theta),v*math.sin(theta)])

    return uprimevector

def euler_step(u, f, dt, *args):
    """
    I wrote out the first one myself, but not this one

    Returns the solution at the next time step using Euler's method.
    
    Parameters
    ----------
    u : numpy.ndarray
        Solution at the previous time step
        as a 1D array of floats.
    f : function
        Function to compute the right-hand side of the system.
    dt : float
        Time-step size.
    args : tuple, optional
        Positional arguments to pass to the function f.
    
    Returns
    -------
    u_new : numpy.ndarray
        The solution at the next time step
        as a 1D array of floats.
    """
    u_new = u + dt*f(u,*args)
    return u_new



dt = 0.01  # time-step size

# Create array to store the solution at each time step.
u = numpy.empty((1, 4))
# Set the initial conditions.
u[0] = numpy.array([v0, theta0, x0, y0])

# Time integration with Euler's method.

# timestep counter
N = 0

while u[N,3]>0:
    u = numpy.vstack([u, euler_step(u[N], uprime, dt, CL, CD, g, vt)])
    N = N+1
    
T = N*dt

x = u[:,2]
y = u[:,3]

# Plot the path of the glider.
pyplot.figure(figsize=(9.0, 4.0))
pyplot.title('Path of the glider (flight distance = {})'.format(x[-1]))
pyplot.xlabel('x')
pyplot.ylabel('y')
pyplot.grid()
pyplot.plot(x, y, color='C0', linestyle='-', linewidth=2)

pyplot.show()



# # Set the list of time-step sizes.
# dt_values = [0.1, 0.05, 0.01, 0.005, 0.001]

# # Create an empty list that will contain the solution of each grid.
# u_values = []

# for dt in dt_values:
#     N = int(T / dt) + 1  # number of time-steps
#     # Create array to store the solution at each time step.
#     u = numpy.empty((N, 4))
#     # Set the initial conditions.
#     u[0] = numpy.array([v0, theta0, x0, y0])
#     # Temporal integration using Euler's method.
#     for n in range(N - 1):
#         u[n + 1] = euler_step(u[n], uprime, dt, CL, CD, g, vt)
#     # Store the solution for the present time-step size
#     u_values.append(u)

# def l1_diff(u_coarse, u_fine, dt):
#     """
#     Returns the difference in the L1-norm between the solution on
#     a coarse grid and the solution on a fine grid.
    
#     Parameters
#     ----------
#     u_coarse : numpy.ndarray
#         Solution on the coarse grid as an array of floats.
#     u_fine : numpy.ndarray
#         Solution on the fine grid as an array of floats.
#     dt : float
#         Time-step size.
    
#     Returns
#     -------
#     diff : float
#         The difference between the two solutions in the L1-norm
#         scaled by the time-step size.
#     """
#     N_coarse = len(u_coarse)
#     N_fine = len(u_fine)
#     ratio = math.ceil(N_fine / N_coarse)
#     diff = dt * numpy.sum(numpy.abs(u_coarse - u_fine[::ratio]))
#     return diff



# # Create an empty list to store the difference in the solution
# # between two consecutive grids.
# diff_values = []

# for i, dt in enumerate(dt_values[:-1]):
#     diff = l1_diff(u_values[i][:, 2], u_values[-1][:, 2], dt)
#     diff_values.append(diff)



# # Plot the difference versus the time-step size.
# pyplot.figure(figsize=(6.0, 6.0))
# pyplot.title('L1-norm difference vs. time-step size')  # set the title
# pyplot.xlabel('$\Delta t$')  # set the x-axis label
# pyplot.ylabel('Difference')  # set the y-axis label
# pyplot.grid()
# pyplot.loglog(dt_values[:-1], diff_values,
#               color='C0', linestyle='--', marker='o')  # log-log plot
# pyplot.axis('equal');  # make axes scale equally

# r = 2  # refinement ratio for the time-step size
# h = 0.001  # base grid size

# dt_values2 = [h, r * h, r**2 * h]
# u_values2 = []

# for dt in dt_values2:
#     N = int(T / dt) + 1  # number of time steps
#     # Create array to store the solution at each time step.
#     u = numpy.empty((N, 4))
#     # Set initial conditions.
#     u[0] = numpy.array([v0, theta0, x0, y0])
#     # Time integration using Euler's method.
#     for n in range(N - 1):
#         u[n + 1] = euler_step(u[n], uprime, dt, CL, CD, g, vt)
#     # Store the solution.
#     u_values2.append(u)

# # Calculate f2 - f1.
# f2_f1 = l1_diff(u_values2[1][:, 2], u_values2[0][:, 2], dt_values2[1])
# # Calculate f3 - f2.
# f3_f2 = l1_diff(u_values2[2][:, 2], u_values2[1][:, 2], dt_values2[2])
# # Calculate the observed order of convergence.
# p = math.log(f3_f2 / f2_f1) / math.log(r)
# print('Observed order of convergence: p = {:.3f}'.format(p))