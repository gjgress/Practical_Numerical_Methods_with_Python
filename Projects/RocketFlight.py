import numpy
from matplotlib import pyplot
import math

# Set the font family and size to use for Matplotlib figures.
pyplot.rcParams['font.family'] = 'serif'
pyplot.rcParams['font.size'] = 16

# initial conditions

ms = 50 # weight of the rocket shell (kg)
g = 9.81 # gravity (m/s^2)
roh = 1.091 # average air density (kg/m^3)
A = math.pi*(.5)**2 # maximum cross sectional area of rocket (r is .5m)
ve = 325 # exhaust speed (m/s)
CD = 0.15 #drag coefficient
mpo = 100 # mp at t = 0, mp is the weight of the rocket propellant
mprime = 20 # burn rate of the fuel, or the change in mp
h0 = 0 # initial height
v0 = 0 # initial velocity


mp = numpy.copy(mpo)

# mp = mpo - integral from 0 to t (mprime * dtau)

# euler's method with timestep dt = .1s, calculates altitude and velocity of rocket from launch until crash

# equation of motion for rocket: dh/dt = v
# (ms + mp) * dv/dt = -(ms + mp)*g + mprime*ve - 1/2 * roh*v*abs(v)*A*CD

# variables

dt = .1 # timestep in seconds
u = numpy.empty((1,2)) # This is our solutions vector, which will contain h and v
u[0] = numpy.array([h0,v0])

N = 0 # number of timesteps

def uprime(u, ms, mp, g, mprime, ve, roh, A, CD): 
    """
    Returns the derivative of our system for the phugoid system of equations

    Parameters
    ----------
    u: list / numpy.ndarray
       this is an array containing our solution from the previous timestep, in the form of a vector
    ms: float
        mass of the shell
    mp: float
        mass of the fuel
    g: float
        gravity
    mprime: float
        burn rate of fuel
    ve: float
        exhaust velocity
    roh: float
        density of air
    A:  float
        cross-sectional area of rocket
    CD: float
        coefficient of drag

    Returns
    -------
    uprime: numpy.ndarray
            This is the iterated derivative of u, so we can multiply by our timestep and use euler's method
    """
    h,v = u
    uprimevector = numpy.array([v, (-(ms + mp)*g + mprime*ve - 1/2 * roh*v*abs(v)*A*CD)/(ms+mp)])
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


while u[N,0]>0 or N<10:


    u = numpy.vstack([u, euler_step(u[N], uprime, dt, ms, mp, g, mprime, ve, roh, A, CD)])
    
    if N*dt<5:
        mp = mp - dt*mprime 
    if mp==0:
        mprime = 0
    
    N = N + 1
    

T = (N-1)*dt
T = math.floor(T*100)/100.0
t = numpy.arange(0,T+dt,dt)

# Plot the path of the rocket.
pyplot.figure(figsize=(T+2, numpy.amax(u[:,0])+10))
pyplot.title('Path of the rocket (flight time = {})'.format(T))
pyplot.xlabel('t')
pyplot.ylabel('h')
pyplot.grid()
pyplot.plot(t, u[:,0], color='C0', linestyle='-', linewidth=2)

print(T)
print(u[:,0].shape)

print("Maximum speed is:" + str(numpy.amax(u[:,1])))
print("The time value that the maximum speed occurs is:" + str(float(numpy.where(u[:,1]==numpy.amax(u[:,1]))[0])*dt))
print("The altitude at this time is:" + str(u[numpy.where(u[:,1]==numpy.amax(u[:,1])),0]))
print("Maximum height is: " + str(numpy.amax(u[:,0])))
print("The maximum height occurs at t = " + str(float(numpy.where(u[:,0]==numpy.amax(u[:,0]))[0])*dt))
print("The rocket hits the ground at t = " + str(T))
print("The velocity of the rocket when it hits the ground is: " + str(u[-1,1]))