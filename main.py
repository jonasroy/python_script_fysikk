import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

#Ekspriment
def read_from_file_to_list(file):
    liste = []
    f = open(file,"r")
    read_liste = f.readlines("") #Legg til hvordan fila ser ut
    f.close()
    for i in range(0,len(read_liste)):
        liste.append(read_liste[i].split())
    return liste



def test_plot(N):
    import random
    liste = []
    for i in range(0,N):
        x = (random.randint(15,45))/100
        liste.append(x)
    return liste

def antall_datapunkt(Lengde,N):
    liste = []
    steglengde = Lengde/N
    x = 0
    for i in range(0,N):
        liste.append(x)
        x += steglengde
    return liste





plt.style.use('bmh')  # Nicer looking plots

# Properties of the rolling object
r = 0.01                 # m      (radius)
rho = 7200             # kg/m^3 (density)
g = 9.81               # m/s^2  (gravitational acceleration)

m = (4/3)*np.pi*r**3*rho # kg     (mass)
c = 2/5
I0 = c*m*r**2            # kg m^2 (moment of inertia)
print(m)
# Properties of the frame
L = 1.4                   # m (length)
yi = [.61, .355, .187, 0.185,
     .302, .375, .197,.0]  # m (y-positions)

N = len(yi)               #   (# of mounts)
xi = np.linspace(0, L, N) # m (x-positions)

# Callable class for the track curve
get_y = interp.CubicSpline(xi, yi, bc_type="natural")

get_dydx = get_y.derivative()  # Callable derivative of track
get_d2ydx2 = get_dydx.derivative()  # Callable double derivative of track

resultat = test_plot(N) ###HAHAHHAHAHA

#Eksprimentele






#Teorietiske

def to_theta(plot_liste,N,Lengde):
    import math
    liste = []
    for i in range(0,len(plot_liste)-1):
        delta_s = math.sqrt((plot_liste[i]-plot_liste[i+1])**2 + (Lengde/N)**2)
        liste.append(math.atan(1/delta_s))
    return liste

resultat = test_plot(10)


def get_theta(x):
    """ Returns the angle of the track. """
    return -np.arctan(get_dydx(x))


def get_R(x):
    """ Returns the radius of the curvature. """
    return (1 + (get_dydx(x)) ** 2) ** 1.5 / get_d2ydx2(x)


def get_curvature(x):
    """ Returns the curvature (1/R). """
    return get_d2ydx2(x) / (1 + (get_dydx(x)) ** 2) ** 1.5


x = np.linspace(xi[0], xi[-1], 200)


# Create figure
#fig, axarr = plt.subplots(3, 1,figsize=(10,8), dpi=400)
#fig.subplots_adjust(hspace=0.02)

# Axes 1:
#axarr[0].plot(x, get_y(x), 'C0', label=r"$y(x)$")
#axarr[0].plot(antall_datapunkt(L,10),resultat, "r") #Sett inn det eksprimentele plottet.
#axarr[0].plot(xi, yi, 'C1o', label="Mounts")
#axarr[0].set_ylabel(r"$y(x)$, [m]", size='15')
#axarr[0].set_aspect('equal')

# Axes 2:
#axarr[1].plot(x, get_theta(x), 'C0')
#axarr[1].plot(antall_datapunkt(L,10), to_theta(resultat,10,L), "r") #Samme her
#axarr[1].set_ylabel(r"$\theta(x)$, [rad]", size='15')

# Axes 2:
#axarr[2].plot(x, get_curvature(x), 'C0')
#axarr[2].plot()
#axarr[2].set_ylabel(r"$\kappa(x)$, [1/m]", size='15')
#axarr[2].set_xlabel(r'$t$, [s]')

plt.show()


#Numerikken

#Put extra kode her

def get_vdot(theta):
    """ Returns the time derivative of the (total) velocity. """
    return g * np.sin(theta) / (1 + c)





def RHS(z):
    """ Evaluates the right hand side of the two coupled
    ODEs given in the text.

    Parameters:
        z : array-like, len(2), float. [v, x]
            The first element is the velocity, the second is the x-position.

    Returns:
        array-like, len(2), float. [a, v]
        The first element is the time derivative of the velocity (acceleration),
        the second is the time derivative of the position (velocity).
    """
    # z = [x, v]
    # w = [vx, a]
    w = np.zeros(2)
    theta = get_theta(z[0])
    w[0] = z[1] * 1 / np.sqrt(1 + np.tan(theta) ** 2) #ACCELERATION
    w[1] = get_vdot(theta)
    return w

#def acce(z):
 #   acc = z * 1 / np.sqrt(1 + np.tan(get_theta(z)) ** 2)
  #  return acc







def rk4step(f, y, h):
    """ Performs one step of the 4th order Runge-Kutta method.

    Parameters:
        f : Callable function with one input parameter.
            The right hand side of the ODE. Note that the
            RHS is in our case not a function of time.
        y : array-like, float. Current position.
        h : float. Time step.
    """
    s1 = f(y)
    s2 = f(y + h * s1 / 2.0)
    s3 = f(y + h * s2 / 2.0)
    s4 = f(y + h * s3)

    return y + h / 6.0 * (s1 + 2.0 * s2 + 2.0 * s3 + s4)

dt = 1e-3 # s
tstop = 5 # s. If the ball has not reached the end within 5 seconds, we stop.
x0 = 0.11# m. Initial x-position  ## 11 cm var første x0 ballen misten kontakt med banen.
##Det vi fant ut var at den realistiske verdien er på 9 cm.
v0 = 0    # m/s. Inital velocity

def get_K(v):
    """ Returns the kinetic energy. """
    return .5*m*(1 + c)*v**2

def get_U(h):
    """ Returns the potential energy. """
    return m*g*h

# Set initial values
x = [x0]        # x-position
v = [v0]        # velocity
h = get_y(x0)   # height
K = [get_K(v0)] # kinetic energy
U = [get_U(h)]  # potential energy

it = 0 # Iterator
itmax = tstop/dt # Maximum number og iterations
while x0 <= L and it < itmax:
    # Perform one step of the Runge-Kutta method
    [x0, v0] = rk4step(RHS, [x0, v0], dt)
    # Append different values to their arrays
    x = np.append(x, x0)
    v = np.append(v, v0)
    h = get_y(x0)
    K = np.append(K, get_K(v0))
    U = np.append(U, get_U(h))
    it += 1

print("Iterations: %i"%(it))
print("Final time: %.2f s"%(it*dt))
dE = (K[0] - K[-1] + U[0] - U[-1])/(K[0] + U[0])
print("Relative change in mechanical energy: %.2e"%(dE))

t = np.linspace(0, it*dt, it + 1)

def acce(v):
    count = 0
    array = [v[1]/dt]
    for i in range(1,len(v)):
        dv = v[i] - v[i-1]
        a = dv/dt
        array.append(a)
    return array


# Create figure
#fig, axarr = plt.subplots(3, 1, sharex=True, figsize=(10, 8), dpi=400)
#fig.subplots_adjust(hspace=0.02)




# Axes 1:
#axarr[0].plot(t, x)
#axarr[0].set_ylabel(r"$x(t)$, [m]")

# Axes 2:
#axarr[1].plot(t, v)
#axarr[1].set_ylabel(r"$v(t)$, [m/s]")

# Axes 2:
#axarr[2].plot(t, U, label="Potential energy, $U$")
#axarr[2].plot(t, K, label="Kinetic energy, $K$")
#axarr[2].plot(t, K + U, label="Mechanical energy, $E=U+K$")
#axarr[2].set_ylabel(r"Energy, [J]")
#axarr[2].set_xlabel(r'$t$, [s]')
#axarr[2].legend()

#plt.show()

theta = get_theta(x)
kappa = get_curvature(x)
N = g*np.cos(theta) + v**2*kappa
index = np.argmax(N < 0)
if index <= 0:
    print("The ball did not loose contact.")
else:
    print("The ball lost contact at x = %.2f m."%(x[index]))
    print("the heigth was", get_y(x0)+0.20)

mu = 0.25
f = mu*N
Fgx = m*g*np.sin(theta)
index = np.argmax(f < Fgx)
if index <= 1:
    print("The ball did not glide.")
else:
    print("The ball began to glide at x = %.2f m."%(x[index]))

def get_Normal(v,x):
    theta = get_theta(x)
    kappa = get_curvature(x)
    array = g * np.cos(theta) + v ** 2 * kappa
    liste = []
    for i in range(len(array)):
        if array[i] < 0:
            liste.append(0)
        else:
            liste.append(array[i])
    return liste



fig, axarr = plt.subplots(2, 1,figsize=(10,8), dpi=400)
fig.subplots_adjust(hspace=0.02)

axarr[0].plot(x, get_y(x),"b")
axarr[0].plot(xi, yi, 'C1o', label="Mounts")
axarr[0].set_ylabel(r"$N(t)$, [m/s^2]")

#axarr[1].plot(x, acce(v),"r")
#axarr[1].set_ylabel(r"$a(x)$, [m/s^2]")




plt.show()

fig, axarr = plt.subplots(2, 1,figsize=(10,8), dpi=400)
fig.subplots_adjust(hspace=0.02)

axarr[0].plot(x, get_Normal(v,x),"r")
axarr[0].set_ylabel(r"$N(x)$, [m/s^2]")

#axarr[1].plot(t, get_Normal(v,x),"r")
#axarr[1].set_ylabel(r"$N(t)$, [m/s^2]")

plt.show()

