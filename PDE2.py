from PDE2_scheme import *
from PDE2_plot import *

# initial velcity
def vel(x):
    #return 1 + np.sin(x)
    #return 1 + np.sin(x)
    return np.zeros(x.shape)
# initial data
def init(x):
    #return 1 + np.sin(x)
    return np.exp(-np.power(x,2)/0.1)
    #return np.zeros(x.shape)

#  u + left(right)_r*ux = left(right)_b
def left_r(t):
    #return np.power(t, 1)
    return 0

def right_r(t):
    #return np.power(t, 1)
    return 0

def left_b(t):
    #return np.power(t, 2)
    return 0

def right_b(t):
    #return np.power(t, 2)
    return 0

def run_heat(dt, nx, x0, x1, D, btype, nt, init, left_b, right_b, left_r, right_r, title = 'heat', mov = False):
    heat = forward_Euler(dt, nx, init, left_b, right_b, x0, x1, btype, D, left_r, right_r)
    u = np.zeros((nx, nt), dtype=float)
    for i in range(nt):
        u[:,i] = heat.step(i*dt)
    plot(u, dt, heat.dx, title)
    if mov:
        movie(u, np.arange(nt)*dt, heat.x[1:heat.nx-1], title)

def run_wave(dt, nx, x0, x1, C2, btype, nt, init, vel, left_b, right_b, left_r, right_r, title = 'wave', mov = False):
    wave = leap_frog(dt, nx, init, left_b, right_b, x0, x1, btype, C2, vel, left_r, right_r)
    u = np.zeros((nx, nt), dtype=float)
    for i in range(nt):
        u[:,i] = wave.step(i*dt)
    plot(u, dt, wave.dx, title)
    if mov:
        m = movie(u, np.arange(nt)*dt, wave.x[1:wave.nx-1], title)

if __name__ == '__main__':
    x0 = -np.pi
    x1 = np.pi
    dt = 0.01
    C2 = 10
    D = 0.2

    btype = np.array([0, 0]) # 0 for dirichlet, 1 for neumann, 2 for robin
    title = 'wave flip'
    nt = 200
    nx = 100
    run_wave(dt, nx, x0, x1, C2, btype, nt, init, vel, left_b, right_b, left_r, right_r, title)

    init_h = lambda x: np.exp(-np.power(x,2)/20.0)
    nt = 2000
    nx = 50
    title = 'heat dissipate'
    run_heat(dt, nx, x0, x1, D, btype, nt, init_h, left_b, right_b, left_r, right_r, title)

    btype = np.array([1, 1]) # 0 for dirichlet, 1 for neumann, 2 for robin
    left_bn = lambda x: 0
    right_bn = lambda x: 0
    nt = 200
    nx = 100
    title = 'wave reflect'
    run_wave(dt, nx, x0, x1, C2, btype, nt, init, vel, left_bn, right_bn, left_r, right_r, title)

    left_bn = lambda x: 0
    right_bn = lambda x: 0
    init_h = lambda x: np.exp(-np.power(x,2)/20.0)
    nt = 2000
    nx = 50
    title = 'heat insolation'
    run_heat(dt, nx, x0, x1, D, btype, nt, init_h, left_bn, right_bn, left_r, right_r, title)

    btype = np.array([1, 1]) # 0 for dirichlet, 1 for neumann, 2 for robin
    left_bn = lambda x: 0
    right_bn = lambda x: 0
    nt = 200
    nx = 100
    dt = dt*2.0
    C2 = 1.01/np.power(dt/((x1-x0)/(nx-1)),2)
    title = 'unstable wave'
    run_wave(dt, nx, x0, x1, C2, btype, nt, init, vel, left_bn, right_bn, left_r, right_r, title)
