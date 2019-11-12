import numpy as np

class PDE2_scheme(object):
    def __init__(self, dt, nx, init, left_b, right_b, x0, x1, btype, nt, left_r = None, right_r = None):
        # size of time step
        self.dt = dt
        # type of condition
        self.btype = btype
        # set left boundary function
        self.set_lb(left_b)
        # set right boundary function
        self.set_rb(right_b)
        # set initial data
        self.set_init(init)
        # array of x, dx
        self.x, self.dx = np.linspace(x0, x1, nx, retstep=True, dtype=float)
        # append x1+dx
        self.x = np.insert(self.x, nx, x1 + self.dx)
        # prepend x0-dx
        self.x = np.insert(self.x, 0, x0 - self.dx)
        # number of discretized x positions
        self.nx = nx + 2
        # u[x] at t=n*dt till t=(n+nt-1)*dt
        self.u = np.zeros((self.nx, nt), dtype=float);
        # set the x domain for stepping scheme
        self.set_xDomain()
    # set function for left boundary
    def set_lb(self, left_b, left_r = None):
        self.left_b = left_b
        self.left_r = left_r
    # set function for right boundary
    def set_rb(self, right_b, right_r = None):
        self.right_b = right_b
        self.right_r = right_r
    # set function for initial data
    def set_init(self, init):
        self.init = init
    # set the x domain for stepping scheme [self.lb,self.rb) in u
    def set_xDomain(self):
        # left
        if self.btype[0] == 0:
            self.lb = 2
        else:
            self.lb = 1
        # right
        if self.btype[1] == 0:
            self.rb = self.nx-2
        else:
            self.rb = self.nx-1
    # evaluate boundary at time t=it*dt
    def eval_boundary(self, t, it):
        # left
        if self.btype[0] == 2:
            #robin
            r = self.left_r(t)
            if r == 0:
                # degenerate to dirichlet
                set_boundary(self.u[:,it], 0, self.btype[0], self.left_b(t),  self.dx)
            else:
                set_boundary(self.u[:,it], 0, self.btype[0], self.left_b(t),  self.dx, self.left_r(t))
        else:
            set_boundary(self.u[:,it], 0, self.btype[0], self.left_b(t),  self.dx)
        # right
        if self.btype[1] == 2:
            #robin
            r = self.right_r(t)
            if r == 0:
                # degenerate to dirichlet
                set_boundary(self.u[:,it], 1, self.btype[1], self.right_b(t), self.dx)
            else:
                set_boundary(self.u[:,it], 1, self.btype[1], self.right_b(t), self.dx, self.right_r(t))
        else:
            set_boundary(self.u[:,it], 1, self.btype[1], self.right_b(t), self.dx)

    def reset_dx(self, dx):
        self.dx = dx
        self.set_coef()

    def reset_dt(self, dt):
        self.dx = dx
        self.set_coef()

    def set_coef(self):
        pass 

class forward_Euler(PDE2_scheme):
    def __init__(self, dt, nx, init, left_b, right_b, x0, x1, btype, D, left_r = None, right_r = None):
        PDE2_scheme.__init__(self, dt, nx, init, left_b, right_b, x0, x1, btype, 2, left_r, right_r)
        # set the diffusion coefficient 
        self.set_D(D)
        # initialize
        self.initialize()
    # initialization
    def initialize(self):
        # set the coefficient in the scheme D * dt / dx^2
        self.set_coef()
        # init stepping direction in u
        self.step_dir = np.array([0, 1])
        # set initial data
        self.u[self.lb:self.rb,0] = self.init(self.x[self.lb:self.rb])
        # print growth factor
        print(f'growth factor = {2*self.coef:.3f}')
    # step forward in time
    def step(self, t):
        # set alias
        d = self.step_dir
        u = self.u
        # set boundaries for u for t=n*dt 
        self.eval_boundary(t, d[0])
        # evolve u to t=(n+1)*dt for x=j*dx, for j in [self.lb,self.rb)
        for i in range(self.lb, self.rb):
            u[i,d[1]] = self.coef * (u[i+1,d[0]] - 2*u[i,d[0]] + u[i-1,d[0]]) + u[i,d[0]]
        # alternate stepping direction for next step
        self.alternate();
        return u[1:self.nx-1,d[0]]
    # alternate stepping direction
    def alternate(self):
        self.step_dir = 1 - self.step_dir

    def set_D(self, D):
        self.D = D
        self.set_coef()

    def set_coef(self):
        self.coef = self.D * self.dt/np.power(self.dx,2)

class leap_frog(PDE2_scheme):
    def __init__(self, dt, nx, init, left_b, right_b, x0, x1, btype, C2, vel, left_r = None, right_r = None):
        PDE2_scheme.__init__(self, dt, nx, init, left_b, right_b, x0, x1, btype, 3, left_r, right_r)
        # set the square of wave speed
        self.set_C2(C2) 
        # set velicty function
        self.set_vel(vel)
        # initialize
        self.initialize()
    # initialize with inital data
    def initialize(self):
        # set the coefficient in the scheme D * dt / dx^2
        self.set_coef()
        # init leap sequence in u
        self.leap = np.array([0, 1, 2])
        # set initial data
        self.u[self.lb:self.rb,0] = self.init(self.x[self.lb:self.rb])
        # set boundary values for t=0
        self.eval_boundary(0.0, 0)
        # alias
        u = self.u
        # set u at t=dt with vel data
        for i in range(self.lb, self.rb):
            u[i,1] = u[i,0] + self.dt * self.vel(self.x[i]) + 0.5*self.coef*(u[i+1,0] - 2*u[i,0] + u[i-1,0])
        # print growth factor
        print(f'growth factor = {self.coef:.3f}')
    # step forward in time
    def step(self, t):
        # alias
        d = self.leap
        u = self.u
        # set boundaries for t=(n+1)*dt
        self.eval_boundary(t, d[1])
        # evolve u to (n+2)*dt for j*dx, j=1,nx-1
        for i in range(self.lb, self.rb):
            u[i,d[2]] = self.coef * (u[i+1,d[1]] - 2*u[i,d[1]] + u[i-1,d[1]]) + 2*u[i,d[1]] - u[i,d[0]]
        # permute leap sequence for next time step
        self.permute_leap()
        return u[1:self.nx-1,d[1]]
    # permute leap sequence 
    def permute_leap(self):
        self.leap = np.roll(self.leap,-1)

    # set function for initial data
    def set_vel(self, vel):
        self.vel = vel 

    def set_C2(self, C2):
        self.C2 = C2
        self.set_coef()

    def set_coef(self):
        self.coef = self.C2 * np.power(self.dt/self.dx,2)


def set_boundary(u, lr, btype, f, dx = 0.0, r = None):
    # set left index
    if lr == 0:
        _ex = 0  # 1 dx outside the boundary
        _on = 1  # on the boundary
        _in = 2  # 1 dx inside the boundary
    # set right index 
    else:
        _ex = -1
        _on = -2
        _in = -3
    # set boundary
    if btype == 0:
        u[_on] = set_dirichlet(f)
    if btype == 1:
        u[_ex] = set_neumann(u[_in], f, dx)
    if btype == 2:
        u[_ex] = set_robin(u[_on], u[_in], f, r, dx)

def set_robin(u_on, u_in, f, r, dx):
    u_ex = (f - u_on)*(2*dx)/r + u_in
    return u_ex
     
def set_neumann(u_in, f, dx):
    u_ex = f*(2*dx) + u_in
    return u_ex

def set_dirichlet(f):
    return f
