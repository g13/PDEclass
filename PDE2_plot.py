import matplotlib.pyplot as plt
from matplotlib import animation, rc 
from IPython.display import HTML
rc('animation', html='html5')

def plot(u, dt, dx, title):
    fig, ax = plt.subplots()
    im = ax.imshow(u)
    plt.imshow(u.transpose(), origin='lower', aspect='auto')
    plt.colorbar()
    ax.set_ylabel(f'dt = {dt:.3f}')
    ax.set_xlabel(f'dx = {dx:.3f}')
    ax.set_title(title)
    plt.show()

def capture(it, line, dt, u, x):
    i = int(round(it*1000/float(dt)))
    line.set_data(x, u[:, i])
    return (line,)

def movie(u, nt, x, title):
    fig, ax = plt.subplots()
    ax.set_xlim((x[0], x[-1]))
    ax.set_ylim((-2, 2))
    line, = ax.plot([], [], lw=2)
    for_it = 10 # constant dont change it
    m = animation.FuncAnimation(fig, capture, frames = nt, interval = 20, fargs = (line, for_it, u, x), blit = True)
    HTML(m.to_html5_video())
    m.save(title+".mp4")
