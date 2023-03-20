import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


class TFMat:
    def __init__(self, d, theta, r, alpha, deg=1):
        if deg:
            theta *= np.pi / 180
            alpha *= np.pi / 180
        self.d = d
        self.theta = theta
        self.r = r
        self.alpha = alpha
        self.T = np.zeros([4, 4])
        self.dhtrans()

    def dhtrans(self):
        cost = np.cos(self.theta)
        sint = np.sin(self.theta)
        cosa = np.cos(self.alpha)
        sina = np.sin(self.alpha)

        self.T = np.array([[cost, -sint * cosa, sint * sina, self.r * cost],
                           [sint, cost * cosa, -cost * sina, self.r * sint],
                           [0, sina, cosa, self.d],
                           [0, 0, 0, 1]])
        return self.T


def plot_line(ax, vec_start, vec_end, color='b'):
    ax.plot([vec_start[0], vec_end[0]], [vec_start[1], vec_end[1]], zs=[vec_start[2], vec_end[2]], color=color)


def plot_tf(ax, tfmat):
    vo = np.array([0, 0, 0, 1])
    vi = np.array([1, 0, 0, 1])
    vj = np.array([0, 1, 0, 1])
    vk = np.array([0, 0, 1, 1])
    uvecs = [vo, vi, vj, vk]
    uvecs_t = []
    for i in range(len(uvecs)):
        uvecs_t.append(tfmat @ uvecs[i])

    colors = ["r", "g", "b"]

    for ii in range(len(uvecs) - 1):
        plot_line(ax, uvecs_t[0], uvecs_t[ii + 1], colors[ii])
    return uvecs_t[0]


def generate_sliders(num, plotfig):
    d_sliders = []
    theta_sliders = []
    r_sliders = []
    alpha_sliders = []

    total_sliders = num * 5 + 4
    step = 1 / total_sliders
    start_coord = step * 2
    plt.figure()

    for iii in range(num):
        d_sliders.append(
            Slider(plt.axes([0.25, start_coord + step * (iii * 5), 0.65, 0.03]), 'D' + str(iii + 1), -10, 10,
                   valinit=iii, valstep=0.5))
        theta_sliders.append(
            Slider(plt.axes([0.25, start_coord + step * (iii * 5 + 1), 0.65, 0.03]), 'theta' + str(iii + 1), -180, 180,
                   valinit=0, valstep=0.5))
        r_sliders.append(
            Slider(plt.axes([0.25, start_coord + step * (iii * 5 + 2), 0.65, 0.03]), 'r' + str(iii + 1), -10, 10,
                   valinit=0, valstep=0.5))
        alpha_sliders.append(
            Slider(plt.axes([0.25, start_coord + step * (iii * 5 + 3), 0.65, 0.03]), 'alpha' + str(iii + 1), -180, 180,
                   valinit=0, valstep=0.5))

    sliders = [d_sliders, theta_sliders, r_sliders, alpha_sliders]
    for iii in range(num):
        d_sliders[iii].on_changed(lambda val: update(val, sliders, plotfig))
        theta_sliders[iii].on_changed(lambda val: update(val, sliders, plotfig))
        r_sliders[iii].on_changed(lambda val: update(val, sliders, plotfig))
        alpha_sliders[iii].on_changed(lambda val: update(val, sliders, plotfig))

    return sliders


def generate_tfchain(ds, thetas, rs, alphas):
    tfmats = []
    initial_tf = np.identity(4)
    final_tf = initial_tf
    tfchain = [initial_tf]
    for i in range(numOuter):
        tfmats.append(TFMat(ds[i], thetas[i], rs[i], alphas[i], 1))

    for ii in range(len(tfmats)):
        final_tf = final_tf @ tfmats[ii].T
        tfchain.append(final_tf)
    return tfchain


def plot_tfchain(tfchain, ax):
    points = []
    for tf in tfchain:
        points.append(plot_tf(ax, tf))

    for i in range(len(points) - 1):
        plot_line(ax, points[i], points[i + 1], "k")
    ax.set_aspect('equal', 'box')


def update(_, sliders, fig):
    ax = fig.axes[0]
    d_sliders = sliders[0]
    theta_sliders = sliders[1]
    r_sliders = sliders[2]
    alpha_sliders = sliders[3]
    ax.clear()
    ds = []
    thetas = []
    rs = []
    alphas = []
    for i in range(len(sliders[0])):
        ds.append(d_sliders[i].val)
        thetas.append(theta_sliders[i].val)
        rs.append(r_sliders[i].val)
        alphas.append(alpha_sliders[i].val)
    tfchain = generate_tfchain(ds, thetas, rs, alphas)
    plot_tfchain(tfchain, ax)
    fig.canvas.draw_idle()


if __name__ == "__main__":
    numOuter = 5
    figOuter = plt.figure()
    axOuter = figOuter.add_subplot(111, projection='3d')
    slidersOuter = generate_sliders(numOuter, figOuter)
    dsOuter = range(numOuter)
    thetasOuter = [0] * numOuter
    rsOuter = [0] * numOuter
    alphasOuter = [0] * numOuter
    tfchainOuter = generate_tfchain(dsOuter, thetasOuter, rsOuter, alphasOuter)
    plot_tfchain(tfchainOuter, axOuter)
    plt.show()
