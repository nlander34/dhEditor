import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

class tfmat:
    def __init__(self,d, theta, r, alpha, deg=1):
        if deg:
            theta *= np.pi/180
            alpha *= np.pi/180
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


def plotLine(ax, VecStart, VecEnd, color='b'):
    ax.plot([VecStart[0], VecEnd[0]], [VecStart[1], VecEnd[1]], zs=[VecStart[2], VecEnd[2]], color=color)


def plotTF(ax, T):
    vo = np.array([0, 0, 0, 1])
    vi = np.array([1, 0, 0, 1])
    vj = np.array([0, 1, 0, 1])
    vk = np.array([0, 0, 1, 1])
    uvecs = [vo, vi, vj, vk]
    uvecsT = []
    for i in range(len(uvecs)):
        uvecsT.append(T @ uvecs[i])

    colors = ["r", "g", "b"]

    for ii in range(len(uvecs)-1):
        plotLine(ax, uvecsT[0], uvecsT[ii+1], colors[ii])
    return uvecsT[0]


def generateSliders(num, plotfig):
    dSliders = []
    thetaSliders = []
    rSliders = []
    alphaSliders = []

    totalSliders = num * 5 + 4
    step = 1 / totalSliders
    startCoord = step * 2
    sliderFig = plt.figure()
    # sliderAx = sliderFig.add_subplot(111)

    for iii in range(num):
        dSliders.append(Slider(plt.axes([0.25, startCoord + step * (iii * 5), 0.65, 0.03]), 'D' + str(iii + 1), -10, 10,
                               valinit=iii, valstep=0.5))
        thetaSliders.append(
            Slider(plt.axes([0.25, startCoord + step * (iii * 5 + 1), 0.65, 0.03]), 'theta' + str(iii + 1), -180, 180,
                   valinit=0, valstep=0.5))
        rSliders.append(
            Slider(plt.axes([0.25, startCoord + step * (iii * 5 + 2), 0.65, 0.03]), 'r' + str(iii + 1), -10, 10,
                   valinit=0, valstep=0.5))
        alphaSliders.append(
            Slider(plt.axes([0.25, startCoord + step * (iii * 5 + 3), 0.65, 0.03]), 'alpha' + str(iii + 1), -180, 180,
                   valinit=0, valstep=0.5))

    sliders = [dSliders, thetaSliders, rSliders, alphaSliders]
    for iii in range(num):
        dSliders[iii].on_changed(lambda val: update(val, sliders, plotfig))
        thetaSliders[iii].on_changed(lambda val: update(val, sliders, plotfig))
        rSliders[iii].on_changed(lambda val: update(val, sliders, plotfig))
        alphaSliders[iii].on_changed(lambda val: update(val, sliders, plotfig))

    return sliders


def generateTfchain(ds, thetas, rs, alphas):
    tfmats = []
    initialTF = np.identity(4)
    finalTF = initialTF
    tfchain = [initialTF]
    for i in range(num):
        tfmats.append(tfmat(ds[i], thetas[i], rs[i], alphas[i], 1))

    for ii in range(len(tfmats)):
        finalTF = finalTF @ tfmats[ii].T
        tfchain.append(finalTF)
    return tfchain


def plotTfchain(tfchain, ax):
    points = []
    for tf in tfchain:
        points.append(plotTF(ax, tf))

    for i in range(len(points)-1):
        plotLine(ax, points[i], points[i+1], "k")
    ax.set_aspect('equal', 'box')


def update(val, sliders, fig):
    ax = fig.axes[0]
    dSliders = sliders[0]
    thetaSliders = sliders[1]
    rSliders = sliders[2]
    alphaSliders = sliders[3]
    ax.clear()
    ds = []
    thetas = []
    rs= []
    alphas = []
    for i in range(len(sliders[0])):
        ds.append(dSliders[i].val)
        thetas.append(thetaSliders[i].val)
        rs.append(rSliders[i].val)
        alphas.append(alphaSliders[i].val)
    tfchain = generateTfchain(ds, thetas, rs, alphas)
    plotTfchain(tfchain, ax)
    fig.canvas.draw_idle()


if __name__ == "__main__":
    num = 5
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sliders = generateSliders(num, fig)
    ds = range(num)
    thetas = [0] * num
    rs = [0] * num
    alphas = [0] * num
    tfchainOuter = generateTfchain(ds, thetas, rs, alphas)
    plotTfchain(tfchainOuter, ax)
    plt.show()
