import numpy as np
# For movies:
#import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def main():

   # For movies:
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    x = np.loadtxt("x")
    y = np.loadtxt("y")

    numpoints = x.shape[1]
    numframes = x.shape[0]
    c = np.random.random((1, numpoints))
    my_list = []
    for _ in range(0,numpoints//2):
        my_list.append("deeppink")

    for _ in range(0,numpoints//2):
        my_list.append("yellow")

    my_list[0] = "lime"
    my_list[numpoints//2] = "red"

    fig = plt.figure()

    scat = plt.scatter(x[0][:], y[0][:], c=my_list, s=2e1)

    plt.axis('equal')
    plt.xlim(-2e18, 2e18)
    plt.ylim(-2e18, 2e18)

    ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes), interval=1, fargs=(x,y, scat))

    # For movies:                            
    #ani.save('animation.mp4', writer=writer)

def update_plot(i, a, b, scat):
    data = np.vstack((a[i][:],b[i][:])).T
    scat.set_offsets(data)
    return scat,

def energy():
    E = np.loadtxt("energy")
    plt.figure()
    plt.plot(E/E[0], 'o-', ms=4, lw=2)
    plt.title('Realtive total energy')
    plt.xlabel('iteration')
    plt.ylabel(r'$E_{tot}(t)/E_{tot}(0)$')


if __name__ == '__main__':
    energy()
    main()
    plt.show()
