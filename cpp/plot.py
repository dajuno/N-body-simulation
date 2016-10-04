import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

def main():

    # Set up formatting for the movie files
   # Writer = animation.writers['ffmpeg']
   # writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    x = np.loadtxt("x")
    y = np.loadtxt("y")
    
    numpoints = x.shape[1]
    numframes = x.shape[0]
    c = np.random.random((1, numpoints))

   # plt.ion()
    fig = plt.figure()
    
    scat = plt.scatter(x[0][:], y[0][:], c=c, s=1e2)
   # scat = plt.scatter(x[500][:], y[500][:], c='blue', s=1e2)
   # scat = plt.scatter(x[999][:], y[999][:], c='k', s=1e2)
    plt.axis('equal')
    plt.xlim(-2e18, 2e18)
    plt.ylim(-2e18, 2e18)
  #  plt.hold(False)
   # for i in range(numpoints):
   #     plt.scatter(x[i, :], y[i, :], c=c, s=1e2)
   #     print(i)
   #     time.sleep(0.3)
   #     plt.title(str(i))
   #     plt.show()
    ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes), interval=1, fargs=(x,y, scat))
                                  
  #  ani.save('animation.mp4', writer=writer)
    plt.show()

def update_plot(i, a, b, scat):
    data = np.vstack((a[i][:],b[i][:])).T
    scat.set_offsets(data)
    #print(i)
    return scat,

main()
