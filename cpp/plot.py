import numpy as np
# For movies:
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def main():

   # For movies:
   # Writer = animation.writers['ffmpeg']
   # writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    x = np.loadtxt("x")
    y = np.loadtxt("y")
    
    numpoints = x.shape[1]
    numframes = x.shape[0]
    c = np.random.random((1, numpoints))

    fig = plt.figure()
    
    scat = plt.scatter(x[0][:], y[0][:], c=c, s=1e1)
    
    plt.axis('equal')
    plt.xlim(-2e18, 2e18)
    plt.ylim(-2e18, 2e18)
    
    ani = animation.FuncAnimation(fig, update_plot, frames=range(numframes), interval=1, fargs=(x,y, scat))
      
    # For movies:                            
	#ani.save('animation.mp4', writer=writer)
    plt.show()

def update_plot(i, a, b, scat):
    data = np.vstack((a[i][:],b[i][:])).T
    scat.set_offsets(data)
    return scat,

main()
