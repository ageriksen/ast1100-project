import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import numpy as np

def plot():
    def update_lines(num, dataLines, lines) :
        for line, data in zip(lines, dataLines) :
            line.set_data(data[0:2, num-1:num])
            line.set_3d_properties(data[2,num-1:num])
        return lines

    # Attach 3D axis to the figure
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    m = 100   # number of frames in the animation
    n = 25   # number of particles you want to animate
    N = 1000 # number of time steps in your data

    data = np.zeros([n, 3, N]) # this is the positions of the particles
    # to be animated. In this code it should be an array of shape (n, 3, N)
    
    # ...Fill out the data array here! data[p, :, i] is particle p, time step i.
    
    # creates animation data for all your different particles
    lines = [i for i in range(n)]
    for i in range(n):
        lines[i] = [ax.plot(data[i][0,0:1], 
        data[i][1,0:1], data[i][2,0:1], 'o')[0]]

    # Set the axes properties
    ax.set_xlim3d([0.0, 1.0])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0, 1.0])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, 1.0])
    ax.set_zlabel('Z')

    ax.set_title('3D random particles')

    # Creating the Animation object
    ani = [i for i in range(n)] 
    for i in range(n):
        ani[i] = animation.FuncAnimation(fig, 
        update_lines, m, fargs=([data[i]], lines[i]),
        interval=50, blit=False)
    plt.show()
    
plot()


from AST1100SolarSystem import AST1100SolarSystem
system = AST1100SolarSystem(1234)
system.massNeededCheck(1,2,3,4,5) 