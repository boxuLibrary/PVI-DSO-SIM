import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')

x = []
y = []
z = []
with open("../bin/all_points.txt", "r") as data:
    for line in data.readlines():
        line = line.strip('\n').split(" ")
        x.append(float(line[0]))
        y.append(float(line[1]))
        z.append(float(line[2]))


position = []
tx_index = 5
position = np.loadtxt("../bin/cam_pose.txt", usecols = (tx_index, tx_index + 1, tx_index + 2))

ax.scatter(x, y, z,c='b')

ax.plot(position[:,0], position[:,1], position[:,2],'g', label='gt')
ax.plot([position[0,0]], [position[0,1]], [position[0,2]], 'r.', label='start')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(-14, 14)
ax.set_ylim(-14, 14)
ax.set_zlim(-1, 5)

    
plt.show()
