from matplotlib import pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
#import matplotlib.ticker as ticker
from mpl_toolkits import mplot3d


def plot_edges(points, edges):
    lc = LineCollection(points[edges])
    fig = plt.figure()
    plt.gca().add_collection(lc)
    plt.xlim(points[:,0].min(), points[:,0].max())
    plt.ylim(points[:,1].min(), points[:,1].max())
    plt.plot(points[:,0], points[:,1], 'ro')
    #fig.savefig('full_figure.png')
    plt.show()


def plot(points, clusters, path):
    plt.figure(figsize=(6, 6))
    #plt.figure(figsize=(6, 12))
    #ax = fig.add_axes([0, 5, 0, 5])
    ax = plt.axes()
    colors = ['green', 'red', 'black', 'purple', 'cyan', 'magenta', 'yellow', 'blue', 'grey']
    for c in range(len(clusters)):
        xdata = [points[p][0] for p in clusters[c]]
        ydata = [points[p][1] for p in clusters[c]]
        #zdata = [points[p][2] for p in clusters[c]]
        ax.scatter(xdata, ydata, color=colors[c], label=str(c))  # , c=xdata, cmap='Greens')
    for segment in path:
        x = [points[segment[0]][0], points[segment[1]][0]]
        y = [points[segment[0]][1], points[segment[1]][1]]
        #z = [points[segment[0]][2], points[segment[1]][2]]
        ax.plot(x, y, color='black')  # label='parametric curve')
    #ax.set_xticks(ax.get_xticks()[::2])
    plt.ylim(-5, 5)
    plt.xlim(-5, 5)
    plt.title('Path')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.yticks(range(-5, 5, 2))  # alterar escala do eixo
    plt.xticks(range(-5, 5, 2))  # alterar escala do eixo
    plt.savefig('saved2D')
    plt.show()
    #x = [a[0] for a in points]
    #y = [a[1] for a in points]
    #plt.plot(x, y, 'ok')
    #plt.show()


def plot3D(points, clusters, path):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    colors = ['green', 'red', 'black', 'purple', 'cyan', 'magenta', 'yellow', 'blue', 'grey']
    for c in range(len(clusters)):
        xdata = [points[p][0] for p in clusters[c]]
        ydata = [points[p][1] for p in clusters[c]]
        zdata = [points[p][2] for p in clusters[c]]
        ax.scatter3D(xdata, ydata, zdata, color=colors[c], label=str(c))  #, c=xdata, cmap='Greens')
    for segment in path:
        x = [points[segment[0]][0], points[segment[1]][0]]
        y = [points[segment[0]][1], points[segment[1]][1]]
        z = [points[segment[0]][2], points[segment[1]][2]]
        ax.plot(x, y, z, color='black')  #label='parametric curve')
    plt.savefig('saved')
    plt.show()


#vertex = np.array([(0, 0), (1, 1), (1, 2), (2, 2), (2, 30)])
#visited = np.array([(0, 3), (0, 4), (1, 2), (1, 3), (2, 4)])

#plot(vertex)
#plot_edges(vertex, visited)

'''
aux1 = np.array(range(5))
for i in range(5):
    aux = aux1[aux1 != i]
    print(aux)
    print([((i, j) if j>i else (j,i)) for j in aux])
    #for j in aux: a = sum((i, j) if j > i else (j, i))
    #print(a)


#print([(a,b) for a in range(5) for b in range(5)])
'''