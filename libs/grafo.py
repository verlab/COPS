#import networkx as nx
import alphashape as alphashape
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
import numpy as np
from descartes import PolygonPatch


class COPS:
    def __init__(self):
        self.list_vertex = []
        self.list_subgroups = []
        self.list_clusters = []
        #self.g = nx.Graph()
        self.dimension = 0
        self.n_subgroups = 0
        self.profit = []

        self.start_cluster = 0
        self.end_cluster = 0
        self.t_max = 0
        self.circular_path = True

    def make_graph(self):
        print('vertex', self.list_vertex)
        print('list_subgroups', self.list_subgroups)
        print('list_clusters', self.list_clusters)
        print('profit', self.profit)
        #self.g.add_nodes_from(range(self.dimension))
        if self.edge_weight_type == 'CEIL_2D':
            self.euclidean_2D()

        #print(self.g.nodes())
        #print(self.g.edges())
        #print(nx.get_edge_attributes(self.g, "weight"))

        #return self.g

    def euclidean_2D(self):
        for i in range(self.dimension-1):
            for j in range(i+1, self.dimension):
                a = np.array(self.list_vertex[i])
                b = np.array(self.list_vertex[j])
                dist = np.linalg.norm(a-b)
                self.matrix_dist[i][j] = dist
                self.matrix_dist[j][i] = dist
                #self.g.add_edge(i, j, weight=dist)

    def draw_2D(self, path=[], legend=[], fill_cluster=False, fill_set=False, name='saved2D', save_img=False):
        sets = self.list_clusters
        clusters = self.list_subgroups
        points = np.array(self.list_vertex)
        profit = self.profit

        """ calculate a nice value for the dilate_hull """
        min_x = np.min(points[:, 0])
        max_x = np.max(points[:, 0])
        dilate_size = (max_x - min_x) / 45

        minrew = min(profit)
        maxrew = max(profit)
        mycmap = plt.cm.get_cmap('RdYlBu_r')
        cNorm = mpl.colors.Normalize(vmin=minrew, vmax=maxrew + 0.1 * (maxrew - minrew))
        mycmapScalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=mycmap)

        fig = plt.figure(figsize=(6,6))
        ax = plt.axes()
        color = iter(cm.rainbow(np.linspace(0, 1, self.n_subgroups)))
        #color = mycmapScalarMap.to_rgba(profit[c])
        # see markers in https://matplotlib.org/3.1.1/api/markers_api.html
        mark = [".","o","v","^","<",">","1","2","3","4","8","s","p","P","*","h","H","+","x","X","D","d","|","_",0,1,2,3,4,5,6,7,8,9,10,11]
        for c in range(len(clusters)):
            xdata = [points[p][0] for p in clusters[c]]
            ydata = [points[p][1] for p in clusters[c]]
            # zdata = [points[p][2] for p in clusters[c]]
            label = f"Cluster {c} v{[i for i in clusters[c]]}"
            #this_color = next(color)
            this_color = mycmapScalarMap.to_rgba(profit[c])
            #mark_index = c % len(mark)
            mark_index = 1
            black = "#080c0f"
            ax.scatter(xdata, ydata, color=black, label=label, marker=mark[mark_index], alpha=0.8, zorder=2)
            if fill_cluster:
                #xdata = [points[p][0] for c in sets[s] for p in clusters[c]]
                #ydata = [points[p][1] for c in sets[s] for p in clusters[c]]

                ppp = ([[points[p][0], points[p][1]] for p in clusters[c]])
                concave_hull = alphashape.alphashape(ppp, 0.)
                color = mycmapScalarMap.to_rgba(profit[c])
                dilated = concave_hull.buffer(dilate_size * 0.6)
                #BLUE = '#6699cc'
                black = "#080c0f"
                white = "#f5f6f7"
                patch1 = PolygonPatch(dilated, fc=this_color, ec=this_color, alpha=0.9, zorder=1)
                #patch1 = PolygonPatch(dilated, fc=white, ec=black, alpha=0.2, zorder=0)
                ax.add_patch(patch1)

        if fill_set:
            for s in range(len(sets)):
                ppp = ([[points[p][0], points[p][1]] for c in sets[s] for p in clusters[c]])
                concave_hull = alphashape.alphashape(ppp, 0.0)
                #color = mycmapScalarMap.to_rgba(profit[c])
                dilated = concave_hull.buffer(dilate_size * 1.2)
                #BLUE = '#6699cc'
                black = "#000000"  #"#080c0f"
                white = "#f5f6f7"
                #patch1 = PolygonPatch(dilated, fc=this_color, ec=this_color, alpha=0.2, zorder=0)
                patch1 = PolygonPatch(dilated, fc=white, ec=black, linewidth=2, alpha=0.3, zorder=0)
                ax.add_patch(patch1)

        for segment in path:
            x = [points[segment[0]][0], points[segment[1]][0]]
            y = [points[segment[0]][1], points[segment[1]][1]]
            # z = [points[segment[0]][2], points[segment[1]][2]]
            ax.plot(x, y, color='blue', linewidth=3)  # label='parametric curve')  #'black'

        if legend[0] == "route":
            ax.legend()
            plt.title(f'Route {path}')
        elif legend[0] == "r":
            pass
        else:
            plt.title(f"{[i for i in legend]}")

        plt.xlabel('X')
        plt.ylabel('Y')
        #
        #plt.yticks(range(-5, 5, 2))  # alterar escala do eixo
        #plt.xticks(range(-5, 5, 2))  # alterar escala do eixo

        if True:
            plt.axis('off')
            nodes_w_rewards = np.zeros((len(profit), 3))
            sc = plt.scatter(nodes_w_rewards[:, 0], nodes_w_rewards[:, 1], c=profit, cmap=mycmap, alpha=1.0,
                             s=1, facecolor='black', lw=0.5)
            tick_font_size = 12
            legend_font_size = 15
            cbar_position = [0.20, 0.049, 0.6, 0.02]
            cbar_ax = fig.add_axes(cbar_position)
            cb = plt.colorbar(sc, cax=cbar_ax, orientation='horizontal')
            cb.ax.tick_params(labelsize=tick_font_size)
            cb.set_label('profit', labelpad=-45.0, y=0.8, fontsize=legend_font_size)

            # offset = 0.08
            #fig.subplots_adjust(left=-0.0035, right=1.035, top=1.07, bottom=0.0)
            #fig.subplots_adjust(left=0.0, right=1, top=1., bottom=0.08)
            ax.set_aspect(1)

        #plt.autoscale(enable=True)
        #fig.set_size_inches(10, 5)#(10.5, 18.5)
        print(name)
        if save_img:
            plt.savefig(name)
        plt.show()
        plt.close('all')

    def read_data(self, dataset):
        file1 = open(dataset, 'r')
        Lines = file1.readlines()

        # Strips the newline character
        Vertices = 0
        Subgroups = 0
        Clusters = 0
        n_clusters = 0
        for line in Lines:
            values = [f for f in line.split()]
            if values[0] == 'DIMENSION:':
                self.dimension = int(values[1])
                self.matrix_dist = np.zeros((self.dimension, self.dimension))
            if values[0] == 'TMAX:':
                self.t_max = float(values[1])
            if values[0] == 'START_CLUSTER:':
                self.start_cluster = int(values[1])
            if values[0] == 'END_CLUSTER:':
                self.end_cluster = int(values[1])
            if values[0] == 'CLUSTERS:':
                n_clusters = int(values[1])
            if values[0] == 'SUBGROUPS:':
                self.n_subgroups = int(values[1])
            if values[0] == 'EDGE_WEIGHT_TYPE:':
                self.edge_weight_type = values[1]

            if 0 < Vertices <= self.dimension:
                self.list_vertex.append([float(values[1]), float(values[2])])
                Vertices += 1
            if values[0] == 'NODE_COORD_SECTION:':
                Vertices = 1

            if 0 < Subgroups <= self.n_subgroups:
                self.list_subgroups.append([int(c) for c in values[2:]])
                self.profit.append(float(values[1]))
                Subgroups += 1
            if values[0] == 'GTSP_SUBGROUP_SECTION:':
                Subgroups = 1

            if 0 < Clusters <= n_clusters:
                self.list_clusters.append([int(s) for s in values[1:]])
                Clusters += 1
            if values[0] == 'GTSP_CLUSTER_SECTION:':
                Clusters = 1

        self.circular_path = (self.start_cluster == self.end_cluster)
        self.make_graph()

        file1.close()

    def print_edges_of_the_route_distances(self, path=[]):
        print("Edges distances")
        total = 0
        for p in path:
            dist = self.matrix_dist[p[0]][p[1]]
            total += dist
            print(f"{p} - {dist}")
            print(f"Total walked {total}")