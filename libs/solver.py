import random

import numpy as np
import copy
from libs.tsp import two_opt, path_distance_for_circular, path_distance_non_circular
from libs.grafo import COPS
from collections import OrderedDict
import time
import fast_tsp


def my_print(msg, index_print=0):
    if index_print == 10:
        print(msg)
    elif index_print == 20:
        print(msg)


class TabuSearchCOPS(COPS):
    def __init__(self, cops_class):
        self.cont_tour = 0
        self.t1 = time.time()
        self.total_tour = 0
        self.all_tested = []
        self.contador_teste_subgrupo = 0

        self.subgrupos_testados = []
        self.subgrupos_maiores_que_tmax = []

        self.tabu_alfa = 2
        self.beta = 4
        self.max_initial_solution_attempts = 10
        self.iterations_without_improvement = 0
        self.max_iterations_without_improvement = 10
        self.iterations_to_change_final_subgroup = 0
        self.time = time.time()

        self.cops = cops_class
        self.array_profit = np.array(cops_class.profit)
        self.matrix_dist = cops_class.matrix_dist
        #self.matrix_dist_int = self.matrix_dist.astype(int)

        self.solution_storage = {}

        # VERTICES
        """ self.array_vertex[v][x]
            [v] -> index of the vertex
            [x] -> [0] -> [x, y]
                   [1] -> if vertex is being visited (1 visited, 0 otherwise) 
                   [2] -> list of subgroups who this vertex belongs to
        """
        self.array_vertex = np.array([np.array([np.array(x), 0, []], dtype=object) for x in cops_class.list_vertex],
                                     dtype=object)
        self.num_vertex = len(self.array_vertex)

        # CLUSTERS
        """  """
        self.array_clusters = np.array([np.array(x) for x in cops_class.list_clusters], dtype=object)
        self.num_clusters = len(self.array_clusters)

        # SUBGROUPS
        """ self.array_subgroups[c][x] 
            [c] -> index of the subgroup 
            [x] -> [0] -> list of vertex who belongs to this subgroup
                   [1] -> list of clusters who this subgroup belongs to
                   [2] -> LONG_TERM_MEMORY ( > 0 in solution) (<= 0 not in solution)
                   [3] -> 0 if only belong to start or end clusters
                   [4] -> Aspiration level (the best profit obtained all times this subgroup was in solution path)
                   [5] -> 1 this subgroup can be inserted in solution 0 otherwise
                   [6] -> index (profit / number_of_subgroups)  
        """
        self.array_subgroups = np.array(
            [np.array([np.array(cops_class.list_subgroups[x]), [], 0, 0, 0, 1,
                       self.array_profit[x] / len(cops_class.list_subgroups[x])], dtype=object) for x in
             range(len(cops_class.list_subgroups))], dtype=object)
        self.num_subgroups = len(self.array_subgroups)
        # Note that the same subgroup can belong to more than one cluster.
        #            This loop says about each subgroup which cluster it belongs to.
        #            Ex: self.array_subgroups[s][1] = [2,3] means that subgroup s belongs to clusters 2 and 3 """
        for s in range(self.num_clusters):
            for c in self.array_clusters[s]:
                self.array_subgroups[c][1].append(s)

        # Note that the same vertex can belong to more than one cluster.
        #            This loop says about each vertex which subgroups it belongs to.
        for j in range(self.num_subgroups):
            for k in self.array_subgroups[j][0]:
                self.array_vertex[k][2].append(j)

        self.start_cluster = cops_class.start_cluster
        self.end_cluster = cops_class.end_cluster

        self.start_subgroup = self.array_clusters[self.start_cluster][0]  # The start subgroup is inside start cluster
        self.array_subgroups[self.start_subgroup][2] = 1  # long term memory (>=1 in solution)

        """ Note that: It's possible that end cluster contain more than one subgroup who finishes the path,
                                but each end subgroup must contain only one subgroup."""
        self.end_subgroup = np.random.choice(self.array_clusters[self.end_cluster])
        self.best_end_subgroup = copy.deepcopy(self.end_subgroup)
        self.array_subgroups[self.end_subgroup][2] = 1  # long term memory (>=1 in solution)

        """ This dictionary contain:
            the subgroups who belongs to any cluster except the start and end clusters """
        # all_subgroups contain all subgroups except the subgroups who belong only to the initial or final clusters.
        # Note: the subgroup who belong to initial or final clusters could belong to another cluster
        # (in this case this subgroup will be in the all_subgroups variable)
        self.client_subgroups = np.array([s for c in range(self.num_clusters) for s in self.array_clusters[c] if
                                          c != self.start_cluster and c != self.end_cluster])
        # the dictionary will eliminate the repeated subgroups in all_subgroups variable
        for i in self.client_subgroups:
            self.array_subgroups[i][3] = 1
        self.subgroups_que_podem_ser_atualizados = np.array(
            [s for s in range(self.num_subgroups) if self.array_subgroups[s][3] == 1])

        """ These variable indicates whether a subgroup or cluster is being visited
            0 -> unvisited  1 -> visited"""
        self.clusters_visited = np.zeros(len(self.array_clusters))
        # self.vertex_visited = np.zeros(len(self.array_vertex))

        self.solution = {"profit": 0,
                         "distance": 0,
                         "route": [],
                         "subgroups_visited": [],
                         }

        self.best_solution = {"profit": 0,
                              "distance": 0,
                              "route": [],
                              "subgroups_visited": [],
                              }

    def insertion_neighborhood(self, neighbor, inserted_subgroup):
        need_a_solution = True
        n_tour, n_distance = self.tour_generation(neighbor)

        """ verify if this neighborhood is feasible """
        if n_distance <= self.cops.t_max:
            """ verify if this neighborhood has the best profit or
                if the neighborhood has the same profit but less distance 
                NOTE: It will choose better paths when the profit is equal and the distance is less """
            n_profit = self.solution["profit"] + self.array_profit[inserted_subgroup]

            if n_profit > self.solution["profit"] or (
                    n_profit == self.solution["profit"] and n_distance < self.solution["distance"]):

                """ storage the solution """
                #self.solution_storage[sorted(neighbor)] = self.solution

                """ update solution """
                self.solution["subgroups_visited"].append(inserted_subgroup)
                self.solution["route"] = n_tour
                self.solution["distance"] = n_distance
                self.solution["profit"] = n_profit

                """ update clustes visited"""
                for i in self.array_subgroups[inserted_subgroup][1]:
                    self.clusters_visited[i] = 1
                    for j in self.array_clusters[i]:
                        """ update subgroup visited """
                        self.array_subgroups[j][5] = 0  # can't be visited

                """ update long_term_memory """
                for j in range(self.num_subgroups):
                    if self.array_subgroups[j][2] > 0:
                        self.array_subgroups[j][2] += 1
                    else:
                        self.array_subgroups[j][2] -= 1
                    self.array_subgroups[inserted_subgroup][2] = 1  # Inserted subgroup should be a value 1

                """ update vertex inserted """
                for v in self.array_subgroups[inserted_subgroup][0]:
                    self.array_vertex[v][1] = 1

                """ update Aspiration level """
                for c in self.solution["subgroups_visited"]:
                    if n_profit > self.array_subgroups[c, 4]:
                        self.array_subgroups[c, 4] = n_profit

                need_a_solution = False
                my_print(f"CHANGE THE SOLUTION NEIGHBORHOOD {self.solution}")
        return need_a_solution

    def tour_update_remove_subgroup(self, removed_subgroup):
        n_tour = []
        n_distance = 0

        """ eliminate a vertex if it belongs ONLY to the subgroup who will be removed """
        eliminated_vertex = []
        for v in self.array_subgroups[removed_subgroup][0]:  # vertex in removed subgroup
            can_eliminate_this_vertex = True
            for c in self.array_vertex[v][
                2]:  # subgroups who this vertex belongs it (remember a vertex could belong to more than one subgroup)
                if self.array_subgroups[c][2] > 0:
                    if c != removed_subgroup:
                        can_eliminate_this_vertex = False
                        break
            if can_eliminate_this_vertex:
                eliminated_vertex.append(v)
        """ eliminate edges """
        new_edge_init = -1
        for t in range(1, len(self.solution["route"])):
            if any(self.solution["route"][t][0] == x for x in eliminated_vertex):
                eliminated_vertex.remove(self.solution["route"][t][0])
                if new_edge_init == -1:
                    new_edge_init = self.solution["route"][t - 1][0]
            else:
                if new_edge_init == -1:
                    init = self.solution["route"][t - 1][0]
                    end = self.solution["route"][t - 1][1]
                    n_tour.append((init, end))
                    n_distance += self.cops.matrix_dist[init][end]
                else:
                    init = new_edge_init
                    new_edge_init = -1
                    end = self.solution["route"][t][0]
                    n_tour.append((init, end))
                    n_distance += self.cops.matrix_dist[init][end]
        # treatment for the last edge
        t = len(self.solution["route"]) - 1
        if new_edge_init == -1:
            init = self.solution["route"][t][0]
            end = self.solution["route"][t][1]
            n_tour.append((init, end))
            n_distance += self.cops.matrix_dist[init][end]
        else:
            init = new_edge_init
            end = self.solution["route"][t][1]
            n_tour.append((init, end))
            n_distance += self.cops.matrix_dist[init][end]

        return n_tour, n_distance

    def removal_neighborhood(self, neighbor, removed_subgroup):
        # need_a_solution = True

        # n_tour, n_distance = self.tour_generation(neighbor)
        n_tour, n_distance = self.tour_update_remove_subgroup(removed_subgroup)

        """ verify if this neighborhood is feasible """
        # if n_distance < self.cops.t_max:
        """ update the aspiration level"""
        self.array_subgroups[removed_subgroup][4] = self.solution["profit"]

        # n_profit = self.solution["profit"] - self.array_profit[removed_subgroup]

        """ update solution """
        self.solution["subgroups_visited"] = neighbor  # .remove(removed_subgroup)
        self.solution["route"] = n_tour
        self.solution["distance"] = n_distance
        self.solution["profit"] -= self.array_profit[removed_subgroup]

        """ update clustes visited"""
        for i in self.array_subgroups[removed_subgroup][1]:
            self.clusters_visited[i] = 0
            for j in self.array_clusters[i]:
                """ update subgroup visited """
                self.array_subgroups[j][5] = 1  # can be visited

        """ update long_term_memory """
        for c in range(self.num_subgroups):
            if self.array_subgroups[c][2] > 0:
                self.array_subgroups[c][2] += 1
            else:
                self.array_subgroups[c][2] -= 1
            self.array_subgroups[removed_subgroup][2] = 0  # removed subgroup should be a value 0

        """ update vertex removed """
        for v in self.array_subgroups[removed_subgroup][0]:
            self.array_vertex[v][1] = 0

        # need_a_solution = False
        my_print(f"CHANGE THE SOLUTION NEIGHBORHOOD {self.solution}")
        """ Always generates a feasible solution"""
        return False  # need_a_solution

    def insertion_criterion(self, index):
        criterion = self.array_subgroups[index, 6] * self.array_subgroups[index, 2]
        # remember: if we want to insert than the subgroup are not in solution (long-term <= 0)
        min_value = np.argmin(criterion)
        chosen_cluster = index[min_value]
        # my_print(f"{chosen_cluster} - {self.array_subgroups[index, 6]} - {self.array_subgroups[index, 2]} - {criterion}", index_print=1)
        return chosen_cluster

    def original_generate_neighborhood(self):
        need_a_neighborhood = True

        visited_subgroups = []
        non_tabu_remove = []
        old_visited = []

        unvisited_subgroups = []
        non_tabu_insertion = []
        tabu_insertion = []

        for i in self.subgroups_que_podem_ser_atualizados:
            long_term_memory = self.array_subgroups[i][2]
            # update the lists for remove
            if long_term_memory > 0:
                visited_subgroups.append(i)
                if long_term_memory > self.tabu_alfa:
                    non_tabu_remove.append(i)
                elif long_term_memory > self.beta:
                    old_visited.append(i)
            else:
                """ it's a possible insertion if this subgroup don't belong to a cluster who is in the solution """
                a = self.array_subgroups[i][1]
                none_cluster_was_visited = True
                for aa in a:
                    if self.clusters_visited[aa] == 1:
                        none_cluster_was_visited = False
                        break
                # update the lists for insertion
                if none_cluster_was_visited:
                    unvisited_subgroups.append(i)
                    if long_term_memory < -self.tabu_alfa:
                        non_tabu_insertion.append(i)
                    else:
                        tabu_insertion.append(i)

        """ Non-Tabu Insertion """
        if any(non_tabu_insertion):
            """ Neighborhoods will be generated by a small modification of the current solution """
            neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
            chosen_cluster = np.random.choice(non_tabu_insertion)
            # chosen_cluster = self.insertion_criterion(non_tabu_insertion)
            # chosen_cluster = non_tabu_insertion[np.argmax([self.array_subgroups[i][4] for i in non_tabu_insertion])]
            neighborhood.append(chosen_cluster)
            need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_cluster)
            if not need_a_neighborhood:
                my_print(f"non_tabu_insertion {chosen_cluster} {non_tabu_insertion}")
            else:
                my_print(f"discarded Non-Tabu Insertion {chosen_cluster} {non_tabu_insertion}")

        """ Old Removal """
        if need_a_neighborhood:
            if any(old_visited):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_cluster = np.random.choice(old_visited)
                neighborhood.remove(chosen_cluster)
                need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_cluster)
                if not need_a_neighborhood:
                    my_print(f"old_Removal {chosen_cluster} {old_visited}")
                else:
                    my_print(f"discarded old_Removal {chosen_cluster} {old_visited}")

        """ Tabu Insertion """
        if need_a_neighborhood:
            if any(tabu_insertion):
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                # choose the subgroup with the aspiration level criterion
                # chosen_cluster = np.random.choice(tabu_insertion)
                # chosen_cluster = self.insertion_criterion(tabu_insertion)
                chosen_cluster = tabu_insertion[np.argmax([self.array_subgroups[i][4] for i in tabu_insertion])]
                neighborhood.append(chosen_cluster)
                need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_cluster)
                if not need_a_neighborhood:
                    my_print(f"tabu_insertion {chosen_cluster} {tabu_insertion}")
                else:
                    my_print(f"discarded tabu_insertion {chosen_cluster} {tabu_insertion}")

        """ Non-Tabu Removal """
        if need_a_neighborhood:
            if any(non_tabu_remove):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_cluster = np.random.choice(non_tabu_remove)
                neighborhood.remove(chosen_cluster)
                need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_cluster)
                if not need_a_neighborhood:
                    my_print(f"Non-Tabu Removal {chosen_cluster} {non_tabu_remove}")
                else:
                    my_print(f"discarded Non-Tabu Removal {chosen_cluster} {non_tabu_remove}")

        """ Random Removal """
        if need_a_neighborhood:
            if any(visited_subgroups):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_cluster = np.random.choice(visited_subgroups)
                neighborhood.remove(chosen_cluster)
                need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_cluster)
                if not need_a_neighborhood:
                    my_print(f"Random Removal {chosen_cluster} {visited_subgroups}")
                else:
                    my_print(f"discarded Random Removal {chosen_cluster} {visited_subgroups}")

        """ Random Insertion """
        if need_a_neighborhood:
            if any(unvisited_subgroups):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_cluster = np.random.choice(unvisited_subgroups)
                neighborhood.append(chosen_cluster)
                need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_cluster)
                if not need_a_neighborhood:
                    my_print(f"Random Insertion {chosen_cluster} {unvisited_subgroups}")
                else:
                    my_print(f"discarded Random Insertion {chosen_cluster} {unvisited_subgroups}")

        if need_a_neighborhood:
            for c in range(self.num_subgroups):
                long_term_memory = self.array_subgroups[i][2]
                if long_term_memory > 0:
                    long_term_memory += 1
                else:
                    long_term_memory -= 1
        else:
            self.choose_best_solution()

    def new_neighborhood(self, lista2):
        for key in self.solution_storage:
            if key == sorted(lista2):
                return False
        return True

    def generate_neighborhood(self):
        need_a_neighborhood = True

        visited_subgroups = []
        non_tabu_remove = []
        old_visited = []

        unvisited_subgroups = []
        non_tabu_insertion = []
        tabu_insertion = []

        for i in self.subgroups_que_podem_ser_atualizados:
            long_term_memory = self.array_subgroups[i][2]
            # update the lists for remove
            if long_term_memory > 0:
                visited_subgroups.append(i)
                if long_term_memory > self.tabu_alfa:
                    non_tabu_remove.append(i)
                elif long_term_memory > self.beta:
                    old_visited.append(i)
            else:
                """ it's a possible insertion if this subgroup don't belong to a cluster who is in the solution """
                a = self.array_subgroups[i][1]
                none_cluster_was_visited = True
                for aa in a:
                    if self.clusters_visited[aa] == 1:
                        none_cluster_was_visited = False
                        break
                # update the lists for insertion
                if none_cluster_was_visited:
                    unvisited_subgroups.append(i)
                    if long_term_memory < -self.tabu_alfa:
                        non_tabu_insertion.append(i)
                    else:
                        tabu_insertion.append(i)

        """ Nearest Insertion """
        if need_a_neighborhood:
            if any(unvisited_subgroups):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                if any(self.solution["subgroups_visited"]):
                    vertices_visited = list(
                        set([k for j in self.solution["subgroups_visited"] for k in self.array_subgroups[j][0]]))
                    aux_selected_client_vertex = [(k, j) for j in unvisited_subgroups for k in
                                                  self.array_subgroups[j][0]]
                    selected_client_vertex = [aux_selected_client_vertex[m][0] for m in
                                              range(len(aux_selected_client_vertex))]
                    # matrix with distances between vertices already visited (lines) and clients vertices that can be visited """
                    m = self.matrix_dist[vertices_visited][:, selected_client_vertex]
                    # index of the vertex with the min distance
                    index_min_dist = np.unravel_index(np.argmin(m, axis=None), m.shape)
                    index_min_dist_vertex = index_min_dist[1]
                    # index_vertex = selected_client_vertex[index_min_dist_vertex]
                    # subgroups that this vertex belongs
                    # subgroupshWichSharesThisVertex = [j for j in self.array_vertex[index_vertex][2]]
                    # subgroup with the higher profit that this vertex belongs
                    # index_ordered = np.argsort([self.array_subgroups[j][6] for j in subgroupshWichSharesThisVertex])
                    chosen_subgroup = aux_selected_client_vertex[index_min_dist_vertex][1]
                    neighborhood.append(chosen_subgroup)
                    need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_subgroup)

                    if not need_a_neighborhood:
                        my_print(f"Nearest Insertion {chosen_subgroup} {unvisited_subgroups}\n",
                                 index_print=2)
                    else:
                        my_print(f"discarded nearest Insertion {chosen_subgroup} {unvisited_subgroups}\n",
                                 index_print=2)

        """ Cheaper Removal """
        if need_a_neighborhood:
            if any(visited_subgroups):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                index = np.argsort([self.array_subgroups[j][6] for j in visited_subgroups])[:5]
                index = np.random.choice(index)
                chosen_subgroup = visited_subgroups[index]
                neighborhood.remove(chosen_subgroup)
                need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_subgroup)
                visited_subgroups.remove(chosen_subgroup)
                if chosen_subgroup in tabu_insertion: tabu_insertion.remove(chosen_subgroup)
                if chosen_subgroup in non_tabu_remove: non_tabu_remove.remove((chosen_subgroup))
                if not need_a_neighborhood:
                    my_print(f"Cheaper Removal {chosen_subgroup} {visited_subgroups}", index_print=1)
                else:
                    my_print(f"discarded Cheaper Removal {chosen_subgroup} {visited_subgroups}", index_print=1)
                need_a_neighborhood = True

        """ Non-Tabu Insertion """
        if need_a_neighborhood:
            if any(non_tabu_insertion):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_subgroup = np.random.choice(non_tabu_insertion)
                # chosen_subgroup = self.insertion_criterion(non_tabu_insertion)
                # chosen_subgroup = non_tabu_insertion[np.argmax([self.array_subgroups[i][4] for i in non_tabu_insertion])]
                neighborhood.append(chosen_subgroup)
                need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_subgroup)
                if not need_a_neighborhood:
                    my_print(f"non_tabu_insertion {chosen_subgroup} {non_tabu_insertion}", index_print=1)
                else:
                    my_print(f"discarded Non-Tabu Insertion {chosen_subgroup} {non_tabu_insertion}", index_print=1)

        """ Old Removal """
        if need_a_neighborhood:
            if any(old_visited):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_subgroup = np.random.choice(old_visited)
                neighborhood.remove(chosen_subgroup)
                need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_subgroup)
                if not need_a_neighborhood:
                    my_print(f"old_Removal {chosen_subgroup} {old_visited}", index_print=1)
                else:
                    my_print(f"discarded old_Removal {chosen_subgroup} {old_visited}", index_print=1)

        """ Tabu Insertion """
        if need_a_neighborhood:
            if random.random() > 0.5:
                if any(tabu_insertion):
                    neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                    # choose the subgroup with the aspiration level criterion
                    # chosen_subgroup = np.random.choice(tabu_insertion)
                    # chosen_subgroup = self.insertion_criterion(tabu_insertion)
                    chosen_subgroup = tabu_insertion[np.argmax([self.array_subgroups[i][4] for i in tabu_insertion])]
                    neighborhood.append(chosen_subgroup)
                    need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_subgroup)
                    if not need_a_neighborhood:
                        my_print(f"tabu_insertion {chosen_subgroup} {tabu_insertion}", index_print=1)
                    else:
                        my_print(f"discarded tabu_insertion {chosen_subgroup} {tabu_insertion}", index_print=1)

        """ Non-Tabu Removal """
        if need_a_neighborhood:
            if random.random() > 0.5:
                if any(non_tabu_remove):
                    """ Neighborhoods will be generated by a small modification of the current solution """
                    neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                    chosen_subgroup = np.random.choice(non_tabu_remove)
                    neighborhood.remove(chosen_subgroup)
                    need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_subgroup)
                    if not need_a_neighborhood:
                        my_print(f"Non-Tabu Removal {chosen_subgroup} {non_tabu_remove}", index_print=1)
                    else:
                        my_print(f"discarded Non-Tabu Removal {chosen_subgroup} {non_tabu_remove}", index_print=1)

        """ Random Removal """
        if need_a_neighborhood:
            if any(visited_subgroups):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_subgroup = np.random.choice(visited_subgroups)
                neighborhood.remove(chosen_subgroup)
                need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_subgroup)
                visited_subgroups.remove(chosen_subgroup)
                if not need_a_neighborhood:
                    my_print(f"Random Removal {chosen_subgroup} {visited_subgroups}", index_print=1)
                else:
                    my_print(f"discarded Random Removal {chosen_subgroup} {visited_subgroups}", index_print=1)
                need_a_neighborhood = True

        """ Random Insertion """
        if need_a_neighborhood:
            if any(unvisited_subgroups):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_subgroup = np.random.choice(unvisited_subgroups)
                neighborhood.append(chosen_subgroup)
                need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_subgroup)
                if not need_a_neighborhood:
                    my_print(f"Random Insertion {chosen_subgroup} {unvisited_subgroups}", index_print=1)
                else:
                    my_print(f"discarded Random Insertion {chosen_subgroup} {unvisited_subgroups}", index_print=1)

        """ most valuable Insertion """
        if need_a_neighborhood:
            if any(unvisited_subgroups):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                index = np.argsort([self.array_subgroups[j][6] for j in unvisited_subgroups])[::-1][:5]
                index = np.random.choice(index)
                chosen_subgroup = unvisited_subgroups[index]
                neighborhood.append(chosen_subgroup)
                need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_subgroup)
                if not need_a_neighborhood:
                    my_print(f"most valuable Insertion {chosen_subgroup} {unvisited_subgroups}", index_print=1)
                else:
                    my_print(f"discarded most valuable Insertion {chosen_subgroup} {unvisited_subgroups}",
                             index_print=1)

        if need_a_neighborhood:
            for i in range(self.num_subgroups):
                long_term_memory = self.array_subgroups[i][2]
                if long_term_memory > 0:
                    long_term_memory += 1
                else:
                    long_term_memory -= 1
        else:
            self.choose_best_solution()

    def formatar_lista(self, lista):
        # Converte todos os np.int64 para int
        return [int(item) for item in lista]

    def formatar_lista_de_listas(self, lista_de_listas):
        return [[int(item) for item in sublista] for sublista in lista_de_listas]

    def main(self):
        """ Compute an initial solution p0 """
        my_print("------ Initial Solution ---------------")
        # self.thirty_initial_solution()
        # self.another_initial_solution()
        self.initial_solution()
        #self.fifth_initial_solution()
        print(f"solution {self.solution['profit']}")

        self.choose_best_solution()
        # return self.best_solution
        my_print("------ Generated Neighborhoods ---------------")
        cont = 0
        #   while cont < 200:
        max_iterations = self.num_vertex * 2
        max_cont = max_iterations * 2
        while self.iterations_without_improvement < max_iterations and cont < max_cont:
            # self.time = time.time()
            #print('cont', cont, 'max_cont', max_cont, self.best_solution['profit'])
            #print(f'iterações sem melhora {self.iterations_without_improvement}')
            #print(f'subgrupos testados {len(self.subgrupos_testados)}')
            #print(f'subgrupos maiores que tmax {len(self.subgrupos_maiores_que_tmax)}')
            #print('--------------')
            self.generate_neighborhood()
            cont += 1
            if not self.cops.circular_path:
                if self.iterations_to_change_final_subgroup > self.max_iterations_without_improvement:
                    self.change_end_subgroup_old()
            '''if self.iterations_without_improvement % 10 == 0:
                print(f"--{cont}-{self.iterations_without_improvement}-{self.solution}")
                print("-----")
                print(f"best {self.best_solution}")'''

            # salvar dados no log
            '''dados = (
                f"Ciclo: {cont} - iterations_without_improvement {self.iterations_without_improvement}\n"
                f"Best profit {self.best_solution['profit']} distance {self.best_solution['distance']}\n"
                f"Runtime Solution: {(time.time() - self.t1)} seconds\n"
                f"Best subgroups_visited {self.formatar_lista(self.best_solution['subgroups_visited'])}\n "
                f"subgrupos_maiores_que_tmax {self.formatar_lista_de_listas(self.subgrupos_maiores_que_tmax)}\n"
                "--------\n"
            )
            with open("log_execucao.txt", 'a') as f:
                f.write(dados)'''

        '''cluster_visi = [i for j in self.best_solution["subgroups_visited"] for i in self.array_subgroups[j][1]]
        print("cluster_visi", np.sort(cluster_visi))
        for i in range(self.num_clusters):
            if i not in cluster_visi:
                print(i, "not visited")'''

        """ add start end end to solution subgroups_visited"""
        self.best_solution["subgroups_visited"].insert(0, self.start_subgroup)
        # self.best_solution["subgroups_visited"].append([c for c in self.array_clusters[self.end_subgroup] if self.array_subgroups[c][2] > 0][0])
        self.best_solution["subgroups_visited"].append(self.best_end_subgroup)
        #########################################

        #print("FINAL ANALISE - total tour = ", self.total_tour)
        #print('self.subgrupos_maiores_que_tmax', self.subgrupos_maiores_que_tmax)
        #print(self.all_tested)
        #for a in self.all_tested:
        #    print(f'{a},')
        return self.best_solution

    def initial_solution(self):
        t1 = time.time()
        cv0 = []
        clusters_visited = []
        """ Order randomly the clusters array, without the initial and final clusters """
        index_cluster = [i for i in range(len(self.array_clusters)) if
                         i != self.start_cluster and i != self.end_cluster]
        np.random.shuffle(index_cluster)

        clusters_visited.append(self.start_cluster)
        profit = 0
        tour = []
        distance = 0
        """ For each cluster chose a subgroup and try to find a plausible path """
        early_stop = self.max_initial_solution_attempts
        cont = 0
        while any(index_cluster) and early_stop > 0:
            #print(cont, early_stop)
            cont += 1
            i = np.random.choice(index_cluster)
            # a = np.random.choice(self.array_clusters[s])  # Chose randomly a subgroup from cluster
            # Chose the subgroup from cluster with the highest profit
            indexMaxProfitPerNunCluster = np.argmax([self.array_subgroups[j][6] for j in self.array_clusters[i]])
            a = self.array_clusters[i][indexMaxProfitPerNunCluster]
            if self.array_subgroups[a][5]:  # if this subgroup can be inserted
                cv0.append(a)

                initial_tour, initial_distance = self.tour_generation(cv0)

                """ test if the solution is plausible """
                if initial_distance > self.cops.t_max:
                    early_stop -= 1
                    cv0.pop(-1)
                    try:
                        index_cluster.remove(i)
                    except ValueError:
                        pass
                else:
                    early_stop = self.max_initial_solution_attempts
                    tour = initial_tour
                    distance = initial_distance
                    profit += self.array_profit[a]
                    for visited in self.array_subgroups[a][1]:  # list of clusters who this subgroup belongs it
                        try:
                            index_cluster.remove(visited)
                            clusters_visited.append(visited)
                            for c in self.array_clusters[visited]:
                                self.array_subgroups[c][5] = 0  # can't be visited
                        except ValueError:
                            pass
                    my_print(f"subgroups in solution {cv0} ")
                # my_print("-------------")

        """ For non_circular_path add the final subgroup and cluster in their respective array """
        if not self.cops.circular_path:
            clusters_visited.append(self.end_cluster)

        self.solution["route"] = tour
        self.solution["subgroups_visited"] = cv0
        self.solution["distance"] = distance
        self.solution["profit"] = profit

        for i in cv0:
            self.array_subgroups[i][2] = 1  # update LONG_TERM_MEMORY
            for v in self.array_subgroups[i][0]:
                self.array_vertex[v][1] = 1  # indicates if a vertex is being visited
        self.array_vertex[self.array_subgroups[self.end_subgroup][0][0]][
            1] = 1  # all end and start cluster have only one vertex
        self.array_vertex[self.array_subgroups[self.start_subgroup][0][0]][1] = 1
        self.clusters_visited[[i for i in clusters_visited]] = 1

        """ update Aspiration level """
        for c in self.solution["subgroups_visited"]:
            if profit > self.array_subgroups[c, 4]:
                self.array_subgroups[c, 4] = profit

        my_print(f"LONG_TERM_MEMORY {[s[2] for s in self.array_subgroups]}")
        my_print(f"clusters_visited {self.clusters_visited}")
        my_print(f"vertex_visited {[v[1] for v in self.array_vertex]}")

        tempoExec = time.time() - t1
        print("Runtime Init Solution: {} seconds".format(tempoExec))

    def greedy_tour(self, subgroups, improvement_threshold=0.001):
        """ tour is generated with 2-opt technic """

        cv0 = subgroups.copy()
        cv0.insert(0, self.start_subgroup)

        """ select only wanted vertex """
        # selected_index = [v for i in cv0 for v in self.array_subgroups[i][0]]
        selected_index = list(OrderedDict.fromkeys(
            [v for i in cv0 for v in self.array_subgroups[i][0]]))  # will eliminate repeated elements

        """ greedy route based on the last route """
        final_selected_index = []
        if any(self.solution["route"]):
            for k in self.solution["route"]:
                index_v = k[0]
                final_selected_index.append(index_v)
                if index_v in selected_index:
                    selected_index.remove(index_v)
            for k in selected_index:
                near_from_k = np.argmin(self.matrix_dist[final_selected_index][:, k])
                final_selected_index.insert(near_from_k, k)

            selected_index = final_selected_index
            #print(selected_index)

        """ Note: the 2-opt solver needs a different treatment for a path that 
                               ends at the same vertex it started """
        if not self.cops.circular_path:
            # cv0.append(self.end_subgroup)
            final_vertex = self.array_subgroups[self.end_subgroup][0][0]
            if final_vertex not in final_selected_index: final_selected_index.append(final_vertex)

        """ define the edges from the found route """
        edges = [(selected_index[i], selected_index[i + 1]) for i in range(len(selected_index) - 1)]
        if self.cops.circular_path:  # increase the last edge for a circular path
            edges.append((selected_index[-1], selected_index[0]))
        d = sum([self.matrix_dist[edge[0]][edge[1]] for edge in edges])

        return edges, d

    def comparar_subgrupos(self, novo_subgrupo):
        lista_existente = self.subgrupos_testados
        lista_grandes = self.subgrupos_maiores_que_tmax

        novo_set = set(novo_subgrupo)
        for subgrupo in lista_existente:
            sub_set = set(subgrupo)
            if sub_set == novo_set:
                return "Igual"

        for subgrupo in lista_grandes:
            sub_set = set(subgrupo)
            if sub_set.issubset(novo_set):
                return "Contém"
        return "Novo"

    def calcula_total_distancia(self, route, selected_index):
        # Mapeia o índice real dos vértices na rota
        real_route = [selected_index[i] for i in route]

        # Define as arestas da rota
        edges = [(real_route[i], real_route[i + 1]) for i in range(len(real_route) - 1)]
        if self.cops.circular_path:
            edges.append((real_route[-1], real_route[0]))

        # Calcula a distância total da rota
        edge_array = np.array(edges)
        total_distance = np.sum(self.matrix_dist[edge_array[:, 0], edge_array[:, 1]])

        return total_distance, edges

    def tour_generation(self, subgroups, improvement_threshold=0.001):
        """Gera um tour usando o algoritmo 2-opt"""

        #start_time = time.time()
        #log = lambda msg: print(f"[{(time.time() - start_time) * 1000:.2f} ms] {msg}")

        # Verifica subgrupo já testado
        ja_testado = self.comparar_subgrupos(subgroups)
        if True:#ja_testado == "Novo":
            #log("Teste de subgrupo.")
            # Cria sequência de subgrupos começando pelo início
            cv0 = [self.start_subgroup] + subgroups

            # Seleciona vértices únicos
            selected_index = list(dict.fromkeys(
                v for i in cv0 for v in self.array_subgroups[i][0]
            ))
            #log("Vértices únicos selecionados.")

            # Use previous solution as starting point for a greedy route

            if any(self.solution["route"]):
                final_selected_index = []
                for k in self.solution["route"]:
                    idx = k[0]
                    final_selected_index.append(idx)
                    if idx in selected_index:
                        selected_index.remove(idx)

                for k in selected_index:
                    distances = self.matrix_dist[final_selected_index][:, k]
                    near_pos = np.argmin(distances)
                    final_selected_index.insert(near_pos+1, k)

                selected_index = final_selected_index
                #log("Rota gulosa construída a partir da solução anterior.")

            # Garante que o destino esteja incluído no final, se necessário
            if not self.cops.circular_path:
                #end_vertex = self.array_subgroups[self.end_subgroup][0][0]
                #if end_vertex not in selected_index:
                #    selected_index.append(end_vertex)
                end_vertex = self.change_end_subgroup(selected_index)
                #end_vertex = self.array_subgroups[self.end_subgroup][0][0]
                selected_index.append(end_vertex)
            #log("Verificação de caminho circular concluída.")

            # Executa o algoritmo 2-opt
            if False:
                """ calc the route using fast_tsp 
                OBS: fast_tsp does not specify the end depot """
                array_selected_index = np.array(selected_index)
                selected_matrix_dist = self.matrix_dist[array_selected_index[:, None], array_selected_index]
                route = fast_tsp.greedy_nearest_neighbor(selected_matrix_dist)
                total_distance, edges = self.calcula_total_distancia(route, selected_index)
                # se distancia do greedy for superior ao tmax então rodar algoritmo mais sofisticado
                if total_distance > self.cops.t_max:
                    route = fast_tsp.find_tour(selected_matrix_dist, duration_seconds=0.5)

            else:
                # Cria vetor de coordenadas dos vértices selecionados
                selected_vertex = np.array([self.array_vertex[i][0] for i in selected_index])
                #log("Coordenadas dos vértices preparadas.")
                route = two_opt(selected_vertex, improvement_threshold, is_a_circular_path=self.cops.circular_path)

            #log("Algoritmo 2-opt finalizado.")

            total_distance, edges = self.calcula_total_distancia(route, selected_index)

            '''# Mapeia o índice real dos vértices na rota
            real_route = [selected_index[i] for i in route]

            # Define as arestas da rota
            edges = [(real_route[i], real_route[i + 1]) for i in range(len(real_route) - 1)]
            if self.cops.circular_path:
                edges.append((real_route[-1], real_route[0]))

            # Calcula a distância total da rota
            edge_array = np.array(edges)
            total_distance = np.sum(self.matrix_dist[edge_array[:, 0], edge_array[:, 1]])'''

            #log("Distância total calculada.")

            #log("Processamento completo.")
            #tempoExec = time.time() - self.t1
            #self.total_tour += (time.time() - start_time)
            #self.cont_tour += 1
            #self.all_tested.append([sorted(subgroups), total_distance, (time.time() - start_time)])
            #print('subgroups_visited', sorted(subgroups), total_distance)
            #print("Runtime Solution: {} seconds".format(tempoExec))
            #print(f"total tour generations {self.cont_tour}")

            if total_distance > self.cops.t_max:
                self.subgrupos_maiores_que_tmax.append(sorted(subgroups))
            else:
                self.subgrupos_testados.append(sorted(subgroups))
            #print('self.iterations_without_improvement', self.iterations_without_improvement)
            return edges, total_distance
        else:
            #log(f"Teste de subgrupo com subgrupo {ja_testado}.")
            #print('self.iterations_without_improvement', self.iterations_without_improvement)
            #self.contador_teste_subgrupo += 1
            #print(f"******qtd de subgrupos não testados pois já havia sido testado", self.contador_teste_subgrupo)
            # salvar dados no log
            """dados = (
                f"**********Teste de subgrupo com subgrupo {ja_testado} {subgroups}."
                "--------\n"
            )
            with open("log_execucao.txt", 'a') as f:
                f.write(dados)"""

            return [], self.cops.t_max + 1

    def choose_best_solution(self):
        if self.solution["profit"] > self.best_solution["profit"] or \
                (self.solution["profit"] == self.best_solution["profit"] and
                 self.solution["distance"] < self.best_solution["distance"]):
            self.best_solution = copy.deepcopy(self.solution)
            self.iterations_without_improvement = 0
            self.iterations_to_change_final_subgroup = 0
            self.best_end_subgroup = self.end_subgroup

            # self.next()
        else:
            self.iterations_without_improvement += 1
            self.iterations_to_change_final_subgroup += 1

    def change_end_subgroup(self, selected_index):
        """
        Selects a new endpoint vertex from the end cluster based on proximity to a the last vertex in the greedy route.

        Parameters
        ----------
        selected_index : list or array-like
            A list of selected vertex on the route. The last element (selected_index[-1] means the last vertex on the route based on the greedy part of the tour generate) is used as the reference 
            point for computing distances.

        Returns
        -------
        int
            The index of the vertex within the end cluster that is closest to the selected reference vertex.
        """
        end_subgroups = self.array_clusters[self.end_cluster]
        end_vertices = [self.array_subgroups[s][0][0] for s in end_subgroups]
        #print("end_vertices", end_vertices)
        #print("selected_index",selected_index)
        #print("self.matrix_dist[end_vertices]", self.matrix_dist[end_vertices])
        distances = self.matrix_dist[end_vertices][:, selected_index[-1]]
        near_pos = np.argmin(distances)
        #posicao = np.unravel_index(near_pos, distances.shape)
        return end_vertices[near_pos]
        

        

    def change_end_subgroup_old(self):
        # my_print(f"-----------CHANGE THE END POINT {self.end_subgroup}", index_print=2)

        # index of the min long-term-memory for the end subgroups
        a = np.argmin([self.array_subgroups[c][2] for c in self.array_clusters[self.end_cluster]])

        # old end subgroup long-term-memory update
        self.array_subgroups[self.end_subgroup][2] = 0  # long term memory (=0 removed from solution)
        self.array_vertex[self.array_subgroups[self.end_subgroup][0][0]][1] = 0  # vertex removed from solution

        # new end subgroup will be the oldest it was in solution
        self.end_subgroup = self.array_clusters[self.end_cluster][a]
        self.array_subgroups[self.end_subgroup][2] = 1  # long term memory (>=1 in solution)
        self.array_vertex[self.array_subgroups[self.end_subgroup][0][0]][1] = 1

        # my_print(f" changed to {self.end_subgroup}", index_print=2)

        # reboot iterations counter
        self.iterations_to_change_final_subgroup = 0
