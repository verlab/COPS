import numpy as np
import copy
from libs.tsp import two_opt, path_distance_for_circular, path_distance_non_circular
from libs.grafo import COPS
from collections import OrderedDict
import time


def my_print(msg, index_print=0):
    if index_print == 1:
        print(msg)
    elif index_print == 2:
        print(msg)


class TabuSearchCOPS(COPS):
    def __init__(self, csop_class):
        self.tabu_alfa = 2
        self.beta = 4
        self.max_initial_solution_attempts = 10
        self.iterations_without_improvement = 0
        self.max_iterations_without_improvement = 10
        self.iterations_to_change_final_set = 0

        self.csop = csop_class
        self.array_profit = np.array(csop_class.profit)
        self.matrix_dist = csop_class.matrix_dist
        self.matrix_dist_int = self.matrix_dist#.astype(int)

        # VERTICES
        """ self.array_vertex[v][x]
            [v] -> index of the vertex
            [x] -> [0] -> [x, y]
                   [1] -> if vertex is being visited (1 visited, 0 otherwise) 
                   [2] -> list of clusters who this vertex belongs it
        """
        self.array_vertex = np.array([np.array([np.array(x), 0, []], dtype=object) for x in csop_class.list_vertex],
                                     dtype=object)

        # SETS
        """  """
        self.array_sets = np.array([np.array(x) for x in csop_class.list_clusters], dtype=object)
        self.num_sets = len(self.array_sets)

        # CLUSTERS
        """ self.array_clusters[c][x] 
            [c] -> index of the cluster 
            [x] -> [0] -> list of vertex who belongs to this cluster
                   [1] -> list of sets who this cluster belongs it
                   [2] -> LONG_TERM_MEMORY ( > 0 in solution) (<= 0 not in solution)
                   [3] -> 0 if only belong to start or end sets
                   [4] -> Aspiration level (the best profit obtained all times this cluster was in solution path)
                   [5] -> 1 this cluster can be inserted in solution 0 otherwise
                   [6] -> index (profit / number_of_cluster)  
        """
        self.array_clusters = np.array(
            [np.array([np.array(csop_class.list_subgroups[x]), [], 0, 0, 0, 1, self.array_profit[x]/len(csop_class.list_subgroups[x])], dtype=object) for x in range(len(csop_class.list_subgroups))], dtype=object)
        #my_print(f"array_cluster {self.array_clusters}", index_print=1)
        self.num_clusters = len(self.array_clusters)
        # Note that the same cluster can belong to more than one set.
        #            This loop says about each cluster which sets it belongs to.
        #            Ex: self.array_clusters[c][1] = [2,3] means that cluster c belongs to sets 2 and 3 """
        for s in range(self.num_sets):
            for c in self.array_sets[s]:
                self.array_clusters[c][1].append(s)
                # self.cluster_match_sets[c].append(s)

        # Note that the same vertex can belong to more than one set.
        #            This loop says about each vertex which clusters it belongs to.
        for c in range(self.num_clusters):
            for v in self.array_clusters[c][0]:
                self.array_vertex[v][2].append(c)

        self.start_set = csop_class.start_cluster
        self.end_set = csop_class.end_cluster

        self.start_cluster = self.array_sets[self.start_set][0]  # The start cluster is inside start set
        self.array_clusters[self.start_cluster][2] = 1  # long term memory (>=1 in solution)

        """ Note that: It's possible that end set contain more than one cluster who finishes the path,
                                but each end cluster must contain only one cluster."""
        self.end_cluster = np.random.choice(self.array_sets[self.end_set])
        self.best_end_cluster = copy.deepcopy(self.end_cluster)
        self.array_clusters[self.end_cluster][2] = 1  # long term memory (>=1 in solution)

        """ This dictionary contain:
            the clusters who belongs to any set except the start and end sets """
        # all_clusters contain all clusters except the clusters who belong only to the initial or final sets.
        # Note: the cluster who belong to initial or final sets could belong to another set
        # (in this case this cluster will be in the all_clusters variable)
        all_clusters = [c for s in range(self.num_sets) for c in self.array_sets[s] if
                        s != self.start_set and s != self.end_set]
        # the dictionary will eliminate the repeated clusters in all_clusters variable
        for c in all_clusters:
            self.array_clusters[c][3] = 1
        self.clusters_que_podem_ser_atualizados = np.array(
            [i for i in range(self.num_clusters) if self.array_clusters[i][3] == 1])

        """ These variable indicates whether a cluster or set is being visited
            0 -> unvisited  1 -> visited"""
        # self.clusters_visited = np.zeros(len(self.array_clusters))
        self.sets_visited = np.zeros(len(self.array_sets))
        # self.vertex_visited = np.zeros(len(self.array_vertex))

        self.solution = {"profit": 0,
                         "distance": 0,
                         "route": [],
                         "subgroups_visited": [],
                         # "sets_visited": [],
                         }

        self.best_solution = {"profit": 0,
                              "distance": 0,
                              "route": [],
                              "subgroups_visited": [],
                              # "sets_visited": [],
                              }

    def insertion_neighborhood(self, neighbor, inserted_cluster):
        need_a_solution = True

        n_tour, n_distance = self.tour_generation(neighbor)

        """ verify if this neighborhood is feasible """
        if n_distance < self.csop.t_max:
            """ verify if this neighborhood has the best profit or
                if the neighborhood has the same profit but less distance 
                NOTE: It will choose better paths when the profit is equal and the distance is less """
            n_profit = self.solution["profit"] + self.array_profit[inserted_cluster]

            if n_profit > self.solution["profit"] or (
                    n_profit == self.solution["profit"] and n_distance < self.solution["distance"]):

                """ update solution """
                self.solution["subgroups_visited"].append(inserted_cluster)
                self.solution["route"] = n_tour
                self.solution["distance"] = n_distance
                self.solution["profit"] = n_profit

                """ update sets visited"""
                for s in self.array_clusters[inserted_cluster][1]:  # self.cluster_match_sets[rand_cluster]:
                    self.sets_visited[s] = 1

                """ update long_term_memory """
                #for c in self.clusters_que_podem_ser_atualizados:
                for c in range(self.num_clusters):
                    if self.array_clusters[c][2] > 0:
                        self.array_clusters[c][2] += 1
                    else:
                        self.array_clusters[c][2] -= 1
                    self.array_clusters[inserted_cluster][2] = 1  # Inserted cluster should be a value 1

                """ update vertex inserted """
                for v in self.array_clusters[inserted_cluster][0]:
                    self.array_vertex[v][1] = 1

                """ update Aspiration level """
                for c in self.solution["subgroups_visited"]:
                    if n_profit > self.array_clusters[c, 4]:
                        self.array_clusters[c, 4] = n_profit

                need_a_solution = False
                my_print(f"CHANGE THE SOLUTION NEIGHBORHOOD {self.solution}")
        return need_a_solution

    def tour_update_remove_cluster(self, removed_cluster):
        n_tour = []
        n_distance = 0
        #my_print("TESTANDO tour_update_remove_cluster")

        """ eliminate a vertex if it belongs ONLY to the cluster who will be removed """
        eliminated_vertex = []
        for v in self.array_clusters[removed_cluster][0]:  # vertex in removed cluster
            can_eliminate_this_vertex = True
            for c in self.array_vertex[v][2]:  # clusters who this vertex belongs it (remember a vertex could belong to more than one cluster)
                if self.array_clusters[c][2] > 0:
                    if c != removed_cluster:
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
                    new_edge_init = self.solution["route"][t-1][0]
            else:
                if new_edge_init == -1:
                    init = self.solution["route"][t - 1][0]
                    end = self.solution["route"][t - 1][1]
                    n_tour.append((init, end))
                    n_distance += self.csop.matrix_dist[init][end]
                else:
                    init = new_edge_init
                    new_edge_init = -1
                    end = self.solution["route"][t][0]
                    n_tour.append((init, end))
                    n_distance += self.csop.matrix_dist[init][end]
        # treatment for the last edge
        t = len(self.solution["route"]) - 1
        if new_edge_init == -1:
            init = self.solution["route"][t][0]
            end = self.solution["route"][t][1]
            n_tour.append((init, end))
            n_distance += self.csop.matrix_dist[init][end]
        else:
            init = new_edge_init
            end = self.solution["route"][t][1]
            n_tour.append((init, end))
            n_distance += self.csop.matrix_dist[init][end]

        return n_tour, n_distance

    def removal_neighborhood(self, neighbor, removed_cluster):
        #need_a_solution = True

        """ melhorar """
        #n_tour, n_distance = self.tour_generation(neighbor)
        n_tour, n_distance = self.tour_update_remove_cluster(removed_cluster)

        """ verify if this neighborhood is feasible """
        #if n_distance < self.csop.t_max:
        """ update the aspiration level"""
        self.array_clusters[removed_cluster][4] = self.solution["profit"]

        #n_profit = self.solution["profit"] - self.array_profit[removed_cluster]

        """ update solution """
        self.solution["subgroups_visited"] = neighbor  # .remove(removed_cluster)
        self.solution["route"] = n_tour
        self.solution["distance"] = n_distance
        self.solution["profit"] -= self.array_profit[removed_cluster]

        """ update sets visited"""
        for s in self.array_clusters[removed_cluster][1]:  # self.cluster_match_sets[rand_cluster]:
            self.sets_visited[s] = 0

        """ update long_term_memory """
        #for c in self.clusters_que_podem_ser_atualizados:
        for c in range(self.num_clusters):
            if self.array_clusters[c][2] > 0:
                self.array_clusters[c][2] += 1
            else:
                self.array_clusters[c][2] -= 1
            self.array_clusters[removed_cluster][2] = 0  # removed cluster should be a value 0

        """ update vertex removed """
        for v in self.array_clusters[removed_cluster][0]:
            self.array_vertex[v][1] = 0

        #need_a_solution = False
        my_print(f"CHANGE THE SOLUTION NEIGHBORHOOD {self.solution}" )
        """ Always generates a feasible solution"""
        return False  #need_a_solution

    def insertion_criterion(self, index):
        #my_print(f"insertion_criterion", index_print=1)
        #print("index", index)
        criterion = self.array_clusters[index, 6] * self.array_clusters[index, 2]
        # remember: if we want to insert than the cluster are not in solution (long-term <= 0)
        min_value = np.argmin(criterion)
        chosen_cluster = index[min_value]
        #my_print(f"{chosen_cluster} - {self.array_clusters[index, 6]} - {self.array_clusters[index, 2]} - {criterion}", index_print=1)
        return chosen_cluster

    def generate_neighborhood(self):
        need_a_neighborhood = True

        visited_clusters = []
        non_tabu_remove = []
        old_visited = []

        unvisited_clusters = []
        non_tabu_insertion = []
        tabu_insertion = []

        for i in self.clusters_que_podem_ser_atualizados:
            long_term_memory = self.array_clusters[i][2]
            # update the lists for remove
            if long_term_memory > 0:
                visited_clusters.append(i)
                if long_term_memory > self.tabu_alfa:
                    non_tabu_remove.append(i)
                elif long_term_memory > self.beta:
                    old_visited.append(i)
            else:
                """ pode melhorar """
                """ it's a possible insertion if this cluster don't belong to a set who is in the solution """
                a = self.array_clusters[i][1]  # self.cluster_match_sets[i]
                none_set_was_visited = True
                for aa in a:
                    if self.sets_visited[aa] == 1:
                        none_set_was_visited = False
                        break
                # update the lists for insertion
                if none_set_was_visited:
                    unvisited_clusters.append(i)
                    if long_term_memory < -self.tabu_alfa:
                        non_tabu_insertion.append(i)
                    else:
                        tabu_insertion.append(i)

        # print("------ generated neighborhoods ---------------")
        # print("visited", visited_clusters)
        # print("unvisited", unvisited_clusters)
        # print("non_tabu_insert", non_tabu_insert)
        #print("LONG_TERM_MEMORY", self.array_clusters[:, 2])  # [c[2] for c in self.array_clusters])

        """ Non-Tabu Insertion """
        if any(non_tabu_insertion):
            """ Neighborhoods will be generated by a small modification of the current solution """
            neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
            chosen_cluster = np.random.choice(non_tabu_insertion)
            #chosen_cluster = self.insertion_criterion(non_tabu_insertion)
            #chosen_cluster = non_tabu_insertion[np.argmax([self.array_clusters[i][4] for i in non_tabu_insertion])]
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
                #choose the cluster with the aspiration level criterion
                #chosen_cluster = np.random.choice(tabu_insertion)
                #chosen_cluster = self.insertion_criterion(tabu_insertion)
                chosen_cluster = tabu_insertion[np.argmax([self.array_clusters[i][4] for i in tabu_insertion])]
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
            if any(visited_clusters):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_cluster = np.random.choice(visited_clusters)
                neighborhood.remove(chosen_cluster)
                need_a_neighborhood = self.removal_neighborhood(neighborhood, chosen_cluster)
                if not need_a_neighborhood:
                    my_print(f"Random Removal {chosen_cluster} {visited_clusters}")
                else:
                    my_print(f"discarded Random Removal {chosen_cluster} {visited_clusters}")

        """ Random Insertion """
        if need_a_neighborhood:
            if any(unvisited_clusters):
                """ Neighborhoods will be generated by a small modification of the current solution """
                neighborhood = copy.deepcopy(self.solution["subgroups_visited"])
                chosen_cluster = np.random.choice(unvisited_clusters)
                neighborhood.append(chosen_cluster)
                need_a_neighborhood = self.insertion_neighborhood(neighborhood, chosen_cluster)
                if not need_a_neighborhood:
                    my_print(f"Random Insertion {chosen_cluster} {unvisited_clusters}")
                else:
                    my_print(f"discarded Random Insertion {chosen_cluster} {unvisited_clusters}")

        if need_a_neighborhood:
            #for i in self.clusters_que_podem_ser_atualizados:
            for c in range(self.num_clusters):
                long_term_memory = self.array_clusters[i][2]
                if long_term_memory > 0:
                    long_term_memory += 1
                else:
                    long_term_memory -= 1
        else:
            self.choose_best_solution()

    def main(self):
        """ Compute an initial solution p0 """
        my_print("------ Initial Solution ---------------")
        self.initial_solution()
        my_print(f"solution {self.solution}")
        self.choose_best_solution()

        ''' melhorar parada - definir estratégia mais esperta '''
        my_print("------ Generated Neighborhoods ---------------")
        cont = 0
        #   while cont < 200:
        max_iterations = (len(self.array_clusters)) * 2
        while self.iterations_without_improvement < max_iterations:
            self.generate_neighborhood()
            cont += 1
            if not self.csop.circular_path:
                if self.iterations_to_change_final_set > self.max_iterations_without_improvement:
                    self.change_end_cluster()

        """ add start end end to solution subgroups_visited"""
        self.best_solution["subgroups_visited"].insert(0, self.start_cluster)
        #self.best_solution["subgroups_visited"].append([c for c in self.array_sets[self.end_set] if self.array_clusters[c][2] > 0][0])
        self.best_solution["subgroups_visited"].append(self.best_end_cluster)
        #########################################
        return self.best_solution

    def initial_solution(self):
        t1 = time.time()
        cv0 = []
        sets_visited = []
        # clusters_visited = []
        """ Order randomly the set array, without the initial and final sets """
        index_set = [i for i in range(len(self.array_sets)) if i != self.start_set and i != self.end_set]
        np.random.shuffle(index_set)

        sets_visited.append(self.start_set)
        # clusters_visited.append(self.start_cluster)
        profit = 0
        tour = []
        distance = 0
        """ For each set chose a cluster and try to find a plausible path """
        early_stop = self.max_initial_solution_attempts
        while any(index_set) and early_stop > 0:
            s = np.random.choice(index_set)
            #a = np.random.choice(self.array_sets[s])  # Chose randomly a cluster from set
            indexMaxProfitPerNunCluster = np.argmax([self.array_clusters[c][6] for c in self.array_sets[s]])  # Chose the cluster from set with the highest profit
            a = self.array_sets[s][indexMaxProfitPerNunCluster]
            if self.array_clusters[a][5]:  # if this cluster can be inserted
                cv0.append(a)

                initial_tour, initial_distance = self.tour_generation(cv0)

                """ test if the solution is plausible """
                if initial_distance > self.csop.t_max:
                    early_stop -= 1
                    cv0.pop(-1)
                    try:
                        index_set.remove(s)
                    except ValueError:
                        pass
                else:
                    early_stop = self.max_initial_solution_attempts
                    tour = initial_tour
                    distance = initial_distance
                    profit += self.array_profit[a]
                    # clusters_visited.append(a)
                    for visited in self.array_clusters[a][1]:  # list of sets who this cluster belongs it
                        try:
                            index_set.remove(visited)
                            sets_visited.append(visited)
                            for c in self.array_sets[visited]:
                                self.array_clusters[c][5] = 0  # can't be visited
                        except ValueError:
                            pass
                    my_print(f"clusters in solution {cv0} ")
                # my_print("-------------")

        """ For non_circular_path add the final cluster and set in their respective array """
        if not self.csop.circular_path:
            sets_visited.append(self.end_set)
            # clusters_visited.append(self.end_cluster)

        self.solution["route"] = tour
        self.solution["subgroups_visited"] = cv0
        # self.solution["subgroups_visited"].sort()
        # self.solution["sets_visited"] = sets_visited
        self.solution["distance"] = distance
        self.solution["profit"] = profit

        for i in cv0:
            self.array_clusters[i][2] = 1  # update LONG_TERM_MEMORY
            for v in self.array_clusters[i][0]:
                self.array_vertex[v][1] = 1  # indicates if a vertex is being visited
        self.array_vertex[self.array_clusters[self.end_cluster][0][0]][1] = 1  # all end and start cluster have only one vertex
        self.array_vertex[self.array_clusters[self.start_cluster][0][0]][1] = 1
        self.sets_visited[[i for i in sets_visited]] = 1

        """ update Aspiration level """
        for c in self.solution["subgroups_visited"]:
            if profit > self.array_clusters[c, 4]:
                self.array_clusters[c, 4] = profit

        my_print(f"LONG_TERM_MEMORY {[c[2] for c in self.array_clusters]}")
        my_print(f"sets_visited {self.sets_visited}")
        my_print(f"vertex_visited {[i[1] for i in self.array_vertex]}")

        tempoExec = time.time() - t1
        print("Tempo de execução Init Solution: {} segundos".format(tempoExec))

    def tour_generation(self, clusters, improvement_threshold=0.001):
        """ tour is generated with 2-opt technic """

        cv0 = clusters.copy()
        cv0.insert(0, self.start_cluster)

        """ select only wanted vertex """
        """ estudar o tempo de execução, com e sem a eliminação de elementos repetidos """
        #selected_index = [v for i in cv0 for v in self.array_clusters[i][0]]
        selected_index = list(OrderedDict.fromkeys([v for i in cv0 for v in self.array_clusters[i][0]]))  # will eliminate repeated elements
        #selected_vertex = np.array([self.array_vertex[i][0] for i in selected_index])

        """ greedy route based on the last route """
        final_selected_index = []
        if any(self.solution["route"]):
            for k in self.solution["route"]:
                index_v = k[0]
                final_selected_index.append(index_v)
                if index_v in selected_index:
                    selected_index.remove(index_v)
            for k in selected_index:
                near_from_k = np.argmin(self.matrix_dist_int[final_selected_index][:, k])
                final_selected_index.insert(near_from_k, k)

            selected_index = final_selected_index
            # print(selected_index)

        """ Note: the 2-opt solver needs a different treatment for a path that 
                                       ends at the same vertex it started """
        if not self.csop.circular_path:
            # cv0.append(self.end_subgroup)
            final_vertex = self.array_subgroups[self.end_subgroup][0][0]
            if final_vertex not in final_selected_index: final_selected_index.append(final_vertex)

        if False:  #self.csop.circular_path:
            """ calc the route using fast_tsp 
            OBS: fast_tsp does not specify the end depot """
            array_selected_index = np.array(selected_index)
            selected_matrix_dist = self.matrix_dist_int[array_selected_index[:, None], array_selected_index]
            route = fast_tsp.find_tour(selected_matrix_dist)
        else:
            """ calc the route from two-opt algorithm """
            selected_vertex = np.array([self.array_vertex[i][0] for i in selected_index])
            route = two_opt(selected_vertex, improvement_threshold, is_a_circular_path=self.csop.circular_path)

        """ calc the route from two-opt algorithm """
        #route = two_opt(selected_vertex, improvement_threshold, is_a_circular_path=self.csop.circular_path)
        real_route = [selected_index[i] for i in route]  # the real route is mapped from route

        """ define the edges from the found route """
        edges = [(real_route[i], real_route[i + 1]) for i in range(len(real_route) - 1)]
        if self.csop.circular_path:  # increase the last edge for a circular path
            edges.append((real_route[-1], real_route[0]))
        #distance = path_distance_for_circular(route, selected_vertex)  # only for conference
        """ calculate the path distance from edge distances matrix """
        d = sum([self.matrix_dist[edge[0]][edge[1]] for edge in edges])

        # print('selected_index', selected_index)
        # print('selected_vertex', selected_vertex)
        # print('route', route, real_route)
        # print('edges', edges)
        # print("another_distance", distance)
        # print("distance", d)

        return edges, d

    def choose_best_solution(self):
        if self.solution["profit"] > self.best_solution["profit"] or \
                (self.solution["profit"] == self.best_solution["profit"] and
                 self.solution["distance"] < self.best_solution["distance"]):
            self.best_solution = copy.deepcopy(self.solution)
            self.iterations_without_improvement = 0
            self.iterations_to_change_final_set = 0
            self.best_end_cluster = self.end_cluster
        else:
            self.iterations_without_improvement += 1
            self.iterations_to_change_final_set += 1

    def change_end_cluster(self):
        #my_print(f"-----------CHANGE THE END POINT {self.end_cluster}", index_print=2)

        # index of the min long-term-memory for the end clusters
        a = np.argmin([self.array_clusters[c][2] for c in self.array_sets[self.end_set]])

        # old end cluster long-term-memory update
        self.array_clusters[self.end_cluster][2] = 0  # long term memory (=0 removed from solution)
        self.array_vertex[self.array_clusters[self.end_cluster][0][0]][1] = 0  # vertex removed from solution

        # new end cluster will be the oldest it was in solution
        self.end_cluster = self.array_sets[self.end_set][a]
        self.array_clusters[self.end_cluster][2] = 1  # long term memory (>=1 in solution)
        self.array_vertex[self.array_clusters[self.end_cluster][0][0]][1] = 1

        #my_print(f" changed to {self.end_cluster}", index_print=2)

        # reboot iterations counter
        self.iterations_to_change_final_set = 0






