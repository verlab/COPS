""" Thanks
https://www.gurobi.com/features/academic-named-user-license/    # install gurobi
https://support.gurobi.com/hc/en-us/community/posts/4404976974865-Orienteering-Problem
"""

# imports
import random
import time
import gurobipy as grb
import numpy as np
from itertools import chain, combinations
from gurobipy import Model
from libs.grafo import COPS


class IlpCOPS(COPS):
    opt_model: Model

    def __init__(self, cops_class):
        self.cops_class = cops_class

        self.matrix_dist = cops_class.matrix_dist
        self.circular_path = cops_class.circular_path
        self.max_distance = cops_class.t_max
        self.profit = np.array([np.array(x) for x in cops_class.profit], dtype=object)

        #self.vertex = np.array([np.array(x) for x in cops_class.list_vertex], dtype=object)
        # VERTICES
        """ self.array_vertex[v][x]
            [v] -> index of the vertex
            [x] -> [0] -> [x, y]
                   [1] -> if vertex is being visited (1 visited, 0 otherwise) 
                   [2] -> list of subgroups who this vertex belongs to
        """
        self.array_vertex = np.array([np.array([np.array(x), 0, []], dtype=object) for x in cops_class.list_vertex],
                                     dtype=object)
        self.qtd_v = len(cops_class.list_vertex)

        #self.subgroups = np.array([np.array(x) for x in cops_class.list_subgroups], dtype=object)
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
                       self.profit[x]], dtype=object) for x in
             range(len(cops_class.list_subgroups))], dtype=object)

        self.qtd_subgroups = len(self.array_subgroups)

        """ Clusters """
        self.clusters = np.array([np.array(x) for x in cops_class.list_clusters], dtype=object)
        self.end_cluster = cops_class.end_cluster
        self.num_clusters = len(self.clusters)

        # Note that the same subgroup can belong to more than one cluster.
        #            This loop says about each subgroup which cluster it belongs to.
        #            Ex: self.array_subgroups[s][1] = [2,3] means that subgroup s belongs to clusters 2 and 3 """
        for s in range(self.num_clusters):
            for c in self.clusters[s]:
                self.array_subgroups[c][1].append(s)

        # Note that the same vertex can belong to more than one cluster.
        #            This loop says about each vertex which subgroups it belongs to.
        for j in range(self.qtd_subgroups):
            for k in self.array_subgroups[j][0]:
                self.array_vertex[k][2].append(j)

        self.clientes = []
        self.aux1 = np.array(range(self.qtd_v))  # 11, 15
        self.end_vertices = [self.array_subgroups[i][0][j] for i in self.clusters[self.end_cluster] for j in
                             range(len(self.array_subgroups[i][0]))]

        self.opt_model = grb.Model()
        self.opt_model.setParam('Threads', 4)  # Define o número de threads desejado, neste caso, 4.
        self.opt_model.params.LogFile = 'log_.log'

        self.greedy_subgroup_tour_dist = np.zeros(self.qtd_subgroups)
        self.non_feasible_subgroups = []

        # Note that the same vertex can belong to more than one cluster.
        #            This loop says about each vertex which subgroups it belongs to.
        for j in range(self.qtd_subgroups):
            for v in self.array_subgroups[j][0]:
                self.array_vertex[v][2].append(j)

        #
        self.possible_end_start_vertex = {int(v) for s in self.clusters[self.end_cluster] for v in
                                          self.array_subgroups[s][0]}
        self.possible_end_start_vertex.add(0)

    def my_callback(self, model, where):
        if where == grb.GRB.Callback.MIPSOL:
            _x_edge = model.cbGetSolution(model._x_edge)
            _y = model.cbGetSolution(model._y)
            #print("_x_edge", _x_edge)
            #print("_y", _y)
            selected = grb.tuplelist((i, j) for i, j in model._x_edge.keys()
                                     if _x_edge[i, j] > 0.5)
            #print("selected", selected)
            incorrect_tour = self.sub_tour(selected)
            #print(incorrect_tour)
            if len(incorrect_tour) > 0:
                combinacoes = []
                #print("self.aux1", self.aux1)
                subset_barr = [x for x in self.aux1 if x not in incorrect_tour]
                #print("subset_barr", subset_barr)
                for i in subset_barr:
                    for j in incorrect_tour:
                        if i < j:
                            combinacoes.append([i, j])
                        elif i > j:
                            combinacoes.append([j, i])
                #print("combinacoes =", combinacoes)
                a = grb.quicksum(model._x_edge[i, j] for i, j in combinacoes)
                model.cbLazy(a >= model._y[incorrect_tour[0]])

            '''for subset in all_subsets(self.clientes):
                subset_barr = np.delete(self.aux1, subset)
                a = grb.quicksum((_x_edge[i, k]) if i < k else (_x_edge[k, i])
                                 for i in subset for k in subset_barr)
                for j in subset:
                    model.cbLazy(a >= _y[j])'''

    def sub_tour(self, edges):
        U = set([edges[i][j] for i in range(len(edges)) for j in range(2) if edges[i][j] in self.clientes])
        aux = set()
        aux.add(0)  # started vertex

        while len(aux) > 0:
            vertex = aux.pop()
            for edge in edges:
                if vertex in edge:
                    if edge[0] != vertex:
                        aux.add(edge[0])
                    else:
                        aux.add(edge[1])
                    edges.remove(edge)

            if vertex in U:
                U.remove(vertex)
        return list(U)  # if U is empty there is no sub_tour

    def create_model(self):
        pass

    def init_greedy_heu(self):
        '''
        Verifica o tsp de cada subgrupo e elimina aqueles que tiverem tsp do subgrupo maior que o tmax

        começa analisando do maior profit e adciona a solução inicial os subgrupos com maiores profits válidos,
        caso a distância do tsp do subgrupo + as duas menores edges deste subgrupo para algum outro do subgrupo

        + menor distancia dos vertices para o start e para o end vértice
        '''

        ''' pode calcular o centro de massa e comecar a buscar pelos clusters mais pertos'''

        for index_c in range(1, len(self.clusters)):
            vertex_init_solution = [0]
            if index_c != self.end_cluster:  # remember: start cluster is index 0, but end_cluster could be different
                list_subgroups_this_cluster = self.clusters[index_c]
                list_max_profit_index = np.flip(self.profit[list_subgroups_this_cluster].argsort())
                list_subgroups_max_profit = list_subgroups_this_cluster[list_max_profit_index]
                print(list_subgroups_max_profit)

                for index_s in list_subgroups_max_profit:
                    index_vertices = self.array_subgroups[index_s][0]
                    print(f'index vertices do subgrupo {index_s}: {index_vertices}')

                    greedy_subgroup_tour_dist = 0
                    visited_vertex = [index_vertices[0]]  # first
                    print('visited_vertex', visited_vertex)

                    current_vertex = index_vertices[0]
                    tour = [current_vertex]

                    while len(visited_vertex) < len(index_vertices):
                        next_vertex = None
                        min_distance = np.inf

                        for v in index_vertices:
                            if v not in visited_vertex and self.matrix_dist[current_vertex][v] < min_distance:
                                min_distance = self.matrix_dist[current_vertex][v]
                                next_vertex = v
                                print('next_vertex', next_vertex)

                        tour.append(next_vertex)
                        visited_vertex.append(next_vertex)
                        current_vertex = next_vertex
                        greedy_subgroup_tour_dist += min_distance

                    if greedy_subgroup_tour_dist < self.max_distance:
                        print('tour', tour, greedy_subgroup_tour_dist)

                        break
                    else:
                        ''' Define this subgroup unfeasible '''
                        print(f'Unfeasible subgroup {index_s}')
                        self.opt_model._z[index_s].setAttr('LB', 0)
                        self.opt_model._z[index_s].setAttr('UB', 0)

                self.greedy_subgroup_tour_dist = greedy_subgroup_tour_dist

    def init_greedy_heu_simple(self):
        """ Heurística de inicialização mais simples
            1) pega os subgrupos menos valiosos de cada cluster e os define como não pertencentes à solução inicial.
            2) não faz nada com o subgrupo mais valioso do cluster, já que o t_max pode ser insuficiente para a solução
            contendo todos os subgrupos mais valiosos de cada cluster.
        """
        for index_c in range(1, len(self.clusters)):
            if index_c != self.end_cluster:  # remember: start cluster is index 0, but end_cluster could be different
                list_subgroups_this_cluster = self.clusters[index_c]
                max_profit_index = self.profit[list_subgroups_this_cluster].argmax()
                subgroup_max_profit = list_subgroups_this_cluster[max_profit_index]

                list_subgroups_smaller_profits = list_subgroups_this_cluster[list_subgroups_this_cluster != subgroup_max_profit]
                for index_s in list_subgroups_smaller_profits:
                    self.opt_model._z[index_s].start = 0

                print(list_subgroups_this_cluster)
                print(max_profit_index)
                print(subgroup_max_profit)
                print(list_subgroups_smaller_profits)
                print('-------------')

    def init_greedy_heu_simple2(self, maxtime=100, display=0):
        """ Heurística de inicialização
            1) verifica a viabilidade de cada subgrupo individualmente
            1.a)  testa apenas subgrupos com mais de 1 vértice, já que os vértices individuais já foram testados.
        """
        # stats
        start_time = time.time()
        optimal = 0
        timelimit = 0
        # change parameters
        origt = self.opt_model.Params.timelimit
        origP = self.opt_model.Params.outputflag
        origG = self.opt_model.Params.mipgap
        self.opt_model.params.timelimit = maxtime
        self.opt_model.params.outputflag = display
        self.opt_model.params.mipgap = 0.005

        print("Testando subgrupos individualmente")

        # setar todos os subgrupos como zero
        for index_s in range(1, self.qtd_subgroups):
            if index_s not in self.clusters[self.end_cluster]:
                self.opt_model._z[index_s].lb = 0
                self.opt_model._z[index_s].ub = 0

        # escolher um subgrupo de cada vez e testar a viabilidade
        for index_s in range(1, self.qtd_subgroups):
            if len(self.array_subgroups[index_s][0]) > 1:
                if index_s not in self.non_feasible_subgroups and index_s not in self.clusters[self.end_cluster]:
                    print(f"Teste subgrupo {index_s}")
                    self.opt_model._z[index_s].lb = 1
                    self.opt_model._z[index_s].ub = 1
                    self.opt_model.update()

                    self.opt_model.optimize(lambda model, where: self.my_callback(model, where))
                    if self.opt_model.Status == grb.GRB.OPTIMAL:
                        optimal += 1
                    if self.opt_model.Status == grb.GRB.TIME_LIMIT:
                        timelimit += 1
                        print(f"time limit cont = {timelimit}")
                    if self.opt_model.getAttr('SolCount') == 0:
                        print(f'Failed in subgroup {index_s} tour')
                        self.non_feasible_subgroups.add(index_s)

                    self.opt_model._z[index_s].lb = 0
                    self.opt_model._z[index_s].ub = 0
                    self.opt_model.update()

        # voltar aos valores iniciais mantendo os subgrupos não viáveis como zero
        max_profit = 0
        index_subgroup_max_profit = -1
        for index_s in range(1, self.qtd_subgroups):
            if index_s not in self.non_feasible_subgroups and index_s not in self.clusters[self.end_cluster]:
                self.opt_model._z[index_s].lb = 0
                self.opt_model._z[index_s].ub = 1
                if self.profit[index_s] > max_profit:
                    max_profit = self.profit[index_s]
                    index_subgroup_max_profit = index_s

        # adicinar subgrupo válido de maior profit na solução inicial
        #print(f"start solution with the subgroup {index_subgroup_max_profit} with the higher profit {max_profit}")
        #if index_subgroup_max_profit > -1:
        #    self.opt_model._z[index_subgroup_max_profit].start = 1

        print(f"Final solution "
              f"optimal {optimal} timelimit {timelimit} {round(time.time() - start_time, 2)} seconds")
        # restore model parameter values
        self.opt_model.params.timelimit = origt
        self.opt_model.params.outputflag = origP
        self.opt_model.Params.mipgap = origG

        self.opt_model.update()

    def try_to_add_the_most_profitable_subgroups(self):
        """ caso """
        pass

    def try_to_add_the_closest_vertex(self, maxtime=100, display=0):
        """ insert a valid subgroup with the closest vertex to the actual rote. """
        # stats
        start_time = time.time()
        optimal = 0
        timelimit = 0
        # change parameters
        origt = self.opt_model.Params.timelimit
        origP = self.opt_model.Params.outputflag
        origG = self.opt_model.Params.mipgap
        self.opt_model.params.timelimit = maxtime
        self.opt_model.params.outputflag = display
        self.opt_model.params.mipgap = 0.005

        # list of vertices in the actual route
        all_client_vertices = set([v for v in range(self.qtd_v)]) - set(self.possible_end_start_vertex)

        invalid_subgroups = []  # for subgroups belonging to clusters already served
        vertices_in_route = [i for i in self.possible_end_start_vertex]
        clients_vertices_NAO_podem_ser_visitados = []
        clients_vertices_PODEM_ser_visitados = set(all_client_vertices) - set(clients_vertices_NAO_podem_ser_visitados)

        # escolher um subgrupo de cada vez e testar a viabilidade
        start_subgroups = []
        while len(clients_vertices_PODEM_ser_visitados) > 0:
            print(f"Testando ...{clients_vertices_PODEM_ser_visitados}")
            # find the closest vertex
            u, v, distancia = menor_aresta(self.possible_end_start_vertex, clients_vertices_PODEM_ser_visitados,
                                           self.matrix_dist)
            #print("print(u, v,distancia)", u, v, distancia)
            #print("lista de subgroups - ", self.array_vertex[v][2])

            index_s = self.array_vertex[v][2][0]  # esta pegando o primeiro subgrupo que o vértice pertence (alterar para o mais valioso)
            print(f"Testando subgrupo {index_s}")
            if index_s not in self.non_feasible_subgroups:
                self.opt_model._z[index_s].lb = 1
                self.opt_model._z[index_s].ub = 1
                self.opt_model.update()

                self.opt_model.optimize(lambda model, where: self.my_callback(model, where))
                if self.opt_model.Status == grb.GRB.OPTIMAL:
                    start_subgroups.append(index_s)
                    print(f"Add vertice {v} within subgroup {index_s}")

                    # remove cluster tested
                    for cluster_to_remove in self.array_subgroups[index_s][1]:
                        for subgroup_to_remove in self.clusters[cluster_to_remove]:
                            print(f"---removing vertices do subgrupo {subgroup_to_remove}")
                            for index_v in self.array_subgroups[subgroup_to_remove][0]:
                                if index_v in clients_vertices_PODEM_ser_visitados:
                                    clients_vertices_PODEM_ser_visitados.remove(index_v)
                if self.opt_model.Status == grb.GRB.TIME_LIMIT:
                    timelimit += 1
                    print(f"time limit cont = {timelimit}")

                    # remove subgroup tested
                    print(f"---removing vertices do subgrupo {index_s}")
                    for index_v in self.array_subgroups[index_s][0]:
                        clients_vertices_PODEM_ser_visitados.remove(index_v)

                if self.opt_model.getAttr('SolCount') == 0:
                    print(f'Failed to Add subgroup {index_s}')
                    self.opt_model._z[index_s].lb = 0
                    self.opt_model._z[index_s].ub = 0

                    # remove subgroup tested
                    print(f"---removing vertices do subgrupo {index_s}")
                    for index_v in self.array_subgroups[index_s][0]:
                        clients_vertices_PODEM_ser_visitados.remove(index_v)
            else:
                # remove subgroup tested
                print(f"---removing vertices do subgrupo {index_s}")
                for index_v in self.array_subgroups[index_s][0]:
                    clients_vertices_PODEM_ser_visitados.remove(index_v)

        # voltar aos valores iniciais mantendo os subgrupos não viáveis como zero
        for index_s in range(1, self.qtd_subgroups):
            if index_s not in self.non_feasible_subgroups and index_s not in self.clusters[self.end_cluster]:
                self.opt_model._z[index_s].lb = 0
                self.opt_model._z[index_s].ub = 1

        # start
        for i in start_subgroups:
            self.opt_model._z[i].start = 1
            print(f"Start subgroups {start_subgroups}")

        # restore model parameter values
        self.opt_model.params.timelimit = origt
        self.opt_model.params.outputflag = origP
        self.opt_model.Params.mipgap = origG
        self.opt_model.update()
        self.opt_model.update()

    def init_greedy_heu_insertion(self):
        """ Heurística de inicialização ...
            1)
        """

        subgrupos_not_used = list(self.non_feasible_subgroups)

        for index_c in range(1, len(self.clusters)):
            if index_c != self.end_cluster:  # remember: start cluster is index 0, but end_cluster could be different
                list_subgroups_this_cluster = self.clusters[index_c]
                possible_vertices = []
                for index_s in list_subgroups_this_cluster:
                    if index_s not in subgrupos_not_used:
                        for index_v in self.array_vertex[index_s]:
                            possible_vertices.append(index_v)

                possible_init_vertices = set(possible_vertices)

    def main(self):
        """MELHORAR"""

        # clientes = []
        for i in range(1, self.qtd_v):
            um_cliente = True
            for j in self.possible_end_start_vertex:
                if i == j:
                    um_cliente = False
                    break
            if um_cliente:
                self.clientes.append(i)

        # 2(i)
        '''w = {}
        for i in range(len(self.clusters)):
            w[i] = self.opt_model.addVar(vtype=grb.GRB.BINARY, name='w' + str(i))
        '''
        # 2(j)
        """ zi ∈ {0, 1}, ∀Si ∈ S """
        z = {}
        for i in range(self.qtd_subgroups):
            z[i] = self.opt_model.addVar(vtype=grb.GRB.BINARY, name='z' + str(i))

        # 2(k) + edge's cost
        """ ∑ texe ≤ Tmax """
        """ xe ∈ {0, 1}, ∀e ∈ E """
        t_cost = {}
        x_edge = {}

        for i in range(self.qtd_v):
            # Delete subgroups
            d = self.matrix_dist[0][i] + np.min(self.matrix_dist[i][self.end_vertices])
            if d > self.max_distance:
                print(f'minimum distance from route contain vertex {i} is {d} -> list subgroups eliminates {self.array_vertex[i][2]}')
                for index_s in self.array_vertex[i][2]:  # list of subgroups this vertex belongs to
                    z[index_s].lb = 0
                    z[index_s].ub = 0
                    self.non_feasible_subgroups.append(index_s)

            for j in range(i + 1, self.qtd_v):
                # if i != j:
                x_edge[i, j] = self.opt_model.addVar(vtype=grb.GRB.BINARY, name='edge_' + str(i) + '_' + str(j))
                t_cost[i, j] = self.matrix_dist[i][j]  # distance(vertex, i, j)

                if t_cost[i, j] == np.inf or t_cost[i, j] > self.max_distance:
                    # Delete this edge
                    #self.opt_model.addConstr(x_edge[i,j] == 0)
                    x_edge[i, j].lb = 0  # lower bound
                    x_edge[i, j].ub = 0  # upper bound
                    t_cost[i, j] = 0   # para evitar erros com valores inf -> não tem outra função já que esta aresta não existe

        self.non_feasible_subgroups = set(self.non_feasible_subgroups)
        print(f"qtd de subgrupos eliminados {len(self.non_feasible_subgroups)}")

        # 2(l)
        """ yj ∈ {0, 1}, ∀vj ∈ V """
        y = {}
        for i in range(self.qtd_v):
            y[i] = self.opt_model.addVar(vtype=grb.GRB.BINARY, name='y' + str(i))

        # constraints
        # 2(a)
        """ y0 = 1 """
        self.opt_model.addConstr(z[0] == 1)  # (2)
        self.opt_model.addConstr(y[0] == 1)  # (2)

        # 2(b)
        """ ∑ xe = 2yj ∀vj ∈ C\{C_start, C_end} """
        # aux1 = np.array(range(qtd_v))
        for j in self.clientes:
            aux = self.aux1[self.aux1 != j]
            self.opt_model.addConstr(grb.quicksum((x_edge[j, k] if k > j else x_edge[k, j]) for k in aux) == 2 * y[j])

        # 2(c)
        """ ∑ xe(s0) ≤ 2"""
        start_vertex = 0
        aux = self.aux1[self.aux1 != start_vertex]
        if self.circular_path:
            # start and end cluster
            self.opt_model.addConstr(
                grb.quicksum(
                    (x_edge[start_vertex, k] if k > start_vertex else x_edge[k, start_vertex]) for k in aux) == 2)
        else:
            # start cluster
            self.opt_model.addConstr(
                grb.quicksum(
                    (x_edge[start_vertex, k] if k > start_vertex else x_edge[k, start_vertex]) for k in aux) == 1)
            # end cluster
            '''aux = [item for item in self.aux1 if item not in self.end_vertices]
            for end_vertex in self.end_vertices:
                self.opt_model.addConstr(
                    grb.quicksum(
                        (x_edge[end_vertex, k] if k > end_vertex else x_edge[k, end_vertex]) for k in aux) <= 1)
'''
        # 2(d)
        """ ∑ texe ≤ Tmax """
        self.opt_model.addConstr(
            grb.quicksum(
                t_cost[i, j] * x_edge[i, j] for i in range(self.qtd_v) for j in range(i + 1, self.qtd_v)) <= self.max_distance)

        # 2(e)
        """ sub-tour constraint """
        """
        for subset in all_subsets(clientes):  # (6)
            qtd = len(subset)
            subset_barr = np.delete(aux1, subset)
            self.opt_model.addConstr(
                grb.quicksum((x_edge[i, j]) if i < j else (x_edge[j, i])
                             for i in subset for j in subset_barr) >= grb.quicksum(y[k] for k in subset) / qtd)
        """
        # 2(f)
        """ zi ≤ yj , ∀Si ∈ S, ∀vj ∈ Si """
        for i in range(1, self.qtd_subgroups):  # (7)
            # self.opt_model.addConstr(y[subgroups[i][j]] for j in range(len(subgroups[i])) >= z[i])
            for j in range(len(self.array_subgroups[i][0])):
                self.opt_model.addConstr(y[self.array_subgroups[i][0][j]] >= z[i])

        # 2(g)
        """ ∑zi ≤ w ∀Cg ∈ C """
        for i in range(len(self.clusters)):  # new constraint
            ''' não pode atender a dois grupos COMPLETOS ao mesmo tempo 
                mas pode atender a um grupo completo e a alguns indivíduos do outro
            '''
            self.opt_model.addConstr(grb.quicksum(z[self.clusters[i][j]] for j in range(len(self.clusters[i]))) <= 1)#w[i])

        # 2(h)
        """ ∑yi ≤ ∑Zi | Si | """
        ''' se atende a um grupo não pode atender a nenhum indivíduo do outro '''
        self.opt_model.addConstr(grb.quicksum(y[v] for v in range(self.qtd_v)) <= grb.quicksum(
            z[s] * len(self.array_subgroups[s][0]) for s in range(self.qtd_subgroups)))

        # 14, 15
        if not self.circular_path:
            """ novo 1 vertice do último superset atendido"""
            self.opt_model.addConstr(grb.quicksum(y[v] for s in self.clusters[self.end_cluster] for v in self.array_subgroups[s][0]) == 1)

            """  """
            self.opt_model.addConstr(grb.quicksum(
                x_edge[min(vvv, end), max(vvv, end)] for vvv in range(self.qtd_v) for s in self.clusters[self.end_cluster]
                for end in
                self.array_subgroups[s][0] if vvv != end) == 1)

        # STARTS
        #start_subgroups = [1,4,7,10,13,16]
        #for i in start_subgroups:
        #    z[i].start = 1

        """ Save additional data within model """
        self.opt_model._x_edge = x_edge  # edges
        self.opt_model._y = y  # vertices
        self.opt_model._z = z  # subgroups
        self.opt_model._t_cost = t_cost

        """ objective max profit """
        objective = grb.quicksum(self.profit[i] * self.opt_model._z[i] for i in range(self.qtd_subgroups))
        self.opt_model.ModelSense = grb.GRB.MAXIMIZE
        self.opt_model.setObjective(objective)
        self.opt_model.Params.lazyConstraints = 1
        self.opt_model.params.outputflag = 0  # verbose false

        usegreedy = True
        if usegreedy:
            self.opt_model.update()  # needed because we access variable attributes in heuristics
            #self.init_greedy_heu()
            #self.init_greedy_heu_simple()
            self.init_greedy_heu_simple2()
            #self.try_to_add_the_closest_vertex()

        # 1
        """ optimize """
        self.opt_model.optimize(lambda model, where: self.my_callback(model, where))
        # opt_model.optimize()

        '''
        status = opt_model.status
        if status == grb.Status.INFEASIBLE:    #opt_model.computeIIS()
            if opt_model.IISMinimal:
                print('IIS is minimal\n')
            else:
                print('IIS is not minimal\n')
                print('\nThe following constraint(s) cannot be satisfied:')
            for c in opt_model.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
            exit()

        '''

        try:
            solution = self.opt_model.getAttr('x', z)
        except grb.GurobiError:
            print("'\nThe following constraint(s) cannot be satisfied:'")
            """ Talvez deva retornar conjunto vazio"""
            exit()

        select = grb.tuplelist((i) for i in z.keys() if z[i].X > 0.5)
        # print('original points =', vertex)
        # print('objective', objective)
        # print('solução encontrada', solution)
        # print('Cost =', np.sum([t_cost[i, j]*x_edge[i, j] for j in range(qtd_v)]))

        visited_edges = []
        total_cost = 0
        select_edges = grb.tuplelist((i) for i in x_edge.keys() if x_edge[i].X > 0.5)
        for i in select_edges:
            # print(i, t_cost[i])
            total_cost += t_cost[i]
            # visited_edges.append(i)

        """ MELHORAR -> organizar rota """
        if any(select_edges):
            len_route = len(select_edges)
            visited_edges.append(select_edges[0])
            new_find = visited_edges[-1][0]
            if new_find == 0:
                new_find = visited_edges[-1][1]
            select_edges.pop(0)
            while len(visited_edges) != len_route:
                for i in select_edges:
                    if i[0] == new_find:
                        visited_edges.append(i)
                        new_find = i[1]
                        select_edges.remove(i)
                    elif i[1] == new_find:
                        visited_edges.append((i[1], i[0]))
                        new_find = i[0]
                        select_edges.remove(i)
        # print('select', select)
        # print('Profit =', opt_model.objVal)  #np.sum([profit[i] for i in select]))
        # print('visited edges =', visited_edges)
        # print('Total cost =', total_cost)

        solution = {"route": visited_edges,
                    "subgroups_visited": select,
                    # "sets_visited": [],
                    "distance": total_cost,
                    "profit": self.opt_model.objVal
                    }

        return solution


def all_subsets(ss):
    return chain(*map(lambda x: combinations(ss, x), range(2, len(ss) + 1)))


def menor_aresta(lista1, lista2, matriz_distancias):
    menor_distancia = float('inf')  # Inicializa com infinito
    vertice_u, vertice_v = None, None  # Para armazenar os vértices correspondentes à menor distância

    for u in lista1:
        for v in lista2:
            distancia = matriz_distancias[u][v]  # Acessa a matriz de distâncias
            if distancia < menor_distancia:
                menor_distancia = distancia
                vertice_u, vertice_v = u, v

    return vertice_u, vertice_v, menor_distancia




