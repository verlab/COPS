""" Thanks
https://www.gurobi.com/features/academic-named-user-license/    # install gurobi
https://support.gurobi.com/hc/en-us/community/posts/4404976974865-Orienteering-Problem
"""

# imports
import random
import gurobipy as grb
import numpy as np
from itertools import chain, combinations


def all_subsets(ss):
    return chain(*map(lambda x: combinations(ss, x), range(2, len(ss) + 1)))


def main(csop_class):
    matrix_dist = csop_class.matrix_dist
    circular_path = csop_class.circular_path
    max_distance = csop_class.t_max
    vertex = np.array([np.array(x) for x in csop_class.list_vertex], dtype=object)
    cluster = np.array([np.array(x) for x in csop_class.list_subgroups], dtype=object)
    profit = np.array([np.array(x) for x in csop_class.profit], dtype=object)
    groups = np.array([np.array(x) for x in csop_class.list_clusters], dtype=object)
    end_set = csop_class.end_cluster

    opt_model = grb.Model()  # "MILP Model")#name="Visit model"
    opt_model.setParam('Threads', 4)  # Define o número de threads desejado, neste caso, 4.

    """MELHORAR"""
    possible_end_start_vertex = [v for c in groups[end_set] for v in cluster[c]]  # cluster[groups[end_cluster]]
    possible_end_start_vertex.append(0)
    print("possible_end_start_vertex", possible_end_start_vertex)
    clientes = []
    for i in range(1, len(vertex)):
        um_cliente = True
        for j in possible_end_start_vertex:
            if i == j:
                um_cliente = False
                break
        if um_cliente:
            clientes.append(i)
        # print("Estes são os clientes", clientes)

    # vertex = [(0, 0), (1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (9, 9), (100, 100)]
    qtd_v = len(vertex)
    qtd_set = len(cluster)

    # groups = [(1, 3), (2, 4)]
    g = {}
    for i in range(len(groups)):
        g[i] = opt_model.addVar(vtype=grb.GRB.BINARY, name='g' + str(i))

    # 11
    """ zi ∈ {0, 1}, ∀Ci ∈ C """
    z = {}
    for i in range(qtd_set):
        z[i] = opt_model.addVar(vtype=grb.GRB.BINARY, name='z' + str(i))

    # 5, 12
    """ ∑ texe ≤ Tmax """
    """ xe ∈ {0, 1}, ∀e ∈ E """
    t_cost = {}
    x_edge = {}
    for i in range(qtd_v):
        for j in range(i + 1, qtd_v):
            # if i != j:
            x_edge[i, j] = opt_model.addVar(vtype=grb.GRB.BINARY, name='edge_' + str(i) + '_' + str(j))
            t_cost[i, j] = matrix_dist[i][j]  # distance(vertex, i, j)
            # print("Todos", i, j)

    # 13
    """ yj ∈ {0, 1}, ∀vj ∈ V """
    y = {}
    for i in range(qtd_v):
        y[i] = opt_model.addVar(vtype=grb.GRB.BINARY, name='y' + str(i))

    # constraints
    # 2
    """ y0 = 1 """
    opt_model.addConstr(z[0] == 1)  # (2)
    opt_model.addConstr(y[0] == 1)  # (2)

    # 3
    """ ∑ xe = 2yj ∀vj ∈ S\{S0, Send} """
    aux1 = np.array(range(qtd_v))
    # for j in range(qtd_v):  #(3)
    for j in clientes:
        aux = aux1[aux1 != j]
        opt_model.addConstr(grb.quicksum((x_edge[j, k] if k > j else x_edge[k, j]) for k in aux) == 2 * y[j])

    # 4 start restriction
    """ ∑ xe(s0) ≤ 2"""
    start_cluster = 0
    aux = aux1[aux1 != start_cluster]
    opt_model.addConstr(
        grb.quicksum((x_edge[start_cluster, k] if k > start_cluster else x_edge[k, start_cluster]) for k in aux) <= 2)

    # 5
    """ ∑ texe ≤ Tmax """
    # for i in range(qtd_v):  #(5)
    opt_model.addConstr(
        grb.quicksum(t_cost[i, j] * x_edge[i, j] for i in range(qtd_v) for j in range(i + 1, qtd_v)) <= max_distance)
    # aux = aux1[aux1 != i]
    # opt_model.addConstr(grb.quicksum(((t_cost[i, j]*x_edge[i, j]) if j > i else (t_cost[j, i]*x_edge[j, i])) for j in aux) <= max_distance)

    # 6
    """ sub-tour constraint """
    """for subset in all_subsets(range(1, qtd_v)):  # (6)
        qtd = len(subset)
        # print('subset', subset, subset[0])
        for t in subset:
            opt_model.addConstr(
                grb.quicksum(x_edge[subset[i], subset[j]] for i in range(qtd) for j in range(i + 1, qtd)) <=
                grb.quicksum(y[k] for k in subset if k != t))"""
    for subset in all_subsets(clientes):  # (6)
        qtd = len(subset)
        subset_barr = np.delete(aux1, subset)
        #print('subset', subset, type(subset))
        #print('subset_barr', subset_barr)
        #print(qtd)
        #for t in subset:
        opt_model.addConstr(
            grb.quicksum((x_edge[i, j]) if i < j else (x_edge[j, i])
                         for i in subset for j in subset_barr) >= grb.quicksum(y[k] for k in subset) / qtd)
        #opt_model.addConstr(
        #    (grb.quicksum((x_edge[i, j]) if i < j else (x_edge[j, i]) for i in subset for j in subset_barr) > 0)
        #    if grb.quicksum(y[k] for k in subset) > 0 else
        #    (grb.quicksum((x_edge[i, j]) if i < j else (x_edge[j, i]) for i in subset for j in subset_barr) == 0))
    # 7
    """ zi ≤ yj , ∀Ci ∈ C, ∀vj ∈ Ci """
    for i in range(1, qtd_set):  # (7)
        # opt_model.addConstr(y[cluster[i][j]] for j in range(len(cluster[i])) >= z[i])
        for j in range(len(cluster[i])):
            # print('aaaaa', i, j, len(cluster[i]))
            opt_model.addConstr(y[cluster[i][j]] >= z[i])

    # 8
    """ ∑zi ≤ sg ∀Sg ∈ S """
    for i in range(len(groups)):  # new constraint
        ''' não pode atender a dois grupos COMPLETOS ao mesmo tempo 
            mas pode atender a um grupo completo e a alguns indivíduos do outro
        '''
        opt_model.addConstr(grb.quicksum(z[groups[i][j]] for j in range(len(groups[i]))) <= g[i])

    # 9
    """ ∑yi ≤ ∑Zi | Ci | """
    ''' se atende a um grupo não pode atender a nenhum indivíduo do outro '''
    opt_model.addConstr(grb.quicksum(y[v] for v in range(len(vertex))) <= grb.quicksum(
        z[c] * len(cluster[c]) for c in range(len(cluster))))

    # 14, 15
    if not circular_path:
        """ novo 1 vertice do último superset atendido"""
        opt_model.addConstr(grb.quicksum(y[v] for c in groups[end_set] for v in cluster[c]) == 1)

        """  """
        opt_model.addConstr(grb.quicksum(
            x_edge[min(vvv, end), max(vvv, end)] for vvv in range(len(vertex)) for c in groups[end_set] for end in
            cluster[c] if vvv != end) == 1)

    # 1
    """ objective max profit """
    objective = grb.quicksum(profit[i] * z[i] for i in range(qtd_set))
    opt_model.ModelSense = grb.GRB.MAXIMIZE
    opt_model.setObjective(objective)
    opt_model.optimize()

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
        solution = opt_model.getAttr('x', z)
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

    """ MELHORAR -> organizar rota"""
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
                "clusters_visited": select,
                # "sets_visited": [],
                "distance": total_cost,
                "profit": opt_model.objVal
                }

    return solution