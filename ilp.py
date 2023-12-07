import numpy as np
import os
import time
import sys

import libs.COPS_ILP as cops
from libs.grafo import COPS



dir = r"datasets/comparacao_ilp_cops/douglas/ilp_teste"
#v = f"{dir}/results_experimento_{num_exp}.csv"
v = f"{dir}/results_experimento_ILP.csv"
f = open(v, "a")
# Lê o conteúdo do file
with open(v, 'r') as arquivo:
    # Ler o conteúdo do file
    conteudo = arquivo.read()
# Verifica se o conteúdo está vazio
if not conteudo:
    f.write(f"name_problem;num_vertices;num_subgroups;num_clusters;time (s);profit;distance;RAM(GB)")  # ;route;clusters_visited")

f.close()


def main(dir, problem, replication=1):
    csop = COPS()

    dataset = fr"{os.getcwd()}/{dir}/{problem}.csop"
    #v = fr"{dir}/results/{replication}_{problem}.csv"
    #f = open(v, "a")
    #v = fr"{dir}/results_ilp/{problem}_ilp.csv"
    f = open(v, "a")
    #f.write(f"name_problem,num_vertices,num_subgroups,num_clusters,time,profit,distance,route,clusters_visited")

    if replication == 1:
        save_img = True
    else:
        save_img = False

    """ Read the dataset """
    csop.read_data(dataset)

    print("------ILP-------")
    t1 = time.time()

    solution = cops.main(csop)

    print("------ Final Solution -------")
    print("solution ILP", solution)

    tempoExec = time.time() - t1
    print("Tempo de execução: {} segundos".format(tempoExec))
    print("Tempo de execução:", time.strftime('%H:%M:%S', time.gmtime(tempoExec)))

    legend = ["r"]
    #csop.draw_2D(path=solution["route"], legend=legend, fill_cluster=True, fill_set=True, name=f"{dir}/img/{problem}_ILP", save_img=save_img)

    f.write(
        f"\n{problem};{len(csop.list_vertex) - 1};{csop.n_subgroups - 1};{len(csop.list_clusters) - 1};{str(round(tempoExec, 3)).replace('.', ',')};{str(solution['profit']).replace('.', ',')}"
        f";{str(round(solution['distance'], 2)).replace('.', ',')};{str(round((peak_ram_usage / (1024 * 1024 * 1024)), 5))}")  # ;{str(solution['route']).replace(',', ' ')};{str(solution['clusters_visited']).replace(',', ' ')}")
    f.close()


if __name__ == '__main__':
    replications = 1
    for nome_arquivo in os.listdir(dir):
        if nome_arquivo.endswith('.csop'):
            experiment = os.path.splitext(nome_arquivo)[0]
            print(experiment)
            for r in range(1, replications + 1):
                main(dir, experiment, r)
                peak_ram_usage = 0
                time.sleep(5)
            print(f'-----------experiment {experiment} has finished -----------------')

    print("TERMINOU O PROCESSO")

