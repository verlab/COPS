import argparse
import os
import time

import libs.COPS_ILP as ilpCOPS
from libs.grafo import COPS


def receive_data():
    # Creating an ArgumentParser object
    parser = argparse.ArgumentParser(description='execution arguments')

    # Adding named arguments
    parser.add_argument('--path', help='address of the file you want to run')
    parser.add_argument('--save_img', help='True or False if you need to save the image result or not, DEFAULT=True')
    # Add more named arguments as needed

    # Parsing the arguments received
    args = parser.parse_args()
    return args


def spreadsheet_header(results_file):
    # write the spreadsheet header
    if os.path.exists(results_file):
        with open(results_file, 'r') as file:
            # read the file
            content = file.read()
        # check if the content is empty
        if not content:
            with open(results_file, 'a') as file:
                file.write(
                    f"name_problem;num_vertices;num_subgroups;num_clusters;time;profit;distance;route;subgroups_visited;method")
    else:
        with open(results_file, 'a') as file:
            file.write(
                f"name_problem;num_vertices;num_subgroups;num_clusters;time;profit;distance;route;subgroups_visited;method")


def main(dir, problem, results_file, save_img):
    f = open(results_file, "a")
    """ Read the dataset """
    dataset = fr"{dir}\{problem}.cops"
    cops = COPS()
    cops.read_data(dataset)

    print("------ILP-------")
    t1 = time.time()
    solution = ilpCOPS.main(cops)

    print("------ Final Solution -------")
    print("solution ILP", solution)

    tempoExec = time.time() - t1
    print("Runtime: {} seconds".format(tempoExec))
    print(f"Runtime: {time.strftime('%H:%M:%S', time.gmtime(tempoExec))}")

    legend = ["r"]
    img_path = fr"{dir}/img"
    if not os.path.exists(img_path):
        # Create a new directory because it does not exist
        os.makedirs(img_path)
    img_saved = fr"{img_path}/{problem}"
    cops.draw_2D(path=solution["route"], legend=legend, fill_cluster=True, fill_set=True, name=img_saved,
                 save_img=save_img)

    f.write(
        f"\n{problem};"
        f"{len(cops.list_vertex) - 1};"
        f"{cops.n_subgroups - 1};"
        f"{len(cops.list_clusters) - 1};"
        f"{str(round(tempoExec, 3)).replace('.', ',')};"
        f"{str(solution['profit']).replace('.', ',')};"
        f"{str(round(solution['distance'], 2)).replace('.', ',')};"
        f"{str(solution['route']).replace(',', ' ')};"
        f"{str(solution['subgroups_visited']).replace(',', ' ')};"
        f"COPS-TABU")
    f.close()


if __name__ == '__main__':
    # Getting parsed problem
    args = receive_data()

    # problem file
    dir = fr"{os.getcwd()}\datasets"
    problem = "example_2"
    if args.path:
        dir = os.path.dirname(args.path)
        problem = os.path.basename(args.path).split('.')[0]
    results_path = fr"{dir}\results"

    # img save or not
    save_img = True  # default true
    if args.save_img == 'False':
        save_img = False

    # result file
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    results_file = fr"{results_path}\{problem}.csv"

    # write spreadsheet header
    spreadsheet_header(results_file)

    # call main
    main(dir, problem, results_file, save_img)
    print("THE PROCESS IS FINISHED")

