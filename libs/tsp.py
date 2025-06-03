""" Traveling salesman problem with 2-opt """
""" Thanks: https://www.keiruaprod.fr/blog/2021/09/15/traveling-salesman-with-2-opt.html
    and https://stackoverflow.com/questions/25585401/travelling-salesman-in-scipy """

# tsp.py
import numpy as np

# Calculate the euclidian distance in n-space of the route r traversing cities c, ending at the path start.
path_distance_for_circular = lambda r, c: np.sum([np.linalg.norm(c[r[p]]-c[r[p-1]]) for p in range(len(r))])  # for circular path
path_distance_non_circular = lambda r, c: np.sum([np.linalg.norm(c[r[p+1]]-c[r[p]]) for p in range(len(r)-1)])  # for non-circular path
# Reverse the order of all elements from element i to element k in array r.
two_opt_swap = lambda r, i, k: np.concatenate((r[0:i],r[k:-len(r)+i-1:-1],r[k+1:len(r)]))


def two_opt(cities, improvement_threshold, is_a_circular_path=True):  # 2-opt Algorithm adapted from https://en.wikipedia.org/wiki/2-opt
    route = np.arange(cities.shape[0])  # Make an array of row numbers corresponding to cities.
    improvement_factor = 1  # Initialize the improvement factor.

    if is_a_circular_path:
        aux = 0  # first point will be the start and the end
        path_distance = path_distance_for_circular
    else:
        aux = 1  # first point will be the start and the last point is the end
        path_distance = path_distance_non_circular
    best_distance = path_distance(route, cities)  # Calculate the distance of the initial path.
    while improvement_factor > improvement_threshold:  # If the route is still improving, keep going!
        distance_to_beat = best_distance  # Record the distance at the beginning of the loop.
        for swap_first in range(1, len(route)-(2+aux)):  # From each city except the first and last,
            for swap_last in range(swap_first+1, len(route)-aux):  # to each of the cities following,
                new_route = two_opt_swap(route, swap_first, swap_last)  # try reversing the order of these cities
                new_distance = path_distance(new_route, cities)  # and check the total distance with this modification.
                if new_distance < best_distance:  # If the path distance is an improvement,
                    route = new_route  # make this the accepted best route
                    best_distance = new_distance  # and update the distance corresponding to this route.
        if distance_to_beat == 0:
            improvement_factor = improvement_threshold
        else:
            improvement_factor = 1 - best_distance / distance_to_beat  # Calculate how much the route has improved.
            #print("improvement_factor", improvement_factor)
    return route  # When the route is no longer improving substantially, stop searching and return the route.