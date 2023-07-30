# general imports
import matplotlib.pyplot as plt

# specific imports
from multi_avg import *

# ------------------------------------------------------------------------------

def inspect_path(path, cost_mat):
    """
    Function that helps in charting out what the optimal path looks like within
    the context of the corresponding cost matrix. Usually used for verification
    and/or introspection.

    Parameters:
    `path`: iterable of tuples.
        The optimal path through the DTW cost matrix. Any iterable will do as
        long as the tuples of indices are valid coordinates in `cost_mat`.

    `cost_mat`: np.matrix
        The DTW cost matrix for the corresponding sequences under analysis.

    Returns nothing, prints the corresponding cost for each tuple of indices in
    `path`.

    Example:
    >>> inspect_path(path, mat)
    matrix(0, 0)	-> 0.0
    matrix(1, 1)	-> 10.0
    matrix(2, 2)	-> 15.0
    matrix(3, 3)	-> 16.0
    matrix(3, 4)	-> 16.0
    matrix(4, 5)	-> 18.0
    matrix(5, 5)	-> 21.0
    matrix(6, 5)	-> 23.0
    matrix(7, 6)	-> 26.0
    """
    for p in path:
        print(f"matrix({p[0]}, {p[1]})\t-> {cost_mat[p[0], p[1]]}")

    return


def plot_averages(seq_list, marker_list = None, label_list = None):
    """
    Tiny function that plots out any given original sequences along with their
    TWP phase average for comparison.

    Parameters:
    `seq_list`: iterable of time sequences
        A collection of time sequences to compute distances between. Any kind of
        iterable will suffice, as long as there is no datetime column, just the
        sequence of values.

    `marker_list`: iterable of characters
        A collection of valid matplotlib marker styles to use for the plot, if
        desired. Refer to the matplotlib documentation for valid characters.

    `label_list`: iterable of strings
        A collection of strings to use for each sequence's plot as its label.
        These labels show up in the plot's legend.

    Returns nothing, plots internally.

    Example:
    >>> plot_averages(
            seq_list = [p1, p1, p2],
            marker_list = ['o', 'o', 'o'],
            label_list = ["p1", "p2", "p3"]
        )
    """
    # MVCE:
    # plt.plot(p1, marker = 'o', label = "p1")
    # plt.plot(p2, marker = 'o', label = "p2")
    # plt.plot(p3, marker = 'o', label = "p3")
    # plt.plot(twp_multi_average([p1, p2, p3]), marker = 'o', label = "avg")
    # plt.grid(True, linestyle = "--")
    # plt.legend();

    # locally-scoped variables
    plt_plot = plt.plot

    # functionality
    for s, m, l in zip(seq_list, marker_list, label_list):
        plt_plot(s, marker = m, label = l)

    plt_plot(twp_multi_average(seq_list), marker = 'o', label = "avg")
    plt.grid(True, linestyle = "--")
    plt.legend()
    plt.show()
    return
