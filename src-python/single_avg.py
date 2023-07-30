# general imports
import numpy as np

# specific imports
from itertools import combinations

# ------------------------------------------------------------------------------

# Single average interfacing

def dtw_cost_matrix(seq_1, seq_2, distance = np.linalg.norm):
    """
    Dynamic Time Warping function that returns the cumulative cost matrix given
    two time sequences, `seq_1` and `seq_2`. The sequences are expected to be
    numpy arrays.

    This function seeds a matrix of `seq_1` rows and `seq_2` columns with np.inf
    and sets coordinate (0, 0) to 0, in order to organise the eventual warping
    path in an easier fashion. After initialisation, each cell is filled with
    the values as the sum of two components: the Euclidean distances between
    each element of both series, and the minimum value of the 3 adjacent,
    preceeding elements.

    Parameters:
    `seq_1`: np.array
        A 1D time sequence (a vector) to average with. There should not be any
        datetime column, just a sequence of values.

    `seq_2`: np.array
        A 1D time sequence (a vector) to average with. There should not be any
        datetime column, just a sequence of values.

    Returns a numpy matrix shaped N x M, where N is the length of `seq_1` and
    M is the length of `seq_2`.

    Example:
    >>> p1 = np.array([1, 2, 4, 9, 6, 5, 6, 9])
    >>> p2 = np.array([14, 12, 9, 8, 9, 8, 12])
    >>> dtw_cost_matrix(p1, p2)
    matrix([[ 0., inf, inf, inf, inf, inf, inf],
            [inf, 10., 17., 23., 30., 36., 46.],
            [inf, 18., 15., 19., 24., 28., 36.],
            [inf, 21., 15., 16., 16., 17., 20.],
            [inf, 27., 18., 17., 19., 18., 23.],
            [inf, 34., 22., 20., 21., 21., 25.],
            [inf, 40., 25., 22., 23., 23., 27.],
            [inf, 43., 25., 23., 22., 23., 26.]])
    """
    # locally-scoped variables
    py_len = len
    py_min = min
    py_range = range
    fn_dist = distance

    # functionality
    dtw_matrix = np.matrix(np.ones((py_len(seq_1), py_len(seq_2))) * np.inf)
    dtw_matrix[0, 0] = 0.0

    for i in py_range(1, py_len(seq_1)):
        for j in py_range(1, py_len(seq_2)):
            cost = distance(seq_1[i] - seq_2[j])
            dtw_matrix[i, j] = cost + py_min(
                dtw_matrix[i-1, j],
                dtw_matrix[i, j-1],
                dtw_matrix[i-1, j-1]
            )

    return dtw_matrix


def dtw_opt_path(cost_matrix):
    """
    Function to chart the optimal warping path through a provided DTW matrix,
    `mat`. Starting from the end of the matrix, i.e. from position (N, M), it
    traverses the matrix backwards looking for the minimum value between
    elements adjacent to the current one in forward-step. E.g. for the matrix:
    1 | 2 | 3
    --+---+--
    4 | 5 | 4
    --+---+--
    3 | 2 | 1

    It begins from position (n, m) with the first value being 1. The second
    value is the minimum of (4, 5, 2) = 2. The third value is the minimum of
    (5, 4, 3) = 3. The fourth value is the minimum of (5, 4) = 4, and the final
    value is the minumum of (5, 2, 1) = 1. The path is therefore
    [(3, 3), (3, 2), (3, 1), (2, 1), (1, 1)].

    Parameters:
    `cost_matrix`: np.matrix
        A DTW cost matrix between two time sequences shaped N x M.

    Returns a list of tuples, where each tuple contains the indices of the
    optimal path to take through `cost_matrix`.

    Example:
    >>> mat = dtw_cost_matrix(p1, p2)
    >>> dtw_opt_path(mat)
    [(0, 0), (1, 1), (2, 2), (3, 3), (3, 4), (4, 5), (5, 5), (6, 5), (7, 6)]
    """
    # locally-scoped variables
    matches = []
    m_append = matches.append
    np_argmin = np.argmin

    # functionality
    (i, j) = cost_matrix.shape
    i -= 1
    j -= 1

    while i > 0 or j > 0:
        m_append((i, j))
        mov = np_argmin([
            cost_matrix[i - 1, j],      # move one row up
            cost_matrix[i - 1, j - 1],  # move diagonally
            cost_matrix[i, j - 1]       # move one col left
        ])

        match mov:
            case 0:
                i -= 1
            case 1:
                i -= 1
                j -= 1
            case 2:
                j -= 1

    m_append((0, 0))
    matches.reverse()
    return matches


def warp_profile(dtw_path):
    """
    Function that calculates the warp profile between two time sequences given
    their optimal DTW path. The DTW path will be in terms of (t1, t2) coordinate
    tuples representing the points that ought to be aligned; the TWP path will
    be in terms of (tau, phi) which represents phase difference over time.

    Parameters:
    `dtw_path`: iterable of tuples
        The optimal DTW path. ANy iterable will work, the only constraint is the
        iterable ought to contain tuples of indices.

    Returns a dictionary containg tau: phi key-value pairs, where tau is the
    time index and phi is the phase difference between both time sequences.

    Example:
    >>> path = dtw_opt_path(mat)
    >>> warp_profile(path)
    {0: 0,
     1: 0,
     2: 0,
     3: 0,
     4: 0,
     5: 0,
     6: 0,
     7: 1,
     8: 1,
     9: 1,
     10: 0,
     11: -1,
     12: -1,
     13: -1}
    """
    # locally-scoped variables
    py_len = len
    np_arr = np.array

    # functionality
    diff_path = [
        np_arr(dtw_path[i + 1]) - np_arr(dtw_path[i])
        for i in range(py_len(dtw_path) - 1)
    ]

    phi = {0: 0}
    t = 0

    for d in diff_path:  # DO NOT ENUMERATE dif WITH {t, d in enum(dif)}!!!
        match np.diff(d):
            case 0:
                phi[t + 1] = phi[t + 2] = phi[t]
                t += 2
            case -1:
                phi[t + 1] = phi[t] - 1
                t += 1
            case 1:
                phi[t + 1] = phi[t] + 1
                t += 1

    return phi


def twp_average(seq_1, seq_2, warp_profile):
    """
    Function that computes the time warp profile's phase average between two
    given time sequences. The TWP's {Tau: Phi} coordinates maintain phase
    alignment between the sequences; these are used to find the phase-aligned
    arithmetic mean between corresponding elements of the two sequences.

    Parameters:
    `seq_1`: np.array
        A 1D time sequence (a vector) to average with. There should not be any
        datetime column, just a sequence of values.

    `seq_2`: np.array
        A 1D time sequence (a vector) to average with. There should not be any
        datetime column, just a sequence of values.

    `warp_profile`: dictionary
        The warp profile dictionary, where keys are Tau (time indices) and
        values are Phi (phase alignment).

    Returns a numpy nd-array sized `np.ceil((L_1 + L_2 - 1) / 2)`, where L_1 and
    L_2 are the lengths of `seq_1`, and `seq_2`.

    Example:
    >>> twp_average(p1, p2, wprof)
    array([7.5, 7. , 6.5, 8.5, 9. , 6.5, 7. ])
    """
    # locally-scoped variables
    py_int = int
    np_floor = np.floor
    wp_half = list(warp_profile.items())[::2]

    # functionality
    s1_i = np_floor([(t - p) / 2 for t, p in wp_half])
    s2_i = np_floor([(t + p) / 2 for t, p in wp_half])

    s1_tau = [seq_1[py_int(i)] for i in s1_i]
    s2_tau = [seq_2[py_int(i)] for i in s2_i]

    return np.mean([s1_tau, s2_tau], axis = 0)
