# specific imports
from single_avg import *

# ------------------------------------------------------------------------------

# Multi average interfacing

def unequal_euc(*seqs):
    """
    Slight modification to the existing Euclidean distance function that allows
    for unequal length sequences, as described in equation 11.

    Parameters:
    `seqs`: iterable of np.array
        A collection of time sequences to compute distances between. Any kind of
        iterable will suffice, as long the time sequences are NumPy arrays.
        There should not be any datetime column, just a sequence of values.

    Returns a scalar value (the squared Euclidean distance).

    Example:
    >>> unequal_euc(p1, p2)
    matrix([[  0.        , 357.71428571],
            [  0.        ,   0.        ]])

    >>> p3 = np.array([7, 10, 9, 8, 11, 12])
    >>> unequal_euc(p1, p2, p3)
    matrix([[  0.        , 347.42857143, 201.33333333],
            [  0.        ,   0.        ,  66.5       ],
            [  0.        ,   0.        ,   0.        ]])
    """
    py_len = len
    py_range = range
    np_sum = np.sum

    n = py_len(seqs)
    d_mat = np.matrix(np.zeros((n, n)))
    combos = list(combinations(range(n), 2))
    min_len = min([py_len(v) for v in seqs]) - 1

    for i, j in combos:
        s_1 = seqs[i]
        s_1_len = py_len(s_1)
        s_2 = seqs[j]
        s_2_len = py_len(s_2)

        arr = [(s_1[i] - s_2[i])**2 for i in py_range(min_len)]
        d_mat[i, j] = (s_1_len / s_2_len) * np_sum(arr)

    return d_mat


def twp_multi_average(seq_list, thresh = 0):
    """
    The iterative algorithm that finds the barycentric phase average of all
    provided time_sequences, in `aligned_t_arr`. The algorithm, as specified,
    is as follows:
    ```
    while max(dist) > thresh {
        1. compute the distance matrix of all time series received
        2. select the pair of sequences where dist(p, q) = max(dist)
        3. compute the phase array (TWP between p, q)
        4. compute the phase average between (p, q)
        5. replace (p, q) in the original matrix with their average
        6. recompute the distance matrix
    }
    7. compute the element-wise average between all replaced sequences
    F. return the barycentric phase average computed in step 7.
    ```

    Returns a numpy array containing the barycentric average of all time series
    that were sent in. Note: given in an n-dimensional heterogeneous matrix of
    multiple time sequences (each of which can be sized differently), only a
    vector will be returned.

    Example:
    >>> twp_multi_average([p1, p2])
    array([7.5, 7. , 6.5, 8.5, 9. , 6.5, 7. ])

    >>> twp_multi_average([p1, p2, p3])
    array([7.33333333, 8.        , 7.33333333, 8.33333333, 8.5       ,
           8.        ])
    """
    np_unrav = np.unravel_index
    np_argmax = np.argmax
    np_arr = np.array

    get_euc = unequal_euc
    get_mat = dtw_cost_matrix
    get_pat = dtw_opt_path
    get_prf = warp_profile
    get_avg = twp_average

    d_mat = get_euc(*seq_list)

    while d_mat.max() > thresh:
        i, j = np_unrav(np_argmax(d_mat), d_mat.shape)
        s1 = seq_list[i]
        s2 = seq_list[j]

        # order of operations:
        # dtw_cost_matrix -> pathway -> warp_profile -> twp_average
        avg = get_avg(
            s1,
            s2,
            get_prf(
                get_pat(
                    get_mat(s1, s2)
                )
            )
        )

        seq_list[i] = seq_list[j] = avg
        d_mat = get_euc(*seq_list)

    return np.mean(seq_list, axis = 0)
