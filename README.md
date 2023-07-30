# Introduction
This repository contains basic implementations of the paper written by Sioros, G., & Nymoen, K. (2021, September 2): Accurate shape and phase averaging of time series through Dynamic Time Warping; available on [arXiv.org](https://arxiv.org/abs/2109.00978v1).

The averaging algorithm, as with most time series averaging algorithms, is based off Dynamic Time Warping to boot, but in contrast to other solutions like DBA, PSA and NLAAF, this approach also considers phase alignment through time resulting in a mean sequence that averages phase and amplitude variation.

Note that these algorithms have been shown in their respective papers to significantly outperform a simple point-wise arithmetic mean between two given temporal vectors; of course, these algorithms are also applicable to unequal length vectors.

# Repository Structure
This repository contains 2 implementations, one in C and another in Python, the latter of which is the more usable one (as of the time of this writing). The structure of the repository is as follows:
```
./src-python
	../__init__.py
	../single_avg.py
	../multi_avg.py
	../utils.py
./src-c
	../Makefile
	../main_single.c
	../math_utils.c
	../nd_array.c
	../phase_average.h
	../twp_core.c
```

## Python
### single_avg
This file contains the base implementation of the TWP phase averaging algorithm. Functions within this file can be used to find the average between two sequences, for example:
```python
>>> from single_avg import *
>>> p1 = np.array([1, 2, 4, 9, 6, 5, 6, 9])
>>> p2 = np.array([14, 12, 9, 8, 9, 8, 12])
>>> mat = dtw_cost_matrix(p1, p2)
>>> mat
matrix([[ 0., inf, inf, inf, inf, inf, inf],
	   [inf, 10., 17., 23., 30., 36., 46.],
	   [inf, 18., 15., 19., 24., 28., 36.],
	   [inf, 21., 15., 16., 16., 17., 20.],
	   [inf, 27., 18., 17., 19., 18., 23.],
	   [inf, 34., 22., 20., 21., 21., 25.],
	   [inf, 40., 25., 22., 23., 23., 27.],
	   [inf, 43., 25., 23., 22., 23., 26.]])
>>> pat = dtw_opt_path(mat)
>>> pat
[(0, 0), (1, 1), (2, 2), (3, 3), (3, 4), (4, 5), (5, 5), (6, 5), (7, 6)]
>>> wprof = warp_profile(pat)
>>> wprof
{0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 1, 8: 1, 9: 1, 10: 0, 11: -1, 12: -1, 13: -1}
>>> twp_average(p1, p2, wprof)
array([7.5, 7. , 6.5, 8.5, 9. , 6.5, 7. ])
```

### multi_avg
This file contains an additional two functions that can be used to find the average between more than two sequences. For example:
```python
>>> from multi_avg import *
>>> p1 = np.array([1, 2, 4, 9, 6, 5, 6, 9])
>>> p2 = np.array([14, 12, 9, 8, 9, 8, 12])
>>> p3 = np.array([7, 10, 9, 8, 11, 12])
>>> twp_multi_average([p1, p2, p3])
array([7.33333333, 8.        , 7.33333333, 8.33333333, 8.5       ,
	  8.        ])
```

### utils
This file contains two functions, one for inspecting the optimal DTW path in relation to the original cost matrix; the other for plotting out all involved time sequences and their computed average. For example:
```python
>>> from single_avg import *
>>> from utils import *
>>> p1 = np.array([1, 2, 4, 9, 6, 5, 6, 9])
>>> p2 = np.array([14, 12, 9, 8, 9, 8, 12])
>>> mat = dtw_cost_matrix(p1, p2)
>>> pat = dtw_opt_path(mat)
>>> inspect_path(pat, mat)
matrix(0, 0)    -> 0.0
matrix(1, 1)    -> 10.0
matrix(2, 2)    -> 15.0
matrix(3, 3)    -> 16.0
matrix(3, 4)    -> 16.0
matrix(4, 5)    -> 18.0
matrix(5, 5)    -> 21.0
matrix(6, 5)    -> 23.0
matrix(7, 6)    -> 26.0

>>> from utils import *
>>> p1 = np.array([1, 2, 4, 9, 6, 5, 6, 9])
>>> p2 = np.array([14, 12, 9, 8, 9, 8, 12])
>>> p3 = np.array([7, 10, 9, 8, 11, 12])
>>> plot_averages([p1, p2, p3], ['o', 'o', 'o'], ["p1", "p2", "p3"])
```
![plot_averages() output](triple_avg.png?raw=true)

### Python Caveats
- There is an accompanying Jupyter notebook in HTML showcasing basic Python functionality.
- NumPy is the only hard dependency. The inclusion of `matplotlib.pyplot` is purely for plotting purposes, of which there is only a single function doing very basic plotting. This can be completely ignored if required. Other modules such as `itertools.combinations`, `time`, etc. are modules contained within Python's standard library.
- It's recommended to have time sequences of >= 5 elements at the minimum. Very large sets of sequences (around 1000 sequences and/or >~ 1550-element arrays) will of course take very long to complete (although the algorithm does indeed converge).
- Equal length sequences sent into the averaging function might not yield an appropriate response signal (it might be shifted in phase, or it might line up perfectly) and sometimes might be worse than a standard arithmetic mean. Gathering a guess, this might be due to the nature of the averaging function which samples points twice when compared to the arithmetic mean.

## C
The directory structure here is pretty straightforward, though usability has a few caveats as of this time:
- The Makefile contains two targets, `all` and `clean`, and uses GCC.
- The C version performs all operations in-place. In `main_single.c`, the struct `res_mat` contains the results of each function in the algorithm. Each function itself creates a local structure that it uses, then copies over its contents into `res_mat` before being freed.
- The C version as of now is a lot less usable than its Python counterpart, utilising two hardcoded arrays and printing a lot of verbose output to `stdout` owing to its POC status. The algorithm itself is correct for all tested array pairs, however multi-average interfacing, accepting input from `stdin` or a file, and a few other features are WIP.

Example output:
```sh
$ main >./twp_avg_single
Array 1:
[ 1.000000, ]
[ 2.000000, ]
[ 4.000000, ]
[ 9.000000, ]
[ 6.000000, ]
[ 5.000000, ]
[ 6.000000, ]
[ 9.000000, ]

Array 2:
[ 14.000000, ]
[ 12.000000, ]
[ 9.000000, ]
[ 8.000000, ]
[ 9.000000, ]
[ 8.000000, ]
[ 12.000000, ]

Cost Matrix:
[ 0.0, inf, inf, inf, inf, inf, inf, ]
[ inf, 18., 15., 19., 24., 28., 36., ]
[ inf, 10., 17., 23., 30., 36., 46., ]
[ inf, 21., 15., 16., 16., 17., 20., ]
[ inf, 27., 18., 17., 19., 18., 23., ]
[ inf, 34., 22., 20., 21., 21., 25., ]
[ inf, 40., 25., 22., 23., 23., 27., ]
[ inf, 43., 25., 23., 22., 23., 26., ]

Optimal pathway:
[ 0.0, 0.0, ]
[ 1.0, 1.0, ]
[ 2.0, 2.0, ]
[ 3.0, 3.0, ]
[ 3.0, 4.0, ]
[ 4.0, 5.0, ]
[ 5.0, 5.0, ]
[ 6.0, 5.0, ]
[ 7.0, 6.0, ]

Differenced path:
[ 1.0, 1.0, ]
[ 1.0, 1.0, ]
[ 1.0, 1.0, ]
[ 0.0, 1.0, ]
[ 1.0, 1.0, ]
[ 1.0, 0.0, ]
[ 1.0, 0.0, ]
[ 1.0, 1.0, ]

Time Warp Profile:
[ 0.0, ]
[ 0.0, ]
[ 0.0, ]
[ 0.0, ]
[ 0.0, ]
[ 0.0, ]
[ 0.0, ]
[ 1.0, ]
[ 1.0, ]
[ 1.0, ]
[ 0.0, ]
[ -1.0, ]
[ -1.0, ]
[ -1.0, ]

Phase averaged sequence:
[ 7.5, ]
[ 7.0, ]
[ 6.5, ]
[ 8.5, ]
[ 9.0, ]
[ 6.5, ]
[ 7.0, ]
```
