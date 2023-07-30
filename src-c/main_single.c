# include "phase_average.h"

/*
 * ------------------------------------ Main -----------------------------------
 */
int main(void)
{
	int err = 0;

	// first time array
	struct nd_array arr_1;
	if (init_array(&arr_1, 8, 8, 1)) {
		fprintf(stderr, "ERROR: main(): error initialising arr_1\n");
		err++;
	}

	float arr_1_contents[] = {1, 2, 4, 9, 6, 5, 6, 9};
	memcpy(arr_1.data, arr_1_contents, sizeof(arr_1_contents));
	print_ndarray(&arr_1, "Array 1:\n");

	// second time array
	struct nd_array arr_2;
	if (init_array(&arr_2, 7, 7, 1)) {
		fprintf(stderr, "ERROR: main(): error initialising arr_2\n");
		err++;
	}

	float arr_2_contents[] = {14, 12, 9, 8, 9, 8, 12};
	memcpy(arr_2.data, arr_2_contents, sizeof(arr_2_contents));
	print_ndarray(&arr_2, "Array 2:\n");

	// cost matrix
	struct nd_array cost_mat;
	if (init_array(&cost_mat, arr_1.size * arr_2.size, arr_1.size, arr_2.size)) {
		fprintf(stderr, "ERROR: main(): error initialising cost_mat\n");
		err++;
	}

	for (uint8_t i = 0; i < cost_mat.size; i++) {
		cost_mat.data[i] = INFINITY;
	}
	cost_mat.data[0] = 0.0;

	if (dtw_cost_matrix(&arr_1, &arr_2, &cost_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating cost matrix\n");
		err++;
	}

	print_ndarray(&cost_mat, "Cost Matrix:\n");
	/* [ 0.0, inf, inf, inf, inf, inf, inf, ]
	 * [ inf, 10., 17., 23., 30., 36., 46., ]
	 * [ inf, 18., 15., 19., 24., 28., 36., ]
	 * [ inf, 21., 15., 16., 16., 17., 20., ]
	 * [ inf, 27., 18., 17., 19., 18., 23., ]
	 * [ inf, 34., 22., 20., 21., 21., 25., ]
	 * [ inf, 40., 25., 22., 23., 23., 27., ]
	 * [ inf, 43., 25., 23., 22., 23., 26., ]
	 */

	if (dtw_path(&cost_mat)) {
		fprintf(stderr, "ERROR: main(): error finding optimal DTW path\n");
		err++;
	}

	print_ndarray(&cost_mat, "Optimal pathway:\n");

	if (difference_path(&cost_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating phase warp path\n");
		err++;
	}

	print_ndarray(&cost_mat, "Differenced path:\n");

	if (warp_profile(arr_1.size, arr_2.size, &cost_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating warp profile\n");
		err++;
	}

	print_ndarray(&cost_mat, "Time Warp Profile:\n");

	if (twp_average(&arr_1, &arr_2, &cost_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating TWP average\n");
		err++;
	}

	print_ndarray(&cost_mat, "Phase averaged sequence:\n");

	free(arr_1.data);
	free(arr_2.data);
	free(cost_mat.data);
	if (err > 0) {
		return 1;
	}
	return 0;
}
