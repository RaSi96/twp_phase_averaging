# include "phase_average.h"

/*
 * ------------------------------------ Main -----------------------------------
 */
int main(void)
{
	// first time array
	float arr_1_contents[] = {1, 2, 4, 9, 6, 5, 6, 9};
	size_t arr_1_size = sizeof(arr_1_contents) / sizeof(arr_1_contents[0]);

	struct nd_array arr_1;
	if (init_array(&arr_1, arr_1_size, arr_1_size, 1)) {
		fprintf(stderr, "ERROR: main(): error initialising arr_1\n");
		return 1;
	}

	memcpy(arr_1.data, arr_1_contents, sizeof(arr_1_contents));
	print_ndarray(&arr_1, "Array 1:\n");

	// second time array
	float arr_2_contents[] = {14, 12, 9, 8, 9, 8, 12};
	size_t arr_2_size = sizeof(arr_2_contents) / sizeof(arr_2_contents[0]);

	struct nd_array arr_2;
	if (init_array(&arr_2, arr_2_size, arr_2_size, 1)) {
		fprintf(stderr, "ERROR: main(): error initialising arr_2\n");
		free(arr_1.data);
		return 1;
	}

	memcpy(arr_2.data, arr_2_contents, sizeof(arr_2_contents));
	print_ndarray(&arr_2, "Array 2:\n");

	// results matrix
	struct nd_array res_mat;
	if (init_array(&res_mat, arr_1.size * arr_2.size, arr_1.size, arr_2.size)) {
		fprintf(stderr, "ERROR: main(): error initialising res_mat\n");
		free(arr_1.data);
		free(arr_2.data);
		return 1;
	}

	for (uint8_t i = 0; i < res_mat.size; i++) {
		res_mat.data[i] = INFINITY;
	}
	res_mat.data[0] = 0.0;

	if (dtw_cost_matrix(&arr_1, &arr_2, &res_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating cost matrix\n");
		free(arr_1.data);
		free(arr_2.data);
		free(res_mat.data);
		return 1;
	}

	print_ndarray(&res_mat, "Cost Matrix:\n");

	if (dtw_path(&res_mat)) {
		fprintf(stderr, "ERROR: main(): error finding optimal DTW path\n");
		free(arr_1.data);
		free(arr_2.data);
		free(res_mat.data);
		return 1;
	}

	print_ndarray(&res_mat, "Optimal pathway:\n");

	if (difference_path(&res_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating phase warp path\n");
		free(arr_1.data);
		free(arr_2.data);
		free(res_mat.data);
		return 1;
	}

	print_ndarray(&res_mat, "Differenced path:\n");

	if (warp_profile(arr_1.size, arr_2.size, &res_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating warp profile\n");
		free(arr_1.data);
		free(arr_2.data);
		free(res_mat.data);
		return 1;
	}

	print_ndarray(&res_mat, "Time Warp Profile:\n");

	if (twp_average(&arr_1, &arr_2, &res_mat)) {
		fprintf(stderr, "ERROR: main(): error calculating TWP average\n");
		free(arr_1.data);
		free(arr_2.data);
		free(res_mat.data);
		return 1;
	}

	print_ndarray(&res_mat, "Phase averaged sequence:\n");

	free(arr_1.data);
	free(arr_2.data);
	free(res_mat.data);
	return 0;
}
