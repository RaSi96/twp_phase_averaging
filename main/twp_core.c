# include "phase_average.h"

/*
 * ------------------- Core functionality subroutines go here ------------------
 */
int dtw_cost_matrix(struct nd_array *s1, struct nd_array *s2, struct nd_array *r)
{
	int N = r->nrow;
	int M = r->ncol;

	struct nd_array min_array;
	if (init_array(&min_array, 3, 3, 1)) {
		free(min_array.data);
		return 1;
	}

	for (int i = 1; i < N; i++) {
		for (int j = 1; j < M; j++) {
			min_array.data[0] = r->data[(i - 1) * M + j];
			min_array.data[1] = r->data[i * M + (j - 1)];
			min_array.data[2] = r->data[(i - 1) * M + (j - 1)];

			r->data[i * M + j] = euclidean(s1->data[i], s2->data[j]) +
			minimum(false, &min_array);
		}
	}

	free(min_array.data);
	return 0;
}


int dtw_path(struct nd_array *r)
{
	/* No reason to crack your head too much on this: the original cost matrix
	 * can be easily visualised as one axis being arr_1, the other axis being
	 * arr_2. This pathway matrix can actually be visualised the same way: one
	 * axis consists of 2 columns, one each for x- and y-coordinates; the other
	 * axis is simply the total number of elements we have since each element
	 * can be described by an (x, y) coordinate. Therefore:
	 * - 2 columns * the total number of elements = size of the entire array;
	 * - nrows = total number of elements;
	 * - ncols = 2.
	 */
	struct nd_array pathway;
	if (init_array(&pathway, (r->size) * 2, r->size, 2)) {
		fprintf(stderr, "ERROR: dtw_path(): error initialising pathway\n");
		return 1;
	} else {
		for (int i = 0; i < pathway.size; i++) {
			pathway.data[i] = INFINITY;
		}
	}

	struct nd_array min_array;
	if (init_array(&min_array, 3, 3, 1)) {
		fprintf(stderr, "ERROR: dtw_path(): error initialising min_array\n");
		free(pathway.data);
		return 1;
	}

	int ele_count = 1;		//remember we always 1-index our nd_arrays
	int N = r->nrow - 1;	//row coordinate (y-axis)
	int M = r->ncol - 1;	//col coordinate (x-axis)
	int i = 0;
	while ((N > 0) || (M > 0)) {
		pathway.data[i * pathway.ncol + 0] = N;
		pathway.data[i * pathway.ncol + 1] = M;

		min_array.data[0] = r->data[(N - 1) * r->ncol + M];
		min_array.data[1] = r->data[(N - 1) * r->ncol + (M - 1)];
		min_array.data[2] = r->data[N * r->ncol + (M - 1)];

		int mov = minimum(true, &min_array);
		switch (mov) {
			case 0:	//move one row up
				N -= 1;
				break;
			case 1:	//move diagonally
				N -= 1;
				M -= 1;
				break;
			case 2:	//move one col left
				M -= 1;
				break;
		}
		i++;
		ele_count++;
	}

	pathway.data[i * pathway.ncol + 0] = 0;
	pathway.data[i * pathway.ncol + 1] = 0;

	/* Need to finally reverse the pathway array */
	r->size = ele_count * 2;
	r->nrow = ele_count;
	r->ncol = 2;
	r->data = realloc(r->data, r->size);
	for (int i = 0; i < r->size; i += 2) {
		/* this indexing jiggery-pokery is necessary since we're dealing with
		 * (x, y) coordinate pairs. We can't directly reverse-index pathway.data
		 * by doing [(r->size - 1) - i] since we'd instead end up with (y, x)
		 * coordinate pairs. It's endianness all over again! :D
		 */
		r->data[i + 0] = pathway.data[(r->size - 1) - (i + 1)];
		r->data[i + 1] = pathway.data[(r->size - 1) - (i + 0)];
	}

	free(pathway.data);
	free(min_array.data);
	return 0;
}


int difference_path(struct nd_array *r)
{
	struct nd_array diff_path;
	if (init_array(&diff_path, r->size, r->nrow, r->ncol)) {
		fprintf(
			stderr,
		  "ERROR: difference_path(): error initialising diff_path\n"
		);
		return 1;
	}

	int ele_count = 1;
	for (int i = 2; i < r->size; i += 2) {
		diff_path.data[i - 2] = r->data[i + 0] - r->data[i - 2];
		diff_path.data[i - 1] = r->data[i + 1] - r->data[i - 1];
		ele_count++;
	}

	/* As usual, need to modify the original array in-place */
	r->size = (ele_count - 1) * 2;
	r->nrow = (ele_count - 1);
	r->data = realloc(r->data, r->size);
	memcpy(r->data, diff_path.data, r->size * sizeof(float));

	free(diff_path.data);
	return 0;
}


int warp_profile(size_t s1_size, size_t s2_size, struct nd_array *r)
{
	size_t profile_size = (s1_size + s2_size) - 1;	//it's deterministic

	struct nd_array warp_path;
	if (init_array(&warp_path, profile_size, profile_size, 1)) {
		fprintf(stderr, "ERROR: warp_profile(): error initialising warp_path\n");
		return 1;
	}

	int j = 0;		//need separate indexers for warp_path.data and diff_path
	int tmp = 0;
	warp_path.data[0] = 0;

	for (int i = 0; i < r->size; i += 2) {
		tmp = r->data[i + 1] - r->data[i + 0];
		switch(tmp) {
			case -1:
				warp_path.data[j + 1] = warp_path.data[j] - 1;
				j++;
				break;
			case 0:
				warp_path.data[j + 1] = warp_path.data[j];
				warp_path.data[j + 2] = warp_path.data[j];
				j += 2;
				break;
			case 1:
				warp_path.data[j + 1] = warp_path.data[j] + 1;
				j++;
				break;
		}
	}

	/* As usual, need to modify the original array in-place */
	r->size = warp_path.size;
	r->nrow = warp_path.nrow;
	r->ncol = warp_path.ncol;
	r->data = realloc(r->data, r->size);
	memcpy(r->data, warp_path.data, r->size * sizeof(float));

	free(warp_path.data);
	return 0;
}


int twp_average(struct nd_array *s1, struct nd_array *s2,
				struct nd_array *profile)
{
	/* Something to keep in mind about this function: the Python implementation
	 * includes some extra variables and techniques like `wp_half`, `s1_i`,
	 * `s2_i`, `s1_tau` and `s2_tau` because it seeks to explicitly leverage the
	 * vectorised nature of NumPy's functions. Here we can take a few shortcuts
	 * that don't require such intermediary variables because C is a scalar
	 * language. That being said, we can also optimise the Python version to be
	 * just as concise.
	 */
	int s1_i;
	int s2_i;
	size_t half_profile = profile->size / 2;
	float twp_avg[half_profile];

	for (int i = 0, j = 0; i < profile->size; i += 2, j++) {
		s1_i = floor((i - profile->data[i]) / 2);
		s2_i = floor((i + profile->data[i]) / 2);
		twp_avg[j] = (s1->data[s1_i] + s2->data[s2_i]) / 2;
	}

	/* As usual, need to modify the original array in-place */
	profile->size = half_profile;
	profile->nrow = half_profile;
	profile->ncol = 1;
	profile->data = realloc(profile->data, profile->size);
	memcpy(profile->data, twp_avg, profile->size * sizeof(float));

	return 0;
}
