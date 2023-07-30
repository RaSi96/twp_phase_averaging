# include "phase_average.h"

/*
 * -------------------- Array interfacing functions go here --------------------
 */
uint8_t init_array(struct nd_array *arr, size_t size, size_t nrow, size_t ncol)
{
	arr->size = size;
	arr->nrow = nrow;
	arr->ncol = ncol;
	arr->data = calloc(arr->size, sizeof(float));

	if (!arr->data) {
		fprintf(
			stderr,
		  "ERROR: init_array(): could not allocate %lu bytes for array (%s)",
				arr->size, strerror(errno)
		);
		return 1;
	}

	return 0;
}


void print_ndarray(struct nd_array *arr, char *flavour)
{
	int N = arr->nrow;
	int M = arr->ncol;
	flavour? printf("%s", flavour) : 0;

	for (int i = 0; i < N; i++) {
		printf("[ ");
		for (int j = 0; j < M; j++) {
			printf("%f, ", arr->data[i * M + j]);
		}
		printf("]\n");
	}
	printf("\n");
}
