# ifndef PHASE_AVERAGE_H_
# define PHASE_AVERAGE_H_

/* DEFINITIONS AND INCLUSIONS */
# include <errno.h>
# include <math.h>
# include <stdbool.h>
# include <stdint.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

/* STRUCTS */
struct nd_array {
	size_t size;	//total number of elements in *data, not number of bytes!
	size_t nrow;
	size_t ncol;
	float *data;
};

/* USERDEF FUNCTION PROTOTYPES */
//nd_array functions
uint8_t init_array(struct nd_array *arr, size_t size, size_t nrow, size_t ncol);
void print_ndarray(struct nd_array *arr, char *flavour);

//math utils
float euclidean(float e1, float e2);
float minimum(bool argmin, struct nd_array *array);

//twp core
int dtw_cost_matrix(struct nd_array *s1, struct nd_array *s2,
					struct nd_array *r);
int dtw_path(struct nd_array *r);
int difference_path(struct nd_array *r);
int warp_profile(size_t s1_size, size_t s2_size, struct nd_array *r);
int twp_average(struct nd_array *s1, struct nd_array *s2,
				struct nd_array *profile);

# endif /* PHASE_AVERAGE_H_ */
