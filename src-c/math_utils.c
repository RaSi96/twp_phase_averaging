# include "phase_average.h"

/*
 * ----------------------- Mathematical functions go here ----------------------
 */
float euclidean(float e1, float e2)
{
	return sqrt(pow((e1 - e2), 2));
}


float minimum(bool argmin, struct nd_array *array)
{
	float min = array->data[0];
	int index = 0;

	for (uint8_t i = 1; i < array->size; i++) {
		float tmp = array->data[i];
		if (tmp < min) {
			min = tmp;
			index = i;
		}
	}

	if (argmin) {
		return index;
	}

	return min;
}
