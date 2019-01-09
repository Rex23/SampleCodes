#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);
int add_two_arrays_example();
__global__ void addKernel(int *c, const int *a, const int *b);
