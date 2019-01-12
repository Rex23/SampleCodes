#include <iostream>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "helper_cuda.h"

#ifndef __CUDACC__  
#define __CUDACC__
#endif
#include <device_functions.h>
#include "cublas_v2.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M 6
#define N 5
#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))

using namespace std;

int cblas_one_based_indexing();