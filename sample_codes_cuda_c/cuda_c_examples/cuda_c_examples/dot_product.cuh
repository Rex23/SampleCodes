#include <iostream>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "helper_cuda.h"

#ifndef __CUDACC__  
#define __CUDACC__
#endif
#include <device_functions.h>

using namespace std;

static void HandleError(cudaError_t err,
	const char *file,
	int line) {
	if (err != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(err),
			file, line);
		exit(EXIT_FAILURE);
	}
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

int dot_product(void);