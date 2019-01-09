//Note: glut64.dll is needed. It is put under C:\Windows\System32
//      glut32.lib is needed. Goto Property Pages -> Linker -> General -> Additionl Library Directories and add $(SolutionDir)cuda_c_examples\lib;

#include <iostream>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "helper_cuda.h"
#include "cpu_bitmap.h"

#ifndef __CUDACC__  
#define __CUDACC__
#endif
#include <device_functions.h>

using namespace std;

#ifndef HANDLE_ERROR( err )
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
#endif


int const_memory(void);