#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include "devicequery.cuh"
#include "kernel.cuh"
#include "dot_product.cuh"
#include "const_memory.cuh"

using namespace std;

int main(int argc, char** argv)
{
	//Sample code 1: Query Device Properties
	devicequery(argc, argv);
	
	//Sample code 2: Add two arrays
	//add_two_arrays_example();
	
	//Sample code 3: Dot product of two vectors. The shared memory is used for each block.
	//dot_product();

	//Sample code 4: Constant memory usage
	const_memory();

	system("pause");

	// finish
	exit(EXIT_SUCCESS);

	return 1;
}

