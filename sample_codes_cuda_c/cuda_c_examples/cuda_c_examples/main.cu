#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include "devicequery.cuh"
#include "kernel.cuh"
#include "dot_product.cuh"
#include "const_memory.cuh"
#include "Event_Create_Record_Synchonize.cuh"

using namespace std;

int main(int argc, char** argv)
{
	int a_case = 4;

	//Sample code: Query Device Properties
	devicequery(argc, argv);

	switch (a_case) {
	case 1:
		//Sample code: Add two arrays
		add_two_arrays_example();
		break;
	case 2:
		//Sample code: Dot product of two vectors. The shared memory is used for each block.
		dot_product();
		break;
	case 3:
		//Sample code: Constant memory usage
		const_memory();
		break;
	case 4:
		//Sample code: Cuda events
		Event_Create_Record_Synchonize();
		break;
	default:
		break;
	}

	system("pause");

	// finish
	exit(EXIT_SUCCESS);

	return 1;
}

