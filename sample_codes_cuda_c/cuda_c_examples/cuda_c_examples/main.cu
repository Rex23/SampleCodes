#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include "devicequery.cuh"
#include "kernel.cuh"
#include "dot_product.cuh"

using namespace std;

int main(int argc, char** argv)
{
	//devicequery(argc, argv);
	//add_two_arrays_example();
	dot_product();

	system("pause");

	// finish
	exit(EXIT_SUCCESS);

	return 1;
}

