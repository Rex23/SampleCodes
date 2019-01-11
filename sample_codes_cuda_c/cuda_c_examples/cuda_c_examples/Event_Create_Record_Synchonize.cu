#include "Event_Create_Record_Synchonize.cuh"

__global__ void saxpy_kernel(const int N, const float a, float* d_x, float* d_y)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i < N) {
		d_y[i] = a * d_x[i] + d_y[i];
	}
}

int Event_Create_Record_Synchonize(void)
{
	int N = 10000;

	float *x, *y, *d_x, *d_y;

	x = (float *)malloc(sizeof(float) * N);
	y = (float *)malloc(sizeof(float) * N);

	for (auto m = 0; m < N; m++) {
		x[m] = 1.0; y[m] = 1.0;
	}

	cudaMalloc(&d_x, sizeof(float)*N);
	cudaMalloc(&d_y, sizeof(float)*N);

	cudaMemcpy(d_x, x, sizeof(float) * N, cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, y, sizeof(float) * N, cudaMemcpyHostToDevice);

	const int ThreadsPerBlock = 256;
	dim3 grid((N + ThreadsPerBlock - 1) / ThreadsPerBlock);
	dim3 block(ThreadsPerBlock);

	cudaEvent_t start, stop;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start);

	saxpy_kernel <<< grid, block >>> (N, 2.0, d_x, d_y);

	cudaEventRecord(stop);

	cudaMemcpy(y, d_y, sizeof(float) * N, cudaMemcpyDeviceToHost);

	//for (auto m = 0; m < 10000; m++)
	//	cout << "Test y: " << y[m] << endl;

	cudaEventSynchronize(stop);

	float milliseconds = 0.0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	cout << "Time Elapsed: " << milliseconds << endl;

	free(x);
	free(y);
	cudaFree(d_x);
	cudaFree(d_y);

	return 1;
}