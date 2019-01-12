/* Copyright (c) 1993-2015, NVIDIA CORPORATION. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*  * Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*  * Neither the name of NVIDIA CORPORATION nor the names of its
*    contributors may be used to endorse or promote products derived
*    from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
* PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
* OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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