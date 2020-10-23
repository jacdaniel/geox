
#include <malloc.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include <cuda_nabla.cuh>

__global__ void cuda_nablax_kernel(double* in, int width, int depth, double* out)
{
	const long x = (long)blockDim.x * blockIdx.x + threadIdx.x;
	const long z = (long)blockDim.y * blockIdx.y + threadIdx.y;
	if (x < width && z < depth)
	{
		int add = width * z + x;
		if (x < width - 1)
			out[add] = in[add + 1] - in[add];
		else
			out[add] = 0.0;
	}
}

__global__ void cuda_nablaz_kernel(double* in, int width, int depth, double* out)
{
	const long x = (long)blockDim.x * blockIdx.x + threadIdx.x;
	const long z = (long)blockDim.y * blockIdx.y + threadIdx.y;
	if (x < width && z < depth)
	{
		int add = width * z + x;
		if (z < depth - 1)
			out[add] = in[add + width] - in[add];
		else
			out[add] = 0.0;
	}
}

void cuda_nablax_compute(double* in, int* size, double* out)
{
	dim3 block(32, 32);
	dim3 grid((size[0] - 1) / block.x + 1, (size[2] - 1) / block.y + 1);
	cuda_nablax_kernel << <grid, block >> > (in, size[0], size[2], out);
	cudaDeviceSynchronize();
}

void cuda_nablaz_compute(double* in, int* size, double* out)
{
	dim3 block(32, 32);
	dim3 grid((size[0] - 1) / block.x + 1, (size[2] - 1) / block.y + 1);
	cuda_nablaz_kernel << <grid, block >> > (in, size[0], size[2], out);
	cudaDeviceSynchronize();
}



void cuda_copy_double_float(long size, double* d_in, float* d_out)
{
	double* in = NULL;
	float* out = NULL;

	in = (double*)calloc(size, sizeof(double));
	out = (float*)calloc(size, sizeof(float));
	cudaMemcpy(in, d_in, size * sizeof(double), cudaMemcpyDeviceToHost);
	for (int i = 0; i < size; i++)
		out[i] = (float)in[i];
	cudaMemcpy(d_out, out, size * sizeof(float), cudaMemcpyHostToDevice);
	free(in);
	free(out);
}

void cuda_copy_float_double(long size, float* d_in, double* d_out)
{
	float* in = NULL;
	double* out = NULL;

	in = (float*)calloc(size, sizeof(float));
	out = (double*)calloc(size, sizeof(double));
	cudaMemcpy(in, d_in, size * sizeof(float), cudaMemcpyDeviceToHost);
	for (int i = 0; i < size; i++)
		out[i] = (double)in[i];
	cudaMemcpy(d_out, out, size * sizeof(double), cudaMemcpyHostToDevice);
	free(in);
	free(out);
}