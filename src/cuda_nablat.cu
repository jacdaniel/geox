
#include <cuda.h>
#include <cuda_runtime.h>

#include <cuda_nablat.cuh>


__global__ void cuda_nablatx_kernel(double* in, int width, int depth, double* out)
{
	const long x = (long)blockDim.x * blockIdx.x + threadIdx.x;
	const long z = (long)blockDim.y * blockIdx.y + threadIdx.y;
	if (x < width && z < depth)
	{
		int add = width * z + x;
		if (x == 0)
		{
			out[add] = -in[add];
		}
		else if (x == width - 1)
		{
			out[add] = in[add - 1];
		}
		else
		{
			out[add] = in[add - 1] - in[add];
		}
	}
}

__global__ void cuda_nablatz_kernel(double* in, int width, int depth, double* out)
{
	const long x = (long)blockDim.x * blockIdx.x + threadIdx.x;
	const long z = (long)blockDim.y * blockIdx.y + threadIdx.y;
	if (x < width && z < depth)
	{
		int add = width * z + x;
		if (z == 0)
		{
			out[add] = -in[add];
		}
		else if (z == depth - 1)
		{
			out[add] = in[add - width];
		}
		else
		{
			out[add] = in[add - width] - in[add];
		}
	}
}

void cuda_nablatx_compute(double* in, int* size, double* out)
{
	dim3 block(32, 32);
	dim3 grid((size[0] - 1) / block.x + 1, (size[2] - 1) / block.y + 1);
	cuda_nablatx_kernel << <grid, block >> > (in, size[0], size[2], out);
	cudaDeviceSynchronize();
}

void cuda_nablatz_compute(double* in, int* size, double* out)
{
	dim3 block(32, 32);
	dim3 grid((size[0] - 1) / block.x + 1, (size[2] - 1) / block.y + 1);
	cuda_nablatz_kernel << <grid, block >> > (in, size[0], size[2], out);
	cudaDeviceSynchronize();
}


