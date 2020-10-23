
#include <cuda.h>
#include <cuda_runtime.h>

#include <cuda_op.cuh>

__global__ void cuda_y_equals_a_plus_b_double_kernel(long size, double* a, double* b, double* y)
{
	const long add = (long)blockDim.x * blockIdx.x + threadIdx.x;

	if ( add < size )
	{
		y[add] = a[add] + b[add];
	}
}




void cuda_y_equals_a_plus_b_double(long size, double* a, double* b, double* y)
{
	dim3 block(1024);
	dim3 grid((size - 1) / block.x + 1);
	cuda_y_equals_a_plus_b_double_kernel << <grid, block >> > (size, a, b, y);
	cudaDeviceSynchronize();
}



__global__ void cuda_y_plus_equals_alpha_x_double_kernel(long size, double* y, double alpha, double* x)
{
	const long add = (long)blockDim.x * blockIdx.x + threadIdx.x;

	if (add < size)
	{
		y[add] += alpha* x[add];
	}
}




void cuda_y_plus_equals_alpha_x_double(long size, double* y, double alpha, double* x)
{
	dim3 block(1024);
	dim3 grid((size - 1) / block.x + 1);
	cuda_y_plus_equals_alpha_x_double_kernel << <grid, block >> > (size, y, alpha, x);
	cudaDeviceSynchronize();
}



__global__ void cuda_y_equal_a_times_b_kernel(long size, double* a, char* b, double* y)
{
	const long add = (long)blockDim.x * blockIdx.x + threadIdx.x;

	if (add < size)
	{
		y[add] = a[add] * (double)b[add];
	}
}


void cuda_y_equal_a_times_b(long size, double* a, char* b, double* y)
{
	dim3 block(1024);
	dim3 grid((size - 1) / block.x + 1);
	cuda_y_equal_a_times_b_kernel << <grid, block >> > (size, a, b, y);
	cudaDeviceSynchronize();
}


__global__ void cuda_y_minus_equals_x_kernel(long size, double* y, double x)
{
	const long add = (long)blockDim.x * blockIdx.x + threadIdx.x;

	if (add < size)
	{
		y[add] -= x;
	}
}




void cuda_y_minus_equals_x(long size, double* y, double x)
{
	dim3 block(1024);
	dim3 grid((size - 1) / block.x + 1);
	cuda_y_minus_equals_x_kernel << <grid, block >> > (size, y, x);
	cudaDeviceSynchronize();
}



__global__ 	void cuda_y_minus_equals_yref_kernel(long size, double* y, double yref)
{
	const long add = (long)blockDim.x * blockIdx.x + threadIdx.x;

	if (add < size)
	{
		y[add] -= yref;
	}
}

void cuda_y_minus_equals_yref(long size, double* y, long add)
{
	dim3 block(1024);
	dim3 grid((size - 1) / block.x + 1);
	double yref = 0.0;
	cudaMemcpy(&yref, y + add, 1 * sizeof(double), cudaMemcpyDeviceToHost);
	cuda_y_minus_equals_yref_kernel << <grid, block >> > (size, y, yref);
	cudaDeviceSynchronize();
}


double cuda_dot0(long size, double* x, double* y)
{
	double* x0 = (double*)calloc(size, sizeof(double));
	double* y0 = (double*)calloc(size, sizeof(double));

	cudaMemcpy(x0, x, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(y0, y, size * sizeof(double), cudaMemcpyDeviceToHost);
	double ret = 0.0;
	for (int i = 0; i < size; i++)
		ret += x0[i] * y0[i];
	free(x0);
	free(y0);
	return ret;
}