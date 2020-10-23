
#include <stdio.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include <cublas_v2.h>


typedef void (*CUDACG_CALLBACK)(double*, double*, void*);

typedef struct _CUDACG
{
	int no, nbiter;
	long size0;
	void* rhs, * x, * d_tmp, * d_z, * d_r, * d_d;
	void* callback, * preconditionner;
	void* callback_data, * preconditionner_data;
	cublasHandle_t handle;
}CUDACG;


void* cuda_cg_init()
{
	CUDACG* data = (CUDACG*)calloc(1, sizeof(CUDACG));
	return data;
}

void* cuda_cg_release(void* _p)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return NULL;
	free(p);
	return NULL;
}

void cuda_cg_set_callback(void* _p, void* callback, void* data)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return;
	p->callback = callback;
	p->callback_data = data;
}

void cuda_cg_set_preconditionner(void* _p, void* f, void* data)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return;
	p->preconditionner = f;
	p->preconditionner_data = data;
}

void cuda_cg_set_size(void* _p, int size)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return;
	p->size0 = size;
}

void cuda_cg_set_rhs(void* _p, void* rhs)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return;
	p->rhs = rhs;
}

void cuda_cg_set_x(void* _p, void* x)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return;
	p->x = x;
}

void cuda_cg_set_nbiter(void* _p, int nbiter)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return;
	p->nbiter = nbiter;
}


// =========================================================================

__global__ void cuda_cg_dot_kernel(double* v1, long size, double* v2, double* d_res)
{
	const long x = (long)blockDim.x * blockIdx.x + threadIdx.x;

	__shared__ double tmp[512], res0;
	double sum = 0.0;
	tmp[threadIdx.x] = 0.0;
	// if (threadIdx.x == 0) sum = *d_res;
	__syncthreads();


	if (x < size)
		tmp[threadIdx.x] = v1[x] * v2[x];
	__syncthreads();


	if (threadIdx.x == 0)
	{
		for (int i = 0; i < 512; i++)
			sum += tmp[i];
		*d_res += sum;
	}
	__syncthreads();
}

__global__ void cuda_cg_daxpy_kernel(double* x, long size, double alpha, double* y)
{
	const long add = (long)blockDim.x * blockIdx.x + threadIdx.x;
	if (add < size)
	{
		// y[add] += alpha * x[add];
		y[add] += alpha * x[add];
	}
}


void cuda_cg_daxpy_compute(double* x, long size, double alpha, double* y)
{
	dim3 block(32);
	dim3 grid((size - 1) / block.x + 1);
	// cuda_cg_daxpy_kernel << <grid, block >> > (x, size, alpha, y);
	// cudaDeviceSynchronize();
	double* h_y, * h_x;
	h_y = (double*)calloc(size, sizeof(double));
	h_x = (double*)calloc(size, sizeof(double));
	cudaMemcpy(h_y, y, size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_x, x, size * sizeof(double), cudaMemcpyDeviceToHost);
	for (long i = 0; i < size; i++)
		h_y[i] += alpha * h_x[i];
	cudaMemcpy(y, h_y, size * sizeof(double), cudaMemcpyHostToDevice);
	free(h_x);
	free(h_y);
}


#include "device_functions.h"
#include "device_launch_parameters.h"
#define BLOCK_SIZE 32
#define NO_SYNC
__global__ void dot(double* v1, int N, double* v2, double* out)
{
	__shared__ float cache[BLOCK_SIZE];
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	cache[threadIdx.x] = 0.f;
	while (i < N) {
		cache[threadIdx.x] += v1[i] * v2[i];
		i += gridDim.x * blockDim.x;
	}
	__syncthreads();  // required because later on the current thread is
					  // accessing data written by another thread
	i = BLOCK_SIZE / 2;
	while (i > 0) {
		if (threadIdx.x < i) cache[threadIdx.x] += cache[threadIdx.x + i];
		__syncthreads();
		i /= 2;  // not sure bitwise operations are actually faster
	}
#ifndef NO_SYNC  // serialized access to shared data;
	if (threadIdx.x == 0) atomicAdd(out, cache[0]);
#else  // no sync, what most likely happens is:
	// 1) all threads read 0
	// 2) all threads write concurrently 16 (local block dot product)
	if (threadIdx.x == 0) *out += cache[0];
#endif
}



void cuda_cg_dot(double* v1, long size, double* v2, double* res)
{
	dim3 block(32);
	dim3 grid((size - 1) / block.x + 1);

	double* d_res = NULL;
	cudaMalloc((void**)&d_res, 1 * sizeof(double));
	cudaMemset(d_res, 0, 1 * sizeof(double));
	// cuda_cg_dot_kernel << <grid, block >> > (v1, size, v2, d_res);
	dot << <grid, block >> > (v1, size, v2, d_res);

	cudaDeviceSynchronize();
	cudaMemcpy(res, d_res, 1 * sizeof(double), cudaMemcpyDeviceToHost);
}






void debug_print(double* x)
{
	double* xx = (double*)calloc(3, sizeof(double));
	cudaMemcpy(xx, x, 3 * sizeof(double), cudaMemcpyDeviceToHost);
	for (int i = 0; i < 3; i++)
		fprintf(stderr, "%f\n", xx[i]);
	free(xx);
}






int cuda_cg_run(void* _p)
{
	CUDACG* p = (CUDACG*)_p;
	if (p == NULL) return 1;
	double db = 0.0, rho0 = 0.0, rho1 = 0.0, beta = 0.0, denom = 1.0, alphax = 1.0;
	CUDACG_CALLBACK callback = (CUDACG_CALLBACK)p->callback;
	CUDACG_CALLBACK preconditionner = (CUDACG_CALLBACK)p->preconditionner;
	int cont = 1;
	double err = 0.0;

	const double minus_one = -1.0;
	const double plus_one = 1.0;

	if (p->d_tmp == NULL)
	{
		cudaMalloc(&p->d_tmp, p->size0 * sizeof(double));
		cudaMalloc(&p->d_z, p->size0 * sizeof(double));
		cudaMalloc(&p->d_r, p->size0 * sizeof(double));
		cudaMalloc(&p->d_d, p->size0 * sizeof(double));
	}
	if (p->handle == NULL)
	{
		cublasCreate_v2(&p->handle);
	}

	cublasDdot_v2(p->handle, p->size0, (const double*)p->rhs, 1, (const double*)p->rhs, 1, &db);
	cudaMemcpy(p->d_r, p->rhs, p->size0 * sizeof(double), cudaMemcpyDeviceToDevice);
	callback((double*)p->x, (double*)p->d_tmp, p->callback_data);
	cublasDaxpy(p->handle, p->size0, &minus_one, (double*)p->d_tmp, 1, (double*)p->d_r, 1);
	preconditionner((double*)p->d_r, (double*)p->d_z, p->preconditionner_data);
	cublasDdot(p->handle, p->size0, (const double*)p->d_r, 1, (const double*)p->d_z, 1, &rho0);
	int iter = 0;
	while (cont == 1)
	{
		if (iter % 50 == 0)
		{
			cublasDdot(p->handle, p->size0, (const double*)p->d_d, 1, (const double*)p->d_d, 1, &err);
			err = fabs(alphax) * sqrt(err / sqrt(p->size0));
			fprintf(stderr, "iter: %d - %d [ %f %f ]\n", iter, p->nbiter, alphax, err);
		}
		if (iter == 0)
		{
			cublasDcopy(p->handle, p->size0, (const double*)p->d_z, 1, (double*)p->d_d, 1);
		}
		else
		{
			beta = rho0 / rho1;
			cublasDcopy(p->handle, p->size0, (const double*)p->d_z, 1, (double*)p->d_tmp, 1);
			cublasDaxpy(p->handle, p->size0, &beta, (double*)p->d_d, 1, (double*)p->d_tmp, 1);
			cublasDcopy(p->handle, p->size0, (const double*)p->d_tmp, 1, (double*)p->d_d, 1);
		}

		callback((double*)p->d_d, (double*)p->d_tmp, p->callback_data);
		cublasDdot(p->handle, p->size0, (const double*)p->d_d, 1, (const double*)p->d_tmp, 1, &denom);
		alphax = rho0 / denom;
		cublasDaxpy_v2(p->handle, p->size0, &alphax, (double*)p->d_d, 1, (double*)p->x, 1);
		alphax = -alphax;
		cublasStatus_t  retx = cublasDaxpy_v2(p->handle, p->size0, &alphax, (double*)p->d_tmp, 1, (double*)p->d_r, 1);
		preconditionner((double*)p->d_r, (double*)p->d_z, p->preconditionner_data);
		rho1 = rho0;
		cublasDdot(p->handle, p->size0, (const double*)p->d_r, 1, (const double*)p->d_z, 1, &rho0);
		iter++;
		if (iter >= p->nbiter) cont = 0;
	}
	return 0;
}

