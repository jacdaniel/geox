
#include <stdio.h>

#include <malloc.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include <cuda_nabla.cuh>
#include <cuda_nablat.cuh>
#include <cuda_op.cuh>
#include <cuda_cg.cuh>
#include <cuda_poisson_fourier.cuh>


void debug_ddata_save(double* d_data, int* size, char* filename)
{
	long size0 = size[0] * size[2];
	double* data = (double*)calloc(size0, sizeof(double));
	cudaMemcpy(data, d_data, size0 * sizeof(double), cudaMemcpyDeviceToHost);
	FILE* pFile = fopen(filename, "wb");
	fwrite(data, sizeof(double), size0, pFile);
	fclose(pFile);
	free(data);
}

void debug_ddata_save(char* d_data, int* size, char* filename)
{
	long size0 = size[0] * size[2];
	char* tmp = (char*)calloc(size0, sizeof(char));
	double* data = (double*)calloc(size0, sizeof(double));
	cudaMemcpy(tmp, d_data, size0 * sizeof(char), cudaMemcpyDeviceToHost);
	for (int i = 0; i < size0; i++) data[i] = tmp[i];
	FILE* pFile = fopen(filename, "wb");
	fwrite(data, sizeof(double), size0, pFile);
	fclose(pFile);
	free(data);
	free(tmp);
}

void debug_ddata_save(float* d_data, int* size, char* filename)
{
	long size0 = size[0] * size[2];
	float* tmp = (float*)calloc(size0, sizeof(float));
	double* data = (double*)calloc(size0, sizeof(double));
	cudaMemcpy(tmp, d_data, size0 * sizeof(float), cudaMemcpyDeviceToHost);
	for (int i = 0; i < size0; i++) data[i] = tmp[i];
	FILE* pFile = fopen(filename, "wb");
	fwrite(data, sizeof(double), size0, pFile);
	fclose(pFile);
	free(data);
	free(tmp);
}



typedef struct _PARAM
{
	int size[3];
	char* d_k;
	double* d_deltatau, * d_tau, *d_px, *d_pz, *d_r0x, *d_r0z, *d_rx, *d_rz, *d_div, *d_ktau, *d_tau0, *d_divx, *d_divz, *d_zk, *d_pk;
	float* d_in_float, * d_out_float;
	void* poisson;
}PARAM;


static PARAM* param_init(int* size)
{
	PARAM* param = (PARAM*)calloc(1, sizeof(PARAM));
	param->size[0] = size[0];
	param->size[1] = size[1];
	param->size[2] = size[2];

	int xc[] = { size[0] / 3 , 2 * size[0] / 3 , 5};
	int zc[] = { size[2] / 3 , 2 * size[2] / 3, 5 };
	float vc[] = { 10.0f, 10.0f, 50.0f };
	int Nc = 3;
	
	int size0 = param->size[0] * param->size[2];
	cudaMalloc((void**)&param->d_k, size0 * sizeof(char));
	cudaMalloc((void**)&param->d_tau, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_deltatau, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_px, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_pz, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_rx, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_rz, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_div, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_r0x, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_r0z, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_divx, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_divz, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_zk, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_ktau, size0 * sizeof(double));
	cudaMalloc((void**)&param->d_pk, size0 * sizeof(double));


	cudaMalloc((void**)&param->d_in_float, size0 * sizeof(float));
	cudaMalloc((void**)&param->d_out_float, size0 * sizeof(float));


	char* k = (char*)calloc(size0, sizeof(char));
	for (int i = 0; i < size0; i++) k[i] = 1;
	for (int n = 0; n < Nc; n++)
		k[size[0] * zc[n] + xc[n]] = 0;
	cudaMemcpy(param->d_k, k, size0 * sizeof(char), cudaMemcpyHostToDevice);
	free(k);
	debug_ddata_save(param->d_k, size, "d:\\k0.raw");


	double* tau0 = (double*)calloc(size0, sizeof(double));
	double* tau = (double*)calloc(size0, sizeof(double));

	for (int n = 0; n < size0; n++)
		tau[n] = vc[0];
	for (int n = 0; n < Nc; n++)
	{
		tau[size[0] * zc[n] + xc[n]] = vc[n];
	}

	for (int n = 0; n < size0; n++)
	{
		tau[n] -= vc[0];
	}


	for (int n = 0; n < size0; n++)
	{
		tau0[n] = tau[n]*0;
	}

	for (int n = 0; n < Nc; n++)
		tau0[size[0] * zc[n] + xc[n]] = vc[n]-vc[0];



	cudaMalloc((void**)&param->d_tau0, size0 * sizeof(double));

	cudaMemcpy(param->d_tau0, tau0, size0 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(param->d_tau, tau, size0 * sizeof(double), cudaMemcpyHostToDevice);

	debug_ddata_save(param->d_tau0, size, "d:\\tau0.raw");
	debug_ddata_save(param->d_tau, size, "d:\\tau.raw");

	cuda_nablax_compute(param->d_tau0, size, param->d_r0x);
	cuda_nablaz_compute(param->d_tau0, size, param->d_r0z);

	debug_ddata_save(param->d_r0x, size, "d:\\r0x.raw");
	debug_ddata_save(param->d_r0z, size, "d:\\r0z.raw");



	int ret = cuda_poisson_fourier_init(&param->poisson, 0, size[0], size[2]);
	return param;
}











typedef struct _CALLBACK_DATA3
{
	int* size;
	void* poisson;
	char* d_k;
	float* d_in_float, * d_out_float;

}CALLBACK_DATA3;


void callback3(double* x, double* y, void* data)
{
	CALLBACK_DATA3* param = (CALLBACK_DATA3*)data;
	int size0 = param->size[0] * param->size[2];

	double* d_tauk = NULL, *d_gx = NULL, *d_gz = NULL, *d_gxt = NULL, *d_gzt = NULL;
	cudaMalloc((void**)&d_tauk, size0 * sizeof(double));
	cudaMalloc((void**)&d_gx, size0 * sizeof(double));
	cudaMalloc((void**)&d_gz, size0 * sizeof(double));
	cudaMalloc((void**)&d_gxt, size0 * sizeof(double));
	cudaMalloc((void**)&d_gzt, size0 * sizeof(double));

	// debug_ddata_save(x, param->size, "d:\\res.raw");

	cuda_y_equal_a_times_b(size0, x, param->d_k, d_tauk);

	// debug_ddata_save(d_tauk, param->size, "d:\\res.raw");

	cuda_nablax_compute(d_tauk, param->size, d_gx);
	cuda_nablaz_compute(d_tauk, param->size, d_gz);

	cuda_nablatx_compute(d_gx, param->size, d_gxt);
	cuda_nablatz_compute(d_gz, param->size, d_gzt);

	cuda_y_plus_equals_alpha_x_double(size0, d_gxt, 1.0, d_gzt);
	cuda_y_equal_a_times_b(size0, d_gxt, param->d_k, y);

	// debug_ddata_save(y, param->size, "d:\\res.raw");



	cudaFree(d_tauk);
	cudaFree(d_gx);
	cudaFree(d_gz);
	cudaFree(d_gxt);
	cudaFree(d_gzt);
}

void preconditionner3(double* x, double* y, void* data)
{
	CALLBACK_DATA3* param = (CALLBACK_DATA3*)data;
	int size0 = param->size[0] * param->size[2];

	// cudaMemcpy(y, x, size0 * sizeof(double), cudaMemcpyDeviceToDevice);	
	
	
	cuda_copy_double_float(size0, x, param->d_in_float);
	debug_ddata_save(param->d_in_float, param->size, "d:\\res.raw");
	cuda_poisson_fourier_do(param->poisson, param->d_in_float, param->size[0], param->size[2], param->d_out_float);
	cuda_copy_float_double(size0, param->d_out_float, y);
	
	
	
	
	debug_ddata_save(param->d_out_float, param->size, "d:\\res.raw");
	

	// cuda_y_minus_equals_yref(size0, y, param->size[0] * (param->size[2] / 3) + (param->size[0] / 3));
	// return;
	double* tmp = NULL;
	tmp = (double*)calloc(size0, sizeof(double));
	cudaMemcpy(tmp, y, size0 * sizeof(double), cudaMemcpyDeviceToHost);
	tmp[param->size[0] * (param->size[2] / 3) + param->size[0] / 3] = 0.0f;
	tmp[param->size[0] * (2 * param->size[2] / 3) + 2 * param->size[0] / 3] = 0.0f;
	tmp[param->size[0] * 5 + 5] = 0.0f;
	cudaMemcpy(y, tmp, size0 * sizeof(double), cudaMemcpyHostToDevice);
	// debug_ddata_save(y, param->size, "d:\\res.raw");

}






void main_constraint_geotime_run(int argc, char** argv)
{
	int size[] = { 100, 101, 102 };
	int size0 = size[0] * size[2];

	PARAM* param = param_init(size);

	// cudaMemcpy(param->d_tau, param->d_tau0, size0 * sizeof(double), cudaMemcpyDeviceToDevice);
	debug_ddata_save(param->d_tau, param->size, "d:\\res.raw");
	cuda_y_equal_a_times_b(size0, param->d_tau, param->d_k, param->d_ktau);
	debug_ddata_save(param->d_ktau, param->size, "d:\\res.raw");

	cuda_nablax_compute(param->d_ktau, size, param->d_rx);
	cuda_y_minus_equals_x(size0, param->d_rx, -1.0);
	debug_ddata_save(param->d_rx, param->size, "d:\\res.raw");

	debug_ddata_save(param->d_r0x, param->size, "d:\\res.raw");
	cuda_y_plus_equals_alpha_x_double(size0, param->d_rx, 1.0, param->d_r0x);
	debug_ddata_save(param->d_rx, param->size, "d:\\res.raw");


	cuda_nablaz_compute(param->d_ktau, size, param->d_rz);
	cuda_y_minus_equals_x(size0, param->d_rz, 2.0);
	cuda_y_plus_equals_alpha_x_double(size0, param->d_rz, 1.0, param->d_r0z);
	debug_ddata_save(param->d_rz, param->size, "d:\\res.raw");

	cuda_nablatx_compute(param->d_rx, size, param->d_divx);
	cuda_nablatz_compute(param->d_rz, size, param->d_divz);

	cuda_y_plus_equals_alpha_x_double(size0, param->d_divx, 1.0, param->d_divz);
	cuda_y_equal_a_times_b(size0, param->d_divx, param->d_k, param->d_divx);
	
	debug_ddata_save(param->d_divx, param->size, "d:\\res.raw");


	CALLBACK_DATA3* param0 = (CALLBACK_DATA3*)calloc(1, sizeof(CALLBACK_DATA3));
	param0->size = size;
	param0->poisson = param->poisson;
	param0->d_k = param->d_k;
	param0->d_in_float = param->d_in_float;
	param0->d_out_float = param->d_out_float;

	double* d_xx = NULL;
	cudaMalloc((void**)&d_xx, size0 * sizeof(double));
	cudaMemset(d_xx, 0, size0 * sizeof(double));

	void* p = cuda_cg_init();
	cuda_cg_set_callback(p, callback3, param0);
	cuda_cg_set_preconditionner(p, preconditionner3, param0);
	cuda_cg_set_size(p, size0);
	cuda_cg_set_rhs(p, param->d_divx);
	cuda_cg_set_x(p, d_xx);
	cuda_cg_set_nbiter(p, 100);
	cuda_cg_run(p);
	p = cuda_cg_release(p);


	cuda_y_plus_equals_alpha_x_double(size0, param->d_tau, -1.0, d_xx);


	debug_ddata_save(param->d_tau, param->size, "d:\\res.raw");






	



}
