
#ifndef __CUDA_POISSON_FOURIER__
#define __CUDA_POISSON_FOURIER__


#ifdef __cplusplus
extern "C" {
#endif

	int cuda_poisson_fourier_init(void** ptr, int type, int width, int depth);

	void* cuda_poisson_fourier_release(void* _data);

	void cuda_poisson_fourier_do(void* _data, float* d_in, int width, int depth, float* d_out);


#ifdef __cplusplus
}
#endif



#endif