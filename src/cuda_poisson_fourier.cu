
#include <stdio.h>
#include <malloc.h>

#include <cufft.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cudaProfiler.h>
#include <device_functions.h>

// #include <cuda_utils.h>
// #include <util.h>

#define SUCCESS 0
#define FAIL    1

#include <cuda_poisson_fourier.cuh>

#if !defined (PI)
#define PI				 3.14159265358979
#endif

#if !defined (PI_2)
#define PI_2  1.570796326794897
#endif

#ifndef CUDAFREE
#define CUDAFREE(ptr) { if ( (ptr) != NULL ) { cudaFree(ptr); (ptr) = NULL; } }
#endif

#define CUFFTDESTROY(plan) { if ( plan ) { cufftDestroy(plan); plan = 0; } }

#ifndef FREE
#define FREE(x){ \
if ( x != NULL ) { \
        free(x); \
        x = NULL;} }
#endif


#define CUDA_ALLOC_SUCCESS 0
#define CUDA_ALLOC_FAIL 1

int cufftplan_chk(cufftResult_t  code, const char* file, int line, int abort)
{
    int ret = CUDA_ALLOC_SUCCESS;
    if (code != cudaSuccess)
    {
        ret = CUDA_ALLOC_FAIL;
        fprintf(stderr, "cufftplan error: %d %s %d\n", code, file, line);
        if (abort) exit(code);
        exit(0);
    }
    return ret;
    // fprintf(stderr,"OK fftplan GPUassert: %d %s %d\n", code, file, line);
}
#define CUFFTPLAN_CHK(ans) cufftplan_chk((ans), __FILE__, __LINE__, 0);


int gpuAssert(cudaError_t code, const char* file, int line, int abort)
{
    int ret = CUDA_ALLOC_SUCCESS;
    if (code != cudaSuccess)
    {
        ret = CUDA_ALLOC_FAIL;
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
        exit(0);
    }
    return ret;
    // fprintf(stderr,"OK GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
}
#define gpuErrchk(ans) gpuAssert((ans), __FILE__, __LINE__, 0);


typedef struct _CUDA_POISSON_FOURIER
{
    cufftHandle cufft_p1, cufft_p2;
    int type, width, depth;
    float* d_dataPre, * d_cos_a, * d_sin_a, * d_cos_b, * d_sin_b, * d_dct, * d_inv_lap, size_poisson, * d_out2;
    cufftComplex* d_freq, * W, * d_fftData;
    cufftReal* d_ifftData, * d_fftData2;
} CUDA_POISSON_FOURIER;


static __global__ void cuda_inverse_laplacian(long x1, long x2, long z1, long z2, long dxt, long dzt, float* filter)
{
    long x, z, dx;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    z = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    dx = x2 + 1 - x1;
    // dz = z2+1-z1;

    if (x >= x1 && x <= x2 && z >= z1 && z <= z2)
        if (x == 0 && z == 0)
            filter[0] = 0.0;
        else
            filter[(z - z1) * dx + x - x1] = (float)(1.0f / (-4.0f + 2.0f * (cosf(2.0f * PI * (float)x / (float)dxt) + cosf(2.0f * PI * (float)z / (float)dzt))));
}


// DCT
static __global__ void cuda_DCTPreprocessing(float* in, long width, long depth, float* out)
{
    long x, z, add;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    //z = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    z = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width && z < depth)
    {
        add = width * z + x;
        if (x % 2 == 0) x = x >> 1; else x = width - ((x + 1) >> 1);
        if (z % 2 == 0) z = z >> 1; else z = depth - ((z + 1) >> 1);
        out[z * width + x] = in[add];
    }
}

static __global__ void cuda_DCTPostprocessingWithNorm(float* in, long width, long depth, float norm, float* out)
{
    long x, z, add;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    // z = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    z = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width && z < depth)
    {
        add = width * z + x;
        if (x % 2 == 0) x = x >> 1; else x = width - ((x + 1) >> 1);
        if (z % 2 == 0) z = z >> 1; else z = depth - ((z + 1) >> 1);
        out[add] = in[z * width + x] / norm;
    }
}

static __global__ void cuda_sincos_generate(long width, long height, float* d_cos_a, float* d_sin_a, float* d_cos_b, float* d_sin_b)
{
    long x, y, add;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    // y = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    y = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width && y < height)
    {
        add = width * y + x;
        d_cos_a[add] = (float)cosf(-(float)PI_2 * x / width - (float)PI_2 * y / height);
        d_sin_a[add] = (float)sinf(-(float)PI_2 * x / width - (float)PI_2 * y / height);
        d_cos_b[add] = (float)cosf(-(float)PI_2 * x / width + (float)PI_2 * y / height);
        d_sin_b[add] = (float)sinf(-(float)PI_2 * x / width + (float)PI_2 * y / height);
    }
}

static __global__ void cuda_fftRecombinaison(float* dct, long width, long height, cufftComplex* W, cufftComplex* fftData)
{
    long x, y, add, width2 = width >> 1;
    float re, im;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    // y = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    y = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width2 + 1 && y < height)
    {
        add = width * y + x;
        re = dct[(width + 1) * y + x] - dct[(width + 1) * (height - y) + width - x];
        im = dct[(width + 1) * y + width - x] + dct[(width + 1) * (height - y) + x];

        fftData[(width2 + 1) * y + x].x = W[add].x * re + W[add].y * im;
        fftData[(width2 + 1) * y + x].y = W[add].y * re - W[add].x * im;
    }
}

static __global__ void cuda_w_create(long width, long height, cufftComplex* data)
{
    long x, y;
    float arg;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    // y = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    y = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width && y < height)
    {
        arg = (float)PI_2 * ((float)x / width + (float)y / height);
        data[width * y + x].x = cos(arg);
        data[width * y + x].y = sin(arg);
    }
}

static __global__ void cuda_DFT2DCT(cufftComplex* in, long width, long height, long width_dft, float* cos_a, float* sin_a, float* cos_b, float* sin_b, float* out)
{
    long x, y, addc, add1, add2, width2 = width >> 1, add_out;
    float re1, im1, re2, im2;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    // y = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    y = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width && y < height)
    {
        addc = width * y + x;
        add_out = (width + 1) * y + x;
        if (y == 0)
        {
            if (x <= width2)
            {
                add1 = x;
                re1 = in[add1].x;
                im1 = in[add1].y;
            }
            else
            {
                add1 = width - x;
                re1 = in[add1].x;
                im1 = -in[add1].y;
            }
            out[add_out] = (float)4.0f * (re1 * cos_a[addc] - im1 * sin_a[addc]);
        }
        else
        {
            if (x <= width2)
            {
                add1 = width_dft * y + x;
                add2 = width_dft * (height - y) + x;
                re1 = in[add1].x; im1 = in[add1].y;
                re2 = in[add2].x; im2 = in[add2].y;
            }
            else
            {
                add1 = width_dft * (height - y) + width - x;
                add2 = width_dft * y + width - x;
                re1 = in[add1].x; im1 = -in[add1].y;
                re2 = in[add2].x; im2 = -in[add2].y;
            }
            out[add_out] = (float)2.0f * (re1 * cos_a[addc] - im1 * sin_a[addc] + re2 * cos_b[addc] - im2 * sin_b[addc]);
        }
    }
}


static __global__ void cuda_LapMult(float* r, float* lap, long size)
{
    long add;

    // add = blockDim.x * blockIdx.x + threadIdx.x;
    add = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;

    if (add < size)
        r[add] *= lap[add];
}

static __global__ void cuda_LapMult2(cufftComplex* r, float* lap, long size)
{
    long add;
    // add = blockDim.x * blockIdx.x + threadIdx.x;

    add = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;

    if (add < size)
    {
        r[add].x *= lap[add];
        r[add].y *= lap[add];
    }
}


static __global__ void cuda_fft_fill(float* in, long width, long depth, float* fftData)
{
    long x, z;
    float val;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    // y = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    z = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width && z < depth)
    {
        val = in[width * z + x];
        fftData[2 * width * z + x] = val;
        fftData[2 * width * (2 * depth - 1 - z) + x] = val;
        fftData[2 * width * z + 2 * width - 1 - x] = val;
        fftData[2 * width * (2 * depth - 1 - z) + 2 * width - 1 - x] = val;
    }
}

static __global__ void cuda_FirstQuarterWithNorm(float* out, int width, int depth, float norm, float* in)
{
    long x, z;

    // x = blockDim.x * blockIdx.x + threadIdx.x;
    // y = blockDim.y * blockIdx.y + threadIdx.y;
    x = __umul24(blockDim.x, blockIdx.x) + threadIdx.x;
    z = __umul24(blockDim.y, blockIdx.y) + threadIdx.y;

    if (x < width && z < depth)
    {
        out[width * z + x] = in[2 * width * z + x] / norm;
    }
}

int cuda_poisson_fourier_init(void** ptr, int type, int width, int depth)
{
    CUDA_POISSON_FOURIER* data = NULL;
    int size, threads_per_block, block_countX, block_countY;
    int ret = SUCCESS;

    if (ptr == NULL) return FAIL;
    *ptr = (CUDA_POISSON_FOURIER*)calloc(1, sizeof(CUDA_POISSON_FOURIER));
    data = (CUDA_POISSON_FOURIER*)(*ptr);

    data->type = type;
    data->width = width;
    data->depth = depth;
    size = width * depth;
    if (type == 0)
    {
        ret = CUFFTPLAN_CHK(cufftPlan2d(&data->cufft_p1, depth, width, CUFFT_R2C)); if (ret != SUCCESS) return FAIL;
        ret = CUFFTPLAN_CHK(cufftPlan2d(&data->cufft_p2, depth, width, CUFFT_C2R)); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_freq, (width / 2 + 1) * depth * sizeof(cufftComplex))); if (ret != SUCCESS) return FAIL;
        // cudaMalloc(&data->d_fftData, (width/2+1)*depth*sizeof(cufftComplex));
        ret = gpuErrchk(cudaMalloc(&data->d_dct, (width + 1) * (depth + 1) * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_dataPre, size * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_ifftData, size * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_inv_lap, (width + 1) * (depth + 1) * sizeof(float))); if (ret != SUCCESS) return FAIL;
        threads_per_block = 32;
        block_countX = (width + threads_per_block - 1) / threads_per_block;
        block_countY = (depth + threads_per_block - 1) / threads_per_block;
        dim3 blockXY(threads_per_block, threads_per_block);
        dim3 gridXY(block_countX, block_countY);
        cuda_inverse_laplacian << <gridXY, blockXY >> > (0, width, 0, depth, 2 * width, 2 * depth, data->d_inv_lap); cudaDeviceSynchronize();

        ret = gpuErrchk(cudaMalloc(&data->d_cos_a, size * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_sin_a, size * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_cos_b, size * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_sin_b, size * sizeof(float))); if (ret != SUCCESS) return FAIL;

        cuda_sincos_generate << <gridXY, blockXY >> > (width, depth, data->d_cos_a, data->d_sin_a, data->d_cos_b, data->d_sin_b); cudaDeviceSynchronize();

        ret = gpuErrchk(cudaMalloc(&data->W, size * sizeof(cufftComplex))); if (ret != SUCCESS) return FAIL;
        cuda_w_create << <gridXY, blockXY >> > (width, depth, data->W); cudaDeviceSynchronize();
        data->size_poisson = (float)(4 * width * depth);
    }
    else
    {
        ret = gpuErrchk(cudaMalloc(&data->d_fftData2, 4 * width * depth * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_freq, (width + 1) * 2 * depth * sizeof(cufftComplex))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_ifftData, 4 * width * depth * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_out2, 4 * width * depth * sizeof(float))); if (ret != SUCCESS) return FAIL;
        ret = gpuErrchk(cudaMalloc(&data->d_inv_lap, (width + 1) * (2 * depth) * sizeof(float))); if (ret != SUCCESS) return FAIL;
        threads_per_block = 32;
        block_countX = (width + 1 + threads_per_block - 1) / threads_per_block;
        block_countY = (2 * depth + 1 + threads_per_block - 1) / threads_per_block;
        dim3 blockXY3(threads_per_block, threads_per_block);
        dim3 gridXY3(block_countX, block_countY);
        cuda_inverse_laplacian << <gridXY3, blockXY3 >> > (0, width, 0, 2 * depth - 1, 2 * width, 2 * depth, data->d_inv_lap); cudaDeviceSynchronize();

        ret = CUFFTPLAN_CHK(cufftPlan2d(&data->cufft_p1, 2 * depth, 2 * width, CUFFT_R2C)); if (ret != SUCCESS) return FAIL;
        ret = CUFFTPLAN_CHK(cufftPlan2d(&data->cufft_p2, 2 * depth, 2 * width, CUFFT_C2R)); if (ret != SUCCESS) return FAIL;
        data->size_poisson = (float)(4 * width * depth);
    }
    return SUCCESS;
}

void* cuda_poisson_fourier_release(void* _data)
{
    CUDA_POISSON_FOURIER* data = (CUDA_POISSON_FOURIER*)_data;

    if (data == NULL) return NULL;
    CUDAFREE(data->d_inv_lap)
        CUDAFREE(data->d_ifftData)
        CUDAFREE(data->d_dct)
        CUDAFREE(data->d_cos_a)
        CUDAFREE(data->d_sin_a)
        CUDAFREE(data->d_cos_b)
        CUDAFREE(data->d_sin_b)
        CUDAFREE(data->W)
        CUDAFREE(data->d_dataPre)
        CUDAFREE(data->d_fftData)
        CUDAFREE(data->d_out2)
        CUDAFREE(data->d_fftData2)
        CUFFTDESTROY(data->cufft_p1)
        CUFFTDESTROY(data->cufft_p2)
        FREE(data)
        return NULL;
}



void cuda_poisson_fourier_do(void* _data, float* d_in, int width, int depth, float* d_out)
{
    CUDA_POISSON_FOURIER* data = (CUDA_POISSON_FOURIER*)_data;
    int threads_per_block, block_count, block_countX, block_countY;

    if (data == NULL) return;

    threads_per_block = 1024;
    block_count = (width * depth + threads_per_block - 1) / threads_per_block;
    dim3 block(threads_per_block);
    dim3 grid(block_count);

    threads_per_block = 32;
    block_countX = (width + threads_per_block - 1) / threads_per_block;
    block_countY = (depth + threads_per_block - 1) / threads_per_block;
    dim3 blockXY(threads_per_block, threads_per_block);
    dim3 gridXY(block_countX, block_countY);

    threads_per_block = 1024;
    block_count = ((width + 1) * (depth + 1) + threads_per_block - 1) / threads_per_block;
    dim3 blockLap(threads_per_block);
    dim3 gridLap(block_count);

    threads_per_block = 1024;
    block_count = ((width + 1) * (2 * depth) + threads_per_block - 1) / threads_per_block;
    dim3 blockLap2(threads_per_block);
    dim3 gridLap2(block_count);

    // data->type = 1;
    if (data->type == 0)
    {
        cudaMemset(data->d_dct, 0, (width + 1) * (depth + 1) * sizeof(float));
        cuda_DCTPreprocessing << <gridXY, blockXY >> > ((float*)d_in, width, depth, data->d_dataPre); cudaDeviceSynchronize(); //cudaThreadSynchronize();
        cufftExecR2C(data->cufft_p1, (cufftReal*)data->d_dataPre, (cufftComplex*)data->d_freq); // cudaThreadSynchronize();
        cuda_DFT2DCT << <gridXY, blockXY >> > ((cufftComplex*)data->d_freq, width, depth, width / 2 + 1, data->d_cos_a, data->d_sin_a, data->d_cos_b, data->d_sin_b, data->d_dct); cudaDeviceSynchronize(); //cudaThreadSynchronize();
        cuda_LapMult << <gridLap, blockLap >> > (data->d_dct, data->d_inv_lap, (width + 1) * (depth + 1)); cudaDeviceSynchronize(); //cudaThreadSynchronize();
        cuda_fftRecombinaison << <gridXY, blockXY >> > (data->d_dct, width, depth, data->W, data->d_freq); cudaDeviceSynchronize(); //cudaThreadSynchronize();
        cufftExecC2R(data->cufft_p2, (cufftComplex*)data->d_freq, (cufftReal*)data->d_ifftData); // cudaThreadSynchronize();
        cuda_DCTPostprocessingWithNorm << <gridXY, blockXY >> > (data->d_ifftData, width, depth, data->size_poisson, d_out); cudaDeviceSynchronize(); //cudaThreadSynchronize();
    }
    else
    {
        cuda_fft_fill << <gridXY, blockXY >> > (d_in, width, depth, data->d_fftData2); cudaDeviceSynchronize(); //cudaThreadSynchronize();
        cufftExecR2C(data->cufft_p1, (cufftReal*)data->d_fftData2, (cufftComplex*)data->d_freq); // cudaThreadSynchronize();
        cuda_LapMult2 << <gridLap2, blockLap2 >> > ((cufftComplex*)data->d_freq, data->d_inv_lap, (width + 1) * (2 * depth)); cudaDeviceSynchronize(); //cudaThreadSynchronize();
        cufftExecC2R(data->cufft_p2, (cufftComplex*)data->d_freq, (cufftReal*)data->d_out2); // cudaThreadSynchronize();
        cuda_FirstQuarterWithNorm << <gridXY, blockXY >> > (d_out, width, depth, data->size_poisson, data->d_out2); cudaDeviceSynchronize(); //cudaThreadSynchronize();
    }

}
