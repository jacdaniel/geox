
#ifdef __LINUX__
// #include <chrono>
// #include <ctime>
#include <sys/time.h>
#include <cstdlib>
#endif
#ifdef __WINDOWS__
#include <windows.h>
#endif

#include <string.h>
#include <malloc.h>
#include <chronos.h>

#ifdef _CUDA_
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#ifdef __ENABLE_CHRONOS__

#ifdef __LINUX__

typedef struct _CHRONOS
{
	timeval start, end;
	double time;
} CHRONOS;


void *chronos_init(int nbre, const char *base)
{
	CHRONOS *chronos = NULL;

	chronos = (CHRONOS*)calloc(nbre, sizeof(CHRONOS));
	return chronos;
}

void *chronos_release(void *_chronos)
{
	// CHRONOS *chronos = (CHRONOS*)_chronos;
	FREE(_chronos)
	return NULL;
}

void chronos_tic(void *_chronos, int no)
{
	CHRONOS *chronos = (CHRONOS*)_chronos;

	if (chronos == NULL) return;
	gettimeofday(&chronos[no].start, NULL);
}

void chronos_toc(void *_chronos, int no)
{
	CHRONOS *chronos = (CHRONOS*)_chronos;
	long seconds = 0, useconds = 0;

	if (chronos == NULL) return;

	gettimeofday(&chronos[no].end, NULL);

	seconds = chronos[no].end.tv_sec - chronos[no].start.tv_sec;
	useconds = chronos[no].end.tv_usec - chronos[no].start.tv_usec;
	chronos[no].time += seconds*1000.0 + useconds / 1000.0;
}

void chronos_reset(void *_chronos, int no)
{
	CHRONOS *chronos = (CHRONOS*)_chronos;

	if (chronos == NULL) return;
	chronos[no].time = 0.0;
}

double chronos_get(void *_chronos, int no)
{
	CHRONOS *chronos = (CHRONOS*)_chronos;

	if (chronos == NULL) return 0.0;
	return chronos[no].time;
}
#endif

#ifdef __WINDOWS__
typedef struct _HR_CLOCK
{
	LARGE_INTEGER counter;
	LARGE_INTEGER frequency;
} HR_CLOCK;


#define DWORD_MAX 4294967296.0


void * attHRClockNew()
{
	HR_CLOCK * hr_clock = NULL;
	LARGE_INTEGER frequency;

	if (!QueryPerformanceFrequency(&frequency)) return NULL;
	hr_clock = (HR_CLOCK *)calloc(1, sizeof(HR_CLOCK));
	hr_clock->frequency = frequency;
	QueryPerformanceCounter(&(hr_clock->counter));
	return hr_clock;
}


void * attHRClockRelease(void * hr_clock)
{
	if (hr_clock == NULL) return NULL;
	free(hr_clock);
	return NULL;
}


static double attHRClockCompute(LARGE_INTEGER * counter1, LARGE_INTEGER * counter2, LARGE_INTEGER * frequency)
{
	double dt, f;

	dt = (double)counter2->HighPart - (double)counter1->HighPart;
	dt *= DWORD_MAX;
	dt += (double)counter2->LowPart - (double)counter1->LowPart;
	f = (double)frequency->HighPart*DWORD_MAX + (double)frequency->LowPart;
	return dt / f*1000.0;
}


double attHRClockDifferenceTime(void * hr_clock_begin, void * hr_clock_end)
{
	HR_CLOCK * hr_clock1, *hr_clock2;

	if (hr_clock_begin == NULL) return 0.0;
	if (hr_clock_end == NULL) return 0.0;
	hr_clock1 = (HR_CLOCK *)hr_clock_begin;
	hr_clock2 = (HR_CLOCK *)hr_clock_end;
	return attHRClockCompute(&(hr_clock1->counter), &(hr_clock2->counter), &(hr_clock1->frequency));
}


double attHRClockElapsedTime(void * hr_clock_begin)
{
	LARGE_INTEGER counter;
	HR_CLOCK * hr_clock1;

	if (hr_clock_begin == NULL) return 0.0;
	QueryPerformanceCounter(&counter);
	hr_clock1 = (HR_CLOCK *)hr_clock_begin;
	return attHRClockCompute(&(hr_clock1->counter), &counter, &(hr_clock1->frequency));
}


typedef struct _CHRONOS
{
	HR_CLOCK ** hr_clock;
	double * tic;
	double * toc;
	double * dt;
} CHRONOS;


void * chronos_init(int nbre, char const *base)
{
	CHRONOS * chr = NULL;
	long n;

	chr = (CHRONOS *)calloc(1, sizeof(CHRONOS));
	chr->hr_clock = (HR_CLOCK **)calloc(100, sizeof(HR_CLOCK *));
	for (n = 0; n<100; n++) chr->hr_clock[n] = (HR_CLOCK*)attHRClockNew();
	chr->tic = (double *)calloc(100, sizeof(double));
	chr->toc = (double *)calloc(100, sizeof(double));
	chr->dt = (double *)calloc(100, sizeof(double));
	return (void *)chr;
}


void * chronos_release(void * chronos)
{
	CHRONOS * chr = (CHRONOS *)chronos;
	long n;

	if (chr == NULL) return NULL;
	for (n = 0; n<100; n++) attHRClockRelease(chr->hr_clock[n]);
	free(chr->hr_clock);
	free(chr->tic);
	free(chr->toc);
	free(chr->dt);
	free(chr);
	return NULL;
}


void chronos_reset(void * chronos, int no)
{
	CHRONOS * chr = (CHRONOS *)chronos;

	if (chr == NULL || no >= 100) return;
	chr->dt[no] = 0.0;
}


void chronos_tic(void * chronos, int no)
{
	CHRONOS * chr = (CHRONOS *)chronos;

	if (chr == NULL || no >= 100) return;
	chr->tic[no] = attHRClockElapsedTime(chr->hr_clock[no]);
}


void chronos_toc(void * chronos, int no)
{
	CHRONOS * chr = (CHRONOS *)chronos;

	if (chr == NULL || no >= 100) return;
	chr->toc[no] = attHRClockElapsedTime(chr->hr_clock[no]);
	chr->dt[no] += chr->toc[no] - chr->tic[no];
}


double chronos_get(void * chronos, int no)
{
	CHRONOS * chr = (CHRONOS *)chronos;

	if (chr == NULL || no >= 100) return 0.0;
	return chr->dt[no];
}
#endif


#ifdef _CUDA_
typedef struct _CUDA_CHRONOS
{
	cudaEvent_t *start, *stop;
	float *dt;
}CUDA_CHRONOS;
#endif

void *cuda_chronos_init()
{
#ifdef _CUDA_
	CUDA_CHRONOS *data = NULL;
	int i;

	data = (CUDA_CHRONOS*)calloc(1, sizeof(CUDA_CHRONOS));
	data->start = (cudaEvent_t*)calloc(100, sizeof(cudaEvent_t));
	data->stop = (cudaEvent_t*)calloc(100, sizeof(cudaEvent_t));
	data->dt = (float*)calloc(100, sizeof(float));
	for (i = 0; i<100; i++)
	{
		cudaEventCreate(&data->start[i]);
		cudaEventCreate(&data->stop[i]);
	}
	return (void*)data;
#else
  return NULL;
#endif
}

void cuda_chronos_tic(void *_chronos, int no)
{
#ifdef _CUDA_
	CUDA_CHRONOS *chronos = (CUDA_CHRONOS*)_chronos;

	if (chronos == NULL) return;
	cudaEventRecord(chronos->start[no], NULL);
#endif
}

void cuda_chronos_toc(void *_chronos, int no)
{
#ifdef _CUDA_
	CUDA_CHRONOS *chronos = (CUDA_CHRONOS*)_chronos;
	float t = 0.0f;

	if (chronos == NULL) return;
	cudaEventRecord(chronos->stop[no], NULL);
	cudaEventSynchronize(chronos->stop[no]);
	cudaEventElapsedTime(&t, chronos->start[no], chronos->stop[no]);
	chronos->dt[no] += t;
#endif
}

float cuda_chronos_get(void *_chronos, int no)
{
#ifdef _CUDA_
	CUDA_CHRONOS *chronos = (CUDA_CHRONOS*)_chronos;
	float t = 0.0f;

	if (chronos == NULL) return 0.0f;
	return chronos->dt[no];
#else
  return 0.0;
#endif
}

void cuda_chronos_reset(void *_chronos, int no)
{
#ifdef _CUDA_
	CUDA_CHRONOS *chronos = (CUDA_CHRONOS*)_chronos;

	if (chronos == NULL) return;
	chronos->dt[no] = 0.0f;
#endif
}

void *cuda_chronos_release(void *_chronos)
{
#ifdef _CUDA_
	CUDA_CHRONOS *chronos = (CUDA_CHRONOS*)_chronos;

	if (chronos == NULL) return NULL;
	FREE(chronos->start);
	FREE(chronos->stop);
	FREE(chronos->dt);
	FREE(chronos)
#endif
	return NULL;
}

#endif