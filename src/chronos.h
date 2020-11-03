
#ifndef __CHRONOS__
#define __CHRONOS__


#ifdef __cplusplus
extern "C" {
#endif

#ifdef __ENABLE_CHRONOS__

	void * attHRClockNew();
	void * attHRClockRelease(void * hr_clock);
	double attHRClockDifferenceTime(void * hr_clock_begin, void * hr_clock_end);
	double attHRClockElapsedTime(void * hr_clock_begin);

	void * chronos_init(int nbre, char const *base);
	void * chronos_release(void * chronos);
	void chronos_reset(void * chronos, int no);
	void chronos_tic(void * chronos, int no);
	void chronos_toc(void * chronos, int no);
	double chronos_get(void * chronos, int no);

	void *cuda_chronos_init();
	void cuda_chronos_tic(void *_chronos, int no);
	void cuda_chronos_toc(void *_chronos, int no);
	float cuda_chronos_get(void *_chronos, int no);
	void cuda_chronos_reset(void *_chronos, int no);
	void *cuda_chronos_release(void *_chronos);


#endif

#ifdef __ENABLE_CHRONOS__
#define CHRONOS_INIT chronos_init(100, "milli")
#else
#define CHRONOS_INIT NULL
#endif

#ifdef __ENABLE_CHRONOS__
#define CHRONOS_TIC(chronos, no) chronos_tic(chronos, no)
#else
#define CHRONOS_TIC(chronos, no)
#endif

#ifdef __ENABLE_CHRONOS__
#define CHRONOS_TOC(chronos, no) chronos_toc(chronos, no)
#else
#define CHRONOS_TOC(chronos, no)
#endif

#ifdef __ENABLE_CHRONOS__
#define CHRONOS_GET(chronos, no) chronos_get(chronos, no)
#else
#define CHRONOS_GET(chronos, no) 0.0
#endif

#ifdef __ENABLE_CHRONOS__
#define CHRONOS_RESET(chronos, no) chronos_reset(chronos, no)
#else
#define CHRONOS_RESET(chronos, no)
#endif

#ifdef __ENABLE_CHRONOS__
#define CHRONOS_RELEASE(chronos) chronos_release(chronos)
#else
#define CHRONOS_RELEASE(chronos) NULL
#endif





#ifdef __ENABLE_CHRONOS__
#define CUDA_CHRONOS_INIT cuda_chronos_init();
#else
#define CHRONOS_INIT NULL
#endif

#ifdef __ENABLE_CHRONOS__
#define CUDA_CHRONOS_TIC(chronos, no) cuda_chronos_tic(chronos, no);
#else
#define CHRONOS_TIC(chronos, no)
#endif

#ifdef __ENABLE_CHRONOS__
#define CUDA_CHRONOS_TOC(chronos, no) cuda_chronos_toc(chronos, no);
#else
#define CHRONOS_TOC(chronos, no)
#endif

#ifdef __ENABLE_CHRONOS__
#define CUDA_CHRONOS_GET(chronos, no) cuda_chronos_get(chronos, no);
#else
#define CHRONOS_GET(chronos, no) 0.0
#endif

#ifdef __ENABLE_CHRONOS__
#define CUDA_CHRONOS_RESET(chronos, no) cuda_chronos_reset(chronos, no);
#else
#define CHRONOS_RESET(chronos, no)
#endif

#ifdef __ENABLE_CHRONOS__
#define CUDA_CHRONOS_RELEASE(chronos) cuda_chronos_release(chronos);
#else
#define CHRONOS_RELEASE(chronos) NULL
#endif



#ifdef __cplusplus
}
#endif


#endif