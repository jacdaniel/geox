
#ifndef __FFTW_UTIL__
#define __FFTW_UTIL__

#define FFTW_DESTROY_PLAN(_plan) {\
	if ( _plan != NULL ) fftw_destroy_plan(_plan); _plan = NULL; }

#define FFTW_FREE(_ptr) { if ( _ptr != NULL ) fftw_free(_ptr); _ptr = NULL; }

#define INIT_PLAN_DST_I_1D(size, in, out) fftw_plan_r2r_1d(size, in, out, FFTW_RODFT00, FFTW_ESTIMATE)
#define INIT_PLAN_DST_I_2D(size2, size1, in, out) fftw_plan_r2r_2d(size2, size1, in, out, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE)

#define INIT_PLAN_DCT_II_1D(size, in, out) fftw_plan_r2r_1d(size, in, out, FFTW_REDFT10, FFTW_ESTIMATE)
#define INIT_PLAN_DCT_II_2D(size2, size1, in, out) fftw_plan_r2r_2d(size2, size1, in, out, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE)

#define INIT_PLAN_DCT_III_1D(size, in, out) fftw_plan_r2r_1d(size, in, out, FFTW_REDFT01, FFTW_ESTIMATE)
#define INIT_PLAN_DCT_III_2D(size2, size1, in, out) fftw_plan_r2r_2d(size2, size1, in, out, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)

#define INIT_PLAN_IDST_I_1D INIT_PLAN_DST_I_1D
#define INIT_PLAN_IDST_I_2D INIT_PLAN_DST_I_2D

#define INIT_PLAN_IDCT_II_1D INIT_PLAN_DCT_III_1D
#define INIT_PLAN_IDCT_II_2D INIT_PLAN_DCT_III_2D

#endif
