#ifndef __CUDA_OP__
#define __CUDA_OP__


void cuda_y_equals_a_plus_b_double(long size, double* a, double* b, double* y);

void cuda_y_plus_equals_alpha_x_double(long size, double* y, double alpha, double* x);

void cuda_y_equal_a_times_b(long size, double* a, char* b, double* y);

void cuda_y_minus_equals_x(long size, double* y, double x);

void cuda_copy_double_float(long size, double*d_in, float *d_out);

void cuda_copy_float_double(long size, float* d_in, double* d_out);

void cuda_y_minus_equals_yref(long size, double *y, long add);

double cuda_dot0(long size, double* x, double* y);




#endif