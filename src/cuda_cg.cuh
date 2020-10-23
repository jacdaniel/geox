#ifndef __CUDA_CG__
#define __CUDA_CG__

void* cuda_cg_init();

void* cuda_cg_release(void* _p);

void cuda_cg_set_callback(void* _p, void* callback, void* data);

void cuda_cg_set_preconditionner(void* _p, void* f, void* data);

void cuda_cg_set_size(void* _p, int size);

void cuda_cg_set_rhs(void* _p, void* rhs);

void cuda_cg_set_x(void* _p, void* x);

void cuda_cg_set_nbiter(void* _p, int nbiter);

int cuda_cg_run(void* _p);



#endif#pragma once
