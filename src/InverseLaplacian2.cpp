// #pragma once

#include <malloc.h>
#include <math.h>
#include <algorithm>

#include <ExtractInsert.h>
#include <InverseLaplacian2.h>

InverseLaplacian2::InverseLaplacian2()
{
	this->param = nullptr;
}

InverseLaplacian2::~InverseLaplacian2()
{

}

void InverseLaplacian2::setSize(int *size)
{
	this->size = size;
}


void InverseLaplacian2::setdataIn(void* dataIn)
{
	this->dataIn = dataIn;
}

void InverseLaplacian2::setDataOut(void* dataOut)
{
	this->dataOut = dataOut;
}

#define PI 3.141592653589793

void InverseLaplacian2::paramInit()
{
	int smax = std::max(size[0], size[1]);
	smax = std::max(smax, size[2]);

	this->param = new InverseLaplacian2::PARAM();
	if (size[2] == 1) param->dim = DIMS::_2D; else param->dim = _3D;
	this->param->size0 = (long)size[0] * size[1] * size[2];

	this->param->in = (void*)fftw_alloc_real(smax*sizeof(double));
	this->param->out = (void**)fftw_alloc_real(smax*sizeof(double));
	this->param->arrayCos = (void**)calloc(3, sizeof(void*));

	for (int n = 0; n < 3; n++)
	{
		this->param->planDirect[n] = fftw_plan_r2r_1d(size[n], (double*)param->in, (double*)param->out, FFTW_REDFT10, FFTW_ESTIMATE);
		this->param->planInverse[n] = fftw_plan_r2r_1d(size[n], (double*)param->in, (double*)param->out, FFTW_REDFT01, FFTW_ESTIMATE);

		param->fftw_norm[n] = (double)(size[n]) * 4.0;

		this->param->arrayCos[n] = (void*)calloc(size[n], sizeof(double));
		for (int add = 0; add < size[n]; add++)
		{
			double wx = 2.0 * PI * (double)add / (double)(2.0 * size[n]);
			((double**)this->param->arrayCos)[n][add] = 2.0 * cos(wx);
		}
	}
	if (param->dim == _2D)
		param->fftw_norm0 = (double)(param->size0) * 4.0;
	else
		param->fftw_norm0 = (double)(param->size0) * 8.0;
}

void InverseLaplacian2::run()
{
	if (this->param == nullptr)
	{
		paramInit();
	}

	long z = 0;
	for (long y = 0; y < size[1]; y++)
	{
		ExtractInsert::extractX((double*)dataIn, size, y, z, (double*)param->in);
		fftw_execute(param->planDirect[0]);
		ExtractInsert::insertX((double*)dataOut, size, y, z, (double*)param->out);
	}

	for (long x = 0; x < size[0]; x++)
	{
		ExtractInsert::extractY((double*)dataOut, size, x, z, (double*)param->in);
		fftw_execute(param->planDirect[1]);
		ExtractInsert::insertY((double*)dataOut, size, x, z, (double*)param->out);
	}

	double c0 = -4.0;
	if (param->dim == _3D) c0 = -6.0;

	for (long z=0; z<size[2]; z++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			for (long x = 0; x < size[0]; x++)
			{
				long add = size[0] * y + x;
				double l1 = -6.0 + ((double**)this->param->arrayCos)[0][x] + ((double**)this->param->arrayCos)[1][y] + ((double**)this->param->arrayCos)[2][z];
				if (l1 == 0.0)
				{
					((double*)dataOut)[add] = 0.0;
				}
				else
				{
					((double*)dataOut)[add] /= l1;
				}
			}
		}
	}

	for (long y = 0; y < size[1]; y++)
	{
		ExtractInsert::extractX((double*)dataOut, size, y, z, (double*)param->in);
		fftw_execute(param->planInverse[0]);
		ExtractInsert::insertX((double*)dataOut, size, y, z, (double*)param->out);
	}

	for (long x = 0; x < size[0]; x++)
	{
		ExtractInsert::extractY((double*)dataOut, size, x, z, (double*)param->in);
		fftw_execute(param->planInverse[1]);
		ExtractInsert::insertY((double*)dataOut, size, x, z, (double*)param->out);
	}

	for (long add = 0; add < param->size0; add++)
	{
		((double*)this->dataOut)[add] /= param->fftw_norm0;
	}
}