
#include <math.h>
#include <malloc.h>
#include <InverseLaplacian.h>


InverseLaplacian::InverseLaplacian()
{
	this->param = nullptr;
}


InverseLaplacian::~InverseLaplacian()
{
	fftw_destroy_plan(param->planDirect);
	fftw_destroy_plan(param->planInverse);
	fftw_free(param->freq);
	free(param->mask);
}

void InverseLaplacian::setSize(int dimx, int dimy, int dimz)
{
	this->dimx = dimx;
	this->dimy = dimy;
	this->dimz = dimz;
}

void InverseLaplacian::setdataIn(void* dataIn)
{
	this->dataIn = dataIn;
}

void InverseLaplacian::setDataOut(void* dataOut)
{
	this->dataOut = dataOut;
}

#define PI 3.141592653589793
void InverseLaplacian::paramInit()
{
	this->param = new PARAM();
	this->param->size0 = (long)dimx * dimy * dimz;

	this->param->freq = fftw_alloc_real(this->param->size0 * sizeof(double));
	this->param->planDirect = fftw_plan_r2r_2d(dimy, dimx, (double*)this->dataIn, (double*)param->freq, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
	this->param->planInverse = fftw_plan_r2r_2d(dimy, dimx, (double*)param->freq, (double*)this->dataOut, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);

	param->mask = calloc(param->size0, sizeof(double));
	for (long y = 0; y < dimy; y++)
	{
		for (long x = 0; x < dimx; x++)
		{
			long add = dimx * y + x;
			double wx = 2.0 * PI * (double)x / (double)(2.0 * dimx);
			double wz = 2.0 * PI * (double)y / (double)(2.0 * dimy);
			double l1 = -4.0 + 2.0 * cos(wx) + 2.0 * cos(wz);
			// freq[add] = l1;
			if (l1 == 0.0) ((double*)param->mask)[add] = 1.0; else ((double*)param->mask)[add] = l1;
		}
	}
	param->fftw_norm = (double)(param->size0) * 4.0;
}

void InverseLaplacian::run()
{
	if (this->param == nullptr)
	{
		paramInit();
	}
	fftw_execute(param->planDirect);
	for (long add = 0; add < param->size0; add++)
	{
		((double*)param->freq)[add] /= ((double*)param->mask)[add];
	}
	fftw_execute(param->planInverse);
	for (long add = 0; add < param->size0; add++)
	{
		((double*)this->dataOut)[add] /= param->fftw_norm;
	}
}