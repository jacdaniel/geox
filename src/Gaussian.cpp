
#include <math.h>

#include "Gaussian.h"

template <typename T> T* Gaussian::smooth1d_(double sigma, int* size)
{
	T* mask = NULL;
	double norm = 0.0, x, den;
	long width = sigmaToSize(sigma), i, width_2;

	width_2 = width / 2;

	den = 2.0 * sigma * sigma;
	mask = new T[width];
	for (i = 0; i < width; i++)
	{
		x = (double)((double)i - width_2);
		mask[i] = (T)exp(-(x * x) / den);
		norm += (double)mask[i];
	}
	for (i = 0; i < width; i++)
		mask[i] /= norm;
	if (size)
		*size = width;
	return mask;
}

template <typename T> T* Gaussian::gradient1d_(double sigma, int* size)
{
	T* mask = NULL;
	double norm = 0.0, x, den;
	long width = sigmaToSize(sigma), i, width_2;

	width_2 = width / 2;

	den = 2.0 * sigma * sigma;
	mask = new T[width];
	for (i = 0; i < width; i++)
	{
		x = (double)((double)i - width_2);
		mask[i] = (T)(-x * exp(-(x * x) / den));
		// norm += fabs(mask[i]);	
		norm += (double)mask[i] * ((double)width / 2 - i);
	}
	for (i = 0; i < width; i++)
		mask[i] /= (T)norm;
	if (size)
		*size = width;
	return mask;
}


long Gaussian::sigmaToSize(double sigma)
{
	return 2 * (long)ceil(3.0 * sigma) + 1;
}


double* Gaussian::smooth1dDouble(double sigma, int* size)
{
	return smooth1d_<double>(sigma, size);
}

float* Gaussian::smooth1dFloat(double sigma, int* size)
{
	return smooth1d_<float>(sigma, size);
}


double* Gaussian::gradient1dDouble(double sigma, int* size)
{
	return gradient1d_<double>(sigma, size);
}


float* Gaussian::gradient1dFloat(double sigma, int* size)
{
	return gradient1d_<float>(sigma, size);
}

void* Gaussian::dataFree(void* data)
{
	delete data;
	return nullptr;
}