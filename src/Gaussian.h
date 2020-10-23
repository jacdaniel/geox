// #pragma once

#ifndef __GAUSSIAN__
#define __GAUSSIAN__


class Gaussian
{
public:
	static long sigmaToSize(double sigma);
	template <typename T> static T* smooth1d_(double sigma, int* size);
	template <typename T> static T* gradient1d_(double sigma, int* size);

	static double* smooth1dDouble(double sigma, int* size);
	static float* smooth1dFloat(double sigma, int* size);
	static double* gradient1dDouble(double sigma, int* size);
	static float* gradient1dFloat(double sigma, int* size);

	static void* dataFree(void* data);
};


#endif