
#include <stdio.h>

#include "OpNabla.h"

template <typename T> int  OpNabla::nablaX(T* f, int* size, T* g)
{
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			for (long x = 0; x < size[0] - 1; x++)
			{
				g[(long)size[0] * size[1] * z + (long)size[0] * y + x] = f[(long)size[0] * size[1] * (long)z + size[0] * y + x + 1] - f[(long)size[0] * size[1] * z + size[0] * y + x];
			}
			long x = size[0] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)0.0;
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablaY(T* f, int* size, T* g)
{
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long x = 0; x < size[0]; x++)
		{
			for (long y = 0; y < size[1] - 1; y++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * z + size[0] * (y + 1) + x] - f[(long)size[0] * size[1] * z + size[0] * y + x];
			}
			long y = size[1] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)0.0;
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablaZ(T* f, int* size, T* g)
{
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long x = 0; x < size[0]; x++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			for (long z = 0; z < size[2] - 1; z++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * (z + 1) + size[0] * y + x] - f[(long)size[0] * size[1] * z + size[0] * y + x];
			}
			long z = size[2] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)0.0;
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablatX(T* f, int* size, T* g)
{
	long x = 0;
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			x = 0;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = -f[(long)size[0] * size[1] * z + size[0] * y + x];
			for (x = 1; x < size[0] - 1; x++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * z + size[0] * y + x - 1] - f[(long)size[0] * size[1] * z + size[0] * y + x];
			}
			x = size[0] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * z + size[0] * y + x - 1];
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablatY(T* f, int* size, T* g)
{
	long y = 0;
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long x = 0; x < size[0]; x++)
		{
			y = 0;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = -f[(long)size[0] * size[1] * z + size[0] * y + x];
			for (y = 1; y < size[1] - 1; y++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * z + size[0] * (y - 1) + x] - f[(long)size[0] * size[1] * z + size[0] * y + x];
			}
			y = size[1] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * z + size[0] * (y - 1) + x];
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablatZ(T* f, int* size, T* g)
{
	long z = 0;
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long x = 0; x < size[0]; x++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			z = 0;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = -f[(long)size[0] * size[1] * z + size[0] * y + x];
			for (z = 1; z < size[2] - 1; z++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * (z - 1) + size[0] * y + x] - f[(long)size[0] * size[1] * z + size[0] * y + x];
			}
			z = size[2] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = f[(long)size[0] * size[1] * (z - 1) + size[0] * y + x];
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablaX(short* f, double norm, int* size, T* g)
{
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			for (long x = 0; x < size[0] - 1; x++)
			{
				g[(long)size[0] * size[1] * z + (long)size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * (long)z + size[0] * y + x + 1] - f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			}
			long x = size[0] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)0.0;
		}
	}
	return 0;
}





template <typename T> int OpNabla::nablaY(short* f, double norm, int* size, T* g)
{
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long x = 0; x < size[0]; x++)
		{
			for (long y = 0; y < size[1] - 1; y++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * z + size[0] * (y + 1) + x] - f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			}
			long y = size[1] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = 0.0f;
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablaZ(short* f, double norm, int* size, T* g)
{
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long x = 0; x < size[0]; x++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			for (long z = 0; z < size[2] - 1; z++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * (z + 1) + size[0] * y + x] - f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			}
			long z = size[2] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)0.0;
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablatX(short* f, double norm, int* size, T* g)
{
	long x = 0;
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			x = 0;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)(-(double)f[(long)size[0] * size[1] * z + size[0] * y + x]/norm);
			for (x = 1; x < size[0] - 1; x++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * z + size[0] * y + x - 1] - f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			}
			x = size[0] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * z + size[0] * y + x - 1])/norm);
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablatY(short* f, double norm, int* size, T* g)
{
	long y = 0;
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long z = 0; z < size[2]; z++)
	{
		for (long x = 0; x < size[0]; x++)
		{
			y = 0;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = -(T)((double)(f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			for (y = 1; y < size[1] - 1; y++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * z + size[0] * (y - 1) + x] - f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			}
			y = size[1] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * z + size[0] * (y - 1) + x])/norm);
		}
	}
	return 0;
}

template <typename T> int OpNabla::nablatZ(short* f, double norm, int* size, T* g)
{
	long z = 0;
	if (f == NULL || size == NULL || g == NULL) return 1;
	for (long x = 0; x < size[0]; x++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			z = 0;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = -(T)((double)(f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			for (z = 1; z < size[2] - 1; z++)
			{
				g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * (z - 1) + size[0] * y + x] - f[(long)size[0] * size[1] * z + size[0] * y + x])/norm);
			}
			z = size[2] - 1;
			g[(long)size[0] * size[1] * z + size[0] * y + x] = (T)((double)(f[(long)size[0] * size[1] * (z - 1) + size[0] * y + x])/norm);
		}
	}
	return 0;
}



int OpNabla::nablaX(float* f, int* size, float* g)
{
	return nablaX<float>(f, size, g);
}

int OpNabla::nablaY(float* f, int* size, float* g)
{
	return nablaY<float>(f, size, g);
}

int OpNabla::nablaZ(float* f, int* size, float* g)
{
	return nablaZ<float>(f, size, g);
}



int OpNabla::nablaX(double* f, int* size, double* g)
{
	return nablaX<double>(f, size, g);
}

int OpNabla::nablaY(double* f, int* size, double* g)
{
	return nablaY<double>(f, size, g);
}

int OpNabla::nablaZ(double* f, int* size, double* g)
{
	return nablaZ<double>(f, size, g);
}



int OpNabla::nablatX(float* f, int* size, float* g)
{
	return nablatX<float>(f, size,  g);
}

int OpNabla::nablatY(float* f, int* size, float* g)
{
	return nablatY<float>(f, size, g);
}

int OpNabla::nablatZ(float* f, int* size, float* g)
{
	return nablatZ<float>(f, size, g);
}


int OpNabla::nablatX(double* f, int* size, double* g)
{
	return nablatX<double>(f, size, g);
}

int OpNabla::nablatY(double* f, int* size, double* g)
{
	return nablatY<double>(f, size, g);
}

int OpNabla::nablatZ(double* f, int* size, double* g)
{
	return nablatZ<double>(f, size, g);
}



int OpNabla::nablaX(short* f, double norm, int* size, double* g)
{
	return nablaX<double>(f, norm, size, g);
}
int OpNabla::nablaY(short* f, double norm, int* size, double* g)
{
	return nablaY<double>(f, norm, size, g);
}

int OpNabla::nablaZ(short* f, double norm, int* size, double* g)
{
	return nablaZ<double>(f, norm, size, g);
}

int OpNabla::nablatX(short* f, double norm, int* size, double* g)
{
	return nablatX<double>(f, norm, size, g);
}

int OpNabla::nablatY(short* f, double norm, int* size, double* g)
{
	return nablatY<double>(f, norm, size, g);
}

int OpNabla::nablatZ(short* f, double norm, int* size, double* g)
{
	return nablatZ<double>(f, norm, size, g);
}