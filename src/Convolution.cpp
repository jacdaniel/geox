#include "ExtractInsert.h"
#include "Convolution.h"



template <typename T> int Convolution::convSame1d_(T* f, long size, T* mask, long maskSize, T* g)
{
	long i, j, size2;
	double r;

	size2 = (maskSize-1) / 2;

	for (i = 0; i < size; i++)
		g[i] = (T)0.0;

	for (i = size2; i < size - size2; i++)
	{
		r = 0.0;
		for (j = -size2; j <= size2; j++)
			r += (double)f[i + j] * (double)mask[maskSize - (size2 + j) - 1];
		g[i] = (T)r;
	}
	return 0;
}

template <typename T> int Convolution::convUtil1d_(T* f, long size, T* mask, long maskSize, T* g)
{
	long i, j, size2;
	double r;

	size2 = (maskSize - 1) / 2;

	// for (i = 0; i < size; i++)
	//	g[i] = (T)0.0;

	for (i = size2; i < size - size2; i++)
	{
		r = 0.0;
		for (j = -size2; j <= size2; j++)
			r += (double)f[i + j] * (double)mask[maskSize - (size2 + j) - 1];
		g[i-size2] = (T)r;
	}
	return 0;
}


template <typename T> int Convolution::convSameX3d_(T* f, int* size, T* mask, int maskSize, T* g)
{
	T* in = new T[size[0]];
	T* out = new T[size[0]];

	for (int z = 0; z < size[2]; z++)
	{
		for (int y = 0; y < size[1]; y++)
		{
			ExtractInsert::extractX(f, size, y, z, in);
			Convolution::convSame1d(in, size[0], mask, maskSize, out);
			ExtractInsert::insertX(g, size, y, z, out);
		}
	}
	delete[] in;
	delete[] out;
	return 0;
}

template <typename T> int Convolution::convSameY3d_(T* f, int* size, T* mask, int maskSize, T* g)
{
	T* in = new T[size[1]];
	T* out = new T[size[1]];

	for (int z = 0; z < size[2]; z++)
	{
		for (int x = 0; x < size[0]; x++)
		{
			ExtractInsert::extractY(f, size, x, z, in);
			Convolution::convSame1d(in, size[1], mask, maskSize, out);
			ExtractInsert::insertY(g, size, x, z, out);
		}
	}
	delete[] in;
	delete[] out;
	return 0;
}

template <typename T> int Convolution::convSameZ3d_(T* f, int* size, T* mask, int maskSize, T* g)
{
	if (size[2] == 1) return 0;
	T* in = new T[size[2]];
	T* out = new T[size[2]];

	for (int y = 0; y < size[1]; y++)
	{
		for (int x = 0; x < size[0]; x++)
		{
			ExtractInsert::extractZ(f, size, x, y, in);
			Convolution::convSame1d(in, size[2], mask, maskSize, out);
			ExtractInsert::insertY(g, size, x, y, out);
		}
	}
	delete[] in;
	delete[] out;
	return 0;
}


template <typename T> int Convolution::convSame3d_(T* f, int* size, T* maskx, int maskSizex, T* masky, int maskSizey, T* maskz, int maskSizez, T* g)
{
	convSameX3d_<T>(f, size, maskx, maskSizex, g);
	convSameY3d_<T>(g, size, masky, maskSizey, g);
	convSameZ3d_<T>(g, size, maskz, maskSizez, g);
	return 0;
}


int Convolution::convSame1d(double* f, long size, double* mask, long maskSize, double* g)
{

	return convSame1d_<double>(f, size, mask, maskSize, g);
}

int Convolution::convSame1d(float* f, long size, float* mask, long maskSize, float* g)
{
	return convSame1d_<float>(f, size, mask, maskSize, g);;
}


int Convolution::convSameX3d(double* f, int* size, double* mask, int maskSize, double* g)
{
	return convSameX3d_<double>(f, size, mask, maskSize, g);

}

int Convolution::convSameY3d(double* f, int* size, double* mask, int maskSize, double* g)
{
	return convSameY3d_<double>(f, size, mask, maskSize, g);
}

int Convolution::convSameZ3d(double* f, int* size, double* mask, int maskSize, double* g)
{
	return convSameZ3d_<double>(f, size, mask, maskSize, g);
}

int Convolution::convSameX3d(float* f, int* size, float* mask, int maskSize, float* g)
{
	return convSameX3d_<float>(f, size, mask, maskSize, g);
}

int Convolution::convSameY3d(float* f, int* size, float* mask, int maskSize, float* g)
{
	return convSameY3d_<float>(f, size, mask, maskSize, g);
}

int Convolution::convSameZ3d(float* f, int* size, float* mask, int maskSize, float* g)
{
	return convSameZ3d_<float>(f, size, mask, maskSize, g);
}


int Convolution::convSame3d(double* f, int* size, double* maskx, int maskSizex, double* masky, int maskSizey, double* maskz, int maskSizez, double* g)
{
	return convSame3d_<double>(f, size, maskx, maskSizex, masky, maskSizey, maskz, maskSizez, g);
}

int Convolution::convSame3d(float* f, int* size, float* maskx, int maskSizex, float* masky, int maskSizey, float* maskz, int maskSizez, float* g)
{
	return convSame3d_<float>(f, size, maskx, maskSizex, masky, maskSizey, maskz, maskSizez, g);
}