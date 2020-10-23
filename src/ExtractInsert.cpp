#include "ExtractInsert.h"

template <typename T1, typename T2> int ExtractInsert::extractX_(T1* f, int* size, int y0, int z0, T2* g)
{
	for (long x = 0; x < size[0]; x++)
	{
		g[x] = (T2)f[(long)size[0] * size[1] * z0 + size[0] * y0 + x];
	}

	return 0;
}

template <typename T1, typename T2> int ExtractInsert::extractY_(T1* f, int* size, int x0, int z0, T2* g)
{
	for (long y = 0; y < size[1]; y++)
	{
		g[y] = (T2)f[(long)size[0] * size[1] * z0 + size[0] * y + x0];
	}

	return 0;
}

template <typename T1, typename T2> int ExtractInsert::extractZ_(T1* f, int* size, int x0, int y0, T2* g)
{
	for (long z = 0; z < size[2]; z++)
	{
		g[z] = (T2)f[(long)size[0] * size[1] * z + size[0] * y0 + x0];
	}

	return 0;
}





template <typename T1, typename T2> int ExtractInsert::insertX_(T1* f, int* size, int y0, int z0, T2* g)
{
	for (long x = 0; x < size[0]; x++)
	{
		f[(long)size[0] * size[1] * z0 + size[0] * y0 + x] = (T1)g[x];
	}

	return 0;
}

template <typename T1, typename T2> int ExtractInsert::insertY_(T1* f, int* size, int x0, int z0, T2* g)
{
	for (long y = 0; y < size[1]; y++)
	{
		f[(long)size[0] * size[1] * z0 + size[0] * y + x0] = (T1)g[y];
	}

	return 0;
}

template <typename T1, typename T2> int ExtractInsert::insertZ_(T1* f, int* size, int x0, int y0, T2* g)
{
	for (long z = 0; z < size[2]; z++)
	{
		f[(long)size[0] * size[1] * z + size[0] * y0 + x0] = (T1)g[z];
	}

	return 0;
}


int ExtractInsert::extractX(double* f, int* size, int y0, int z0, double* g)
{
	return extractX_<double, double>(f, size, y0, z0, g);
}

int ExtractInsert::extractY(double* f, int* size, int x0, int z0, double* g)
{
	return extractY_<double, double>(f, size, x0, z0, g);
}
int ExtractInsert::extractZ(double* f, int* size, int x0, int y0, double* g)
{
	return extractZ_<double, double>(f, size, x0, y0, g);
}

int ExtractInsert::extractX(float* f, int* size, int y0, int z0, float* g)
{
	return extractX_<float, float>(f, size, y0, z0, g);
}

int ExtractInsert::extractY(float* f, int* size, int x0, int z0, float* g)
{
	return extractY_<float, float>(f, size, x0, z0, g);
}

int ExtractInsert::extractZ(float* f, int* size, int x0, int y0, float* g)
{
	return extractZ_<float, float>(f, size, x0, y0, g);
}



int ExtractInsert::insertX(double* f, int* size, int y0, int z0, double* g)
{
	return insertX_<double, double>(f, size, y0, z0, g);
}

int ExtractInsert::insertY(double* f, int* size, int x0, int z0, double* g)
{
	return insertY_<double, double>(f, size, x0, z0, g);
}
int ExtractInsert::insertZ(double* f, int* size, int x0, int y0, double* g)
{
	return insertZ_<double, double>(f, size, x0, y0, g);
}

int ExtractInsert::insertX(float* f, int* size, int y0, int z0, float* g)
{
	return insertX_<float, float>(f, size, y0, z0, g);
}

int ExtractInsert::insertY(float* f, int* size, int x0, int z0, float* g)
{
	return insertY_<float, float>(f, size, x0, z0, g);
}

int ExtractInsert::insertZ(float* f, int* size, int x0, int y0, float* g)
{
	return insertZ_<float, float>(f, size, x0, y0, g);
}