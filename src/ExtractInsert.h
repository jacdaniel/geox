//#pragma once

#ifndef __EXTRACTINSERT__
#define __EXTRACTINSERT__


class ExtractInsert
{
public:
	template <typename T1, typename T2> static int extractX_(T1* f, int* size, int y0, int z0, T2* g);
	template <typename T1, typename T2> static int extractY_(T1* f, int* size, int x0, int z0, T2* g);
	template <typename T1, typename T2> static int extractZ_(T1* f, int* size, int x0, int y0, T2* g);

	template <typename T1, typename T2> static int insertX_(T1* f, int* size, int y0, int z0, T2* g);
	template <typename T1, typename T2> static int insertY_(T1* f, int* size, int x0, int z0, T2* g);
	template <typename T1, typename T2> static int insertZ_(T1* f, int* size, int x0, int y0, T2* g);

	static int extractX(double* f, int* size, int y0, int z0, double* g);
	static int extractY(double* f, int* size, int x0, int z0, double* g);
	static int extractZ(double* f, int* size, int x0, int y0, double* g);

	static int extractX(float* f, int* size, int y0, int z0, float* g);
	static int extractY(float* f, int* size, int x0, int z0, float* g);
	static int extractZ(float* f, int* size, int x0, int y0, float* g);

	static int insertX(double* f, int* size, int y0, int z0, double* g);
	static int insertY(double* f, int* size, int x0, int z0, double* g);
	static int insertZ(double* f, int* size, int x0, int y0, double* g);

	static int insertX(float* f, int* size, int y0, int z0, float* g);
	static int insertY(float* f, int* size, int x0, int z0, float* g);
	static int insertZ(float* f, int* size, int x0, int y0, float* g);
};



#endif
