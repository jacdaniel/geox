// #pragma once

#ifndef __CONVOLUTION__
#define __CONVOLUTION__

class Convolution
{
private:
	template <typename T> static int convSame1d_(T* f, long size, T* mask, long maskSize, T* g);
	template <typename T> static int convUtil1d_(T* f, long size, T* mask, long maskSize, T* g);

	template <typename T> static int convSameX3d_(T* f, int* size, T* mask, int maskSize, T* g);
	template <typename T> static int convSameY3d_(T* f, int* size, T* mask, int maskSize, T* g);
	template <typename T> static int convSameZ3d_(T* f, int* size, T* mask, int maskSize, T* g);
	template <typename T> static int convSame3d_(T* f, int* size, T* maskx, int maskSizex, T* masky, int maskSizey, T* maskz, int maskSizez, T* g);



public:
	static int convSame1d(double* f, long size, double* mask, long maskSize, double* g);
	static int convSame1d(float* f, long size, float* mask, long maskSize, float* g);

	static int convSameX3d(double* f, int* size, double* mask, int maskSize, double* g);
	static int convSameY3d(double* f, int* size, double* mask, int maskSize, double* g);
	static int convSameZ3d(double* f, int* size, double* mask, int maskSize, double* g);

	static int convSameX3d(float* f, int* size, float* mask, int maskSize, float* g);
	static int convSameY3d(float* f, int* size, float* mask, int maskSize, float* g);
	static int convSameZ3d(float* f, int* size, float* mask, int maskSize, float* g);

	static int convSame3d(double* f, int* size, double* maskx, int maskSizex, double* masky, int maskSizey, double* maskz, int maskSizez, double* g);
	static int convSame3d(float* f, int* size, float* maskx, int maskSizex, float* masky, int maskSizey, float* maskz, int maskSizez, float* g);



};


#endif

