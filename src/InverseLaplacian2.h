// #pragma once


// #pragma once

#ifndef __INVERSELAPLACIAN2__
#define __INVERSELAPLACIAN2__

#include <fftw3.h>

class InverseLaplacian2
{
public:
	enum DIMS {_1D, _2D, _3D};
	class PARAM
	{
	public:
		long size0, dim;
		fftw_plan planDirect[3], planInverse[3];
		void* in, * out;
		void** arrayCos;
		double fftw_norm[3], fftw_norm0;
	};
public:
	InverseLaplacian2();
	~InverseLaplacian2();
	void setSize(int* size);
	void setdataIn(void* dataIn);
	void setDataOut(void* dataOut);
	void run();

private:
	void paramInit();
	PARAM* param;
	int* size;
	void* dataIn, * dataOut;
};




#endif

