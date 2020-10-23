// #pragma once

#ifndef __INVERSELAPLACIAN__
#define __INVERSELAPLACIAN__

#include <fftw3.h>

class InverseLaplacian 
{
	class PARAM
	{
	public:
		long size0;
		fftw_plan planDirect, planInverse;
		void* freq, *mask;
		double fftw_norm;
	};
public:
	InverseLaplacian();
	~InverseLaplacian();
	void setSize(int dimx, int dimy, int dimz);
	void setdataIn(void* dataIn);
	void setDataOut(void* dataOut);
	void run();

private:
	void paramInit();
	PARAM* param;
	int dimx, dimy, dimz;
	void* dataIn, *dataOut;


};




#endif
