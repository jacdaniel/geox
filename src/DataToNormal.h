// #pragma once

#ifndef __DATATONORMAL__
#define __DATATONORMAL__


class DataToNormal
{
private:
	class PARAM
	{
	public:
		long size0;
		float* gx, * gy, * Txx, * Txy, * Tyy, *smooth, *deriv, *smooth2;
		int sizeMaskD, sizeMask2;
	};
	float sigmaG, sigmaT;
	int* size;
	PARAM* param;
public:
	DataToNormal();
	~DataToNormal();
	int setSigmaGradient(float sigma);
	int setSigmaTensor(float sigma);
	int setSize(int* size);

	int dataToDip(float* data, int* size, short* dip);
};




#endif
