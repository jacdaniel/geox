// #pragma once

#ifndef __DIPPROCESSING__
#define __DIPPROCESSING__

class DipProcessing
{
private:
	void* dataIn, * dipXY, * dipXZ;
	int* sizeIn, *blockSize;
	double sigmaGrad, sigmaTens;

public:
	void setDataIn(void* in);
	void setSize(int* size);
	void setBlocSize(int* size);
	void setSigmaGradient(double sigma);
	void setSigmaTensor(double sigma);
	void setDipxXY(void* data);
	void setDipXZ(void* data);
	void run();
};



#endif