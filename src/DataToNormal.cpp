
#include <Convolution.h>
#include <Gaussian.h>
#include <eigen.h>

#include <DataToNormal.h>

DataToNormal::DataToNormal()
{
	this->param = nullptr;
	this->size = nullptr;
	setSigmaGradient(1.0f);
	setSigmaTensor(1.5f);	
}

DataToNormal::~DataToNormal()
{
	if (this->param != nullptr) delete this->param;
}

int DataToNormal::setSigmaGradient(float sigma)
{
	this->sigmaG = sigma;
	return 0;
}

int DataToNormal::setSigmaTensor(float sigma)
{
	this->sigmaT = sigma;
	return 0;
}

int DataToNormal::setSize(int* size)
{
	this->size = size;
	return 0;
}



int DataToNormal::dataToDip(float* data, int* size, short* dip)
{
	long size0 = this->size[0] * this->size[1] * this->size[2];
	if (this->param == nullptr)
	{
		this->param = new DataToNormal::PARAM();
		param->gx = new float[size0];
		param->gy = new float[size0];
		param->Txx = new float[size0];
		param->Txy = new float[size0];
		param->Tyy = new float[size0];
		param->smooth = Gaussian::smooth1dFloat(this->sigmaG, &param->sizeMaskD);
		param->deriv = Gaussian::gradient1dFloat(this->sigmaG, &param->sizeMaskD);
		param->smooth2 = Gaussian::smooth1dFloat(this->sigmaT, &param->sizeMask2);
	}

	/*
	long size0 = size[0] * size[1];
	float* gx = new float[size0];
	float* gy = new float[size0];
	float* Txx = new float[size0];
	float* Txy = new float[size0];
	float* Tyy = new float[size0];

	int sizeMaskD = 0, sizeMask2 = 0;
	float* smooth = Gaussian::smooth1dFloat(1.0, &sizeMaskD);
	float* deriv = Gaussian::gradient1dFloat(1.0, &sizeMaskD);
	float* smooth2 = Gaussian::smooth1dFloat(1.5, &sizeMask2);
	*/

	Convolution::convSame3d(data, this->size, param->deriv, param->sizeMaskD, param->smooth, param->sizeMaskD, param->smooth, param->sizeMaskD, param->gx);
	Convolution::convSame3d(data, this->size, param->smooth, param->sizeMaskD, param->deriv, param->sizeMaskD, param->smooth, param->sizeMaskD, param->gy);

	for (long add = 0; add < size0; add++) param->Txx[add] = param->gx[add] * param->gx[add];
	for (long add = 0; add < size0; add++) param->Txy[add] = param->gx[add] * param->gy[add];
	for (long add = 0; add < size0; add++) param->Tyy[add] = param->gy[add] * param->gy[add];

	Convolution::convSameX3d(param->Txx, size, param->smooth2, param->sizeMask2, param->Txx); Convolution::convSameY3d(param->Txx, size, param->smooth2, param->sizeMask2, param->Txx);
	Convolution::convSameX3d(param->Txy, size, param->smooth2, param->sizeMask2, param->Txy); Convolution::convSameY3d(param->Txy, size, param->smooth2, param->sizeMask2, param->Txy);
	Convolution::convSameX3d(param->Tyy, size, param->smooth2, param->sizeMask2, param->Tyy); Convolution::convSameY3d(param->Tyy, size, param->smooth2, param->sizeMask2, param->Tyy);

	float ux, uy;
	for (int add = 0; add < size0; add++)
	{
		Eigen::PrincipalVector2X2(param->Tyy[add], param->Txy[add], param->Txx[add], &ux, &uy);
		double dip_tmp = -ux / uy;
		if (dip_tmp > 2.0) dip_tmp = 2.0;
		if (dip_tmp < -2.0) dip_tmp = -2.0;
		dip[add] = (short)(1000.0 * dip_tmp);
	}
	
	return 0;
}