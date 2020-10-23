
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <malloc.h>
#include <OpNabla.h>
#include <ConjugateGradient.h>
#include <DipProcessing.h>
#include <InverseLaplacian.h>
#include <fftw3.h>
#include <stdlib.h>

void debug_data_save(double* data, int* size, char* filename)
{
	long size0 = size[0] * size[1] * size[2];
	FILE* pFile = fopen(filename, "wb");
	fwrite(data, sizeof(double), size0, pFile);
	fclose(pFile);
}

void debug_data_float_save(float* data, int* size, char* filename)
{
	long size0 = size[0] * size[1] * size[2];
	FILE* pFile = fopen(filename, "wb");
	fwrite(data, sizeof(float), size0, pFile);
	fclose(pFile);
}

#define PI 3.141592653589793
void inverseLaplacian(void* in, int *size, void* out)
{
	double* freq = (double*)calloc(size[0] * size[1], sizeof(double));
	fftw_plan p1, p2;

	p1 = fftw_plan_r2r_2d(size[1], size[0], (double*)in, (double*)freq, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);
	p2 = fftw_plan_r2r_2d(size[1], size[0], (double*)freq, (double*)out, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);

	fftw_execute(p1);
	debug_data_save((double*)in, size, "d:\\text.raw");
	debug_data_save(freq, size, "d:\\text.raw");
	for (long y = 0; y < size[1]; y++)
	{
		for (long x = 0; x < size[0]; x++)
		{
			long add = size[0] * y + x;
			double wx = 2.0 * PI * (double)x / (double)(2.0*size[0]);
			double wz = 2.0 * PI * (double)y / (double)(2.0*size[1]);
			double l1 = -4.0 + 2.0 * cos(wx) + 2.0 * cos(wz);
			// freq[add] = l1;
			if (l1 == 0.0) freq[add] = freq[add]; else freq[add] /= l1;
		}
	}
	debug_data_save(freq, size, "d:\\text.raw");
	fftw_execute(p2);
	debug_data_save((double*)out, size, "d:\\text.raw");
	double norm = (double)(size[0] * size[1])*4.0;
	for (long add = 0; add < size[0] * size[1]; add++)
	{
		((double*)out)[add] /= norm;
	}

	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	free(freq);
}

void debug_data_short_save(short* data, int* size, char* filename)
{
	long size0 = size[0] * size[1] * size[2];
	FILE* pFile = fopen(filename, "wb");
	fwrite(data, sizeof(short), size0, pFile);
	fclose(pFile);
}


//
void dipEstimation2D(short* data, int* size, double sigmaG, double sigmaT, short* dip)
{
	// int sizeB[] = { 16, 16, 1 };
	int sizeB[] = { 256, 256, 1 };
	DipProcessing *p = new DipProcessing();
	p->setDataIn(data);
	p->setSize(size);
	p->setBlocSize(sizeB);
	p->setSigmaGradient(sigmaG);
	p->setSigmaTensor(sigmaT);
	p->setDipxXY(dip);
	p->setDipXZ(nullptr);
	p->run();
}

void rhsInit(double* data, int* size, double* dipxy)
{
	int size0 = size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));

	OpNabla::nablatX(dipxy, size, data);
	for (long add = 0; add < size0; add++)
		data[add] *= -dipxy[add];
	OpNabla::nablatY(dipxy, size, tmp); 
	for (long add = 0; add < size0; add++)
		data[add] -= tmp[add];
	free(tmp);
}


/*
void rhsInit(double* data, int* size, short* dipxy)
{
	int size0 = size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));

	debug_data_short_save(dipxy, size, "d:\\text.raw");

	OpNabla::nablatX(dipxy, 1000, size, data);
	debug_data_save(data, size, "d:\\text.raw");



	for (long add = 0; add < size0; add++)
		data[add] *= -(double)dipxy[add]/1000.0;

	debug_data_save(data, size, "d:\\text.raw");



	OpNabla::nablatY(dipxy, 1000, size, tmp);

	debug_data_save(tmp, size, "d:\\text.raw");


	for (long add = 0; add < size0; add++)
		data[add] -= tmp[add];

	debug_data_save(data, size, "d:\\text.raw");
	free(tmp);
}
*/


void rhsInit(double* data, int* size, short* dipxy)
{
	int size0 = size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));
	double* in = (double*)calloc(size0, sizeof(double));

	for (int i = 0; i < size0; i++) in[i] = -(double)dipxy[i]/1000.0 * (double)dipxy[i]/1000.0;

	OpNabla::nablatX(in, size, data);
	debug_data_save(data, size, "d:\\text.raw");

	OpNabla::nablatY(dipxy, 1000, size, tmp);

	debug_data_save(tmp, size, "d:\\text.raw");

	for (long add = 0; add < size0; add++)
		data[add] -= tmp[add];

	debug_data_save(data, size, "d:\\text.raw");
	free(tmp);
	free(in);

}


void opDirect(double* in, double* dipxy, int* size, double epsilon, double* data)
{
	long size0 = (long)size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));
	double* tmp2 = (double*)calloc(size0, sizeof(double));

	OpNabla::nablaX(in, size, data);
	for (long add = 0; add < size0; add++) data[add] *= -dipxy[add];
	OpNabla::nablaY(in, size, tmp2);
	for (long add = 0; add < size0; add++) data[add] -= tmp2[add];

	OpNabla::nablatX(data, size, tmp);
	for (long add = 0; add < size0; add++) tmp[add] *= -dipxy[add];
	OpNabla::nablatY(data, size, tmp2);
	for (long add = 0; add < size0; add++) tmp[add] -= tmp2[add];

	for (long add = 0; add < size0; add++) data[add] = tmp[add];

	OpNabla::nablaX(in, size, tmp);
	OpNabla::nablatX(tmp, size, tmp2);
	for (long add = 0; add < size0; add++) data[add] += epsilon * tmp2[add];

	free(tmp);
	free(tmp2);
}

/*
void opDirect(double* in, short* dipxy, int* size, double epsilon, double* data)
{
	long size0 = (long)size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));
	double* tmp2 = (double*)calloc(size0, sizeof(double));

	OpNabla::nablaX(in, size, data);
	for (long add = 0; add < size0; add++) data[add] *= -(double)dipxy[add]/1000.0;
	OpNabla::nablaY(in, size, tmp2);
	for (long add = 0; add < size0; add++) data[add] -= tmp2[add];

	OpNabla::nablatX(data, size, tmp);
	for (long add = 0; add < size0; add++) tmp[add] *= -(double)dipxy[add]/1000.0;
	OpNabla::nablatY(data, size, tmp2);
	for (long add = 0; add < size0; add++) tmp[add] -= tmp2[add];

	for (long add = 0; add < size0; add++) data[add] = tmp[add];

	OpNabla::nablaX(in, size, tmp);
	OpNabla::nablatX(tmp, size, tmp2);
	for (long add = 0; add < size0; add++) data[add] += epsilon * tmp2[add];

	free(tmp);
	free(tmp2);
}
*/

void opDirect(double* in, short* dipxy, int* size, double epsilon, double* data)
{
	long size0 = (long)size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));
	double* tmp2 = (double*)calloc(size0, sizeof(double));
	double* data2 = (double*)calloc(size0, sizeof(double));

	OpNabla::nablaX(in, size, data);
	for (long add = 0; add < size0; add++) data[add] *= -(double)dipxy[add] / 1000.0;
	OpNabla::nablaY(in, size, tmp2);
	for (long add = 0; add < size0; add++) data[add] -= tmp2[add];

	for (int i = 0; i < size0; i++) data2[i] = -data[i] * (double)dipxy[i]/1000.0;

	OpNabla::nablatX(data2, size, tmp);
	// for (long add = 0; add < size0; add++) tmp[add] *= -(double)dipxy[add] / 1000.0;
	OpNabla::nablatY(data, size, tmp2);
	for (long add = 0; add < size0; add++) tmp[add] -= tmp2[add];

	for (long add = 0; add < size0; add++) data[add] = tmp[add];

	OpNabla::nablaX(in, size, tmp);
	OpNabla::nablatX(tmp, size, tmp2);
	for (long add = 0; add < size0; add++) data[add] += epsilon * tmp2[add];

	free(tmp);
	free(tmp2);
	free(data2);
}


class myCallBack2 : public ConjugateGradient::CALLBACK
{
private:
	class PARAM
	{
	public:
		InverseLaplacian* lap;
	};
public:
	myCallBack2();
	int *size, size0;
	short* dipxy;
	double epsilon;
	void setSize(int *size);
	void CallBack(void* in, void* out);
	void Preconditionner(void* in, void* out);
	void setDipxy(short* dipxy);
	void setEpsilon(double epsilon);
private:
	PARAM* param;
};

myCallBack2::myCallBack2()
{
	param = nullptr;
}

void myCallBack2::setSize(int *size)
{
	this->size = size;
}

void myCallBack2::setDipxy(short* dipxy)
{
	this->dipxy = dipxy;
}

void myCallBack2::setEpsilon(double epsilon)
{
	this->epsilon = epsilon;
}


void myCallBack2::CallBack(void* in, void* out)
{
	opDirect((double*)in, dipxy, size, epsilon, (double*)out);
}

void myCallBack2::Preconditionner(void* in, void* out)
{
	// inverseLaplacian(in, size, out);
	if (this->param == nullptr)
	{
		this->param = new myCallBack2::PARAM();
		param->lap = new InverseLaplacian();
		param->lap->setdataIn(in);
		param->lap->setDataOut(out);
		param->lap->setSize(this->size[0], this->size[1], this->size[2]);
	}
	param->lap->run();
	// memcpy(out, in, (size_t)size[0] * size[1] * size[2] * sizeof(double));
}


#define YSWAP(x) ((x & 0xff) << 8) | ((x & 0xff00) >> 8)

int main_geotimeConstraint(int argc, char** argv)
{
	int size[] = { 350, 2000, 1 };
	long size0 = (long)size[0] * size[1] * size[2];

	short* dipxy = (short*)calloc(size0, sizeof(short));
	for (int y = 0; y < size[1]; y++)
		for (int x = 0; x < size[0]; x++)
		{
			dipxy[y * size[0] + x] = 200;// (short)(.5 + (double)y / size[1] * 1.5 * 1000.0);
		}

	short* seismic = new short[static_cast<unsigned __int64>(size[0]) * size[1]];
	for (int i = 0; i < (long)size[0] * size[1]; i++) seismic[i] = i;

	char* filename = "D:\\TOTAL\\PLI\\seismic.xt";
	FILE* pf = nullptr;
	pf = fopen(filename, "rb");
	int ret = fseek(pf, (long)5120+200*350*2000*2, SEEK_SET);
	ret = fread(seismic, sizeof(short), size[0]*size[1], pf);
	fclose(pf);
	for (int i = 0; i < size[0]*size[1]; i++)
		seismic[i] = YSWAP(seismic[i]);


	dipEstimation2D(seismic, size, 1.0, 1.5, dipxy);
	// return 0;

	double* tau = (double*)calloc(size0, sizeof(double));	
	double* rhs = (double*)calloc(size0, sizeof(double));
	rhsInit(rhs, size, dipxy);
/*
	inverseLaplacian(rhs, size, tau);
	for (long y = 0; y < size[1]; y++)
	{
		for (long x = 0; x < size[0]; x++)
		{
			tau[size[0] * y + x] -= x;
		}
	}
	*/


	debug_data_save(rhs, size, "d:\\text.raw");
	debug_data_save(tau, size, "d:\\text.raw");

	myCallBack2* c = new myCallBack2();
	c->setSize(size);
	c->setDipxy(dipxy);
	c->setEpsilon(0.005);

	ConjugateGradient* p = new ConjugateGradient();
	p->setCallback(c, nullptr);
	p->setSize(size0);
	p->setRhs(rhs);
	p->setX(tau);
	p->setNbiter(20);
	p->run();

	delete p;

	for (long y = 0; y < size[1]; y++)
	{
		for (long x = 0; x < size[0]; x++)
		{
			tau[size[0] * y + x] += x;
		}
	}
	debug_data_save(tau, size, "d:\\text.raw");

	fprintf(stderr, "%s ok\n", __FILE__);
	return 0;
}