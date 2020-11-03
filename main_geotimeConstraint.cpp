
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <malloc.h>

#include <iostream>
#include <chrono>
#include <chronos.h>
#include <OpNabla.h>
#include <ConjugateGradient.h>
#include <DipProcessing.h>
#include <InverseLaplacian.h>
#include <InverseLaplacian2.h>
#include <fftw3.h>
#include <xt_file.h>
#include <stdlib.h>

void* CHRONOS;
#define CHRONOS_OPDIRECT 1
#define CHRONOS_PRECOND  2


typedef struct _CONSTRAINT
{
	int* no, * x0, * y0, * z0, nbre;
}CONSTRAINT;

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
	// int sizeB[] = { 512, 512, 1 };

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

void dipEstimationSynthetic3D(int* size, short* dipxy, short* dipxz)
{
	long size0 = size[0] * size[1] * size[2];

	for (long add = 0; add < size0; add++)
	{
		dipxy[add] = 100;
		dipxz[add] = 200;
	}

	for (long z = 0; z < size[2]; z++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			for (long x = 0; x < size[0]; x++)
			{
				long add = size[0] * size[1] * z + size[0] * y + x;
				dipxy[add] = 300.0 * (double)x / size[0] + 100;
				dipxz[add] = 200.0 * (double)x / size[2] - 100;
			}
		}
	}
}

void rhsInit(double* data, int* size, short* dipxy, short *dipxz)
{
	int size0 = size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));
	double* in = (double*)calloc(size0, sizeof(double));

	for (int i = 0; i < size0; i++) in[i] = (double)dipxy[i] / 1000.0 * (double)dipxy[i] / 1000.0;

	OpNabla::nablatX(in, size, data);
	OpNabla::nablatY(dipxy, 1000, size, tmp);
	for (long add = 0; add < size0; add++)
		data[add] += tmp[add];

	if ( size[2] > 1 && dipxz != nullptr )
	{
		for (int i = 0; i < size0; i++) in[i] = (double)dipxz[i] / 1000.0 * (double)dipxz[i] / 1000.0;
		OpNabla::nablatX(in, size, tmp);
		for (long add = 0; add < size0; add++) data[add] += tmp[add];
		OpNabla::nablatZ(dipxz, 1000, size, tmp);
		for (long add = 0; add < size0; add++)
			data[add] += tmp[add];
	}

	for (long add = 0; add < size0; add++) data[add] = -data[add];
	debug_data_save(data, size, "d:\\text.raw");
	free(tmp);
	free(in);
}

void opDirect(double* in, short* dipxy, short *dipxz, int* size, double epsilon, double* data)
{
	long size0 = (long)size[0] * size[1] * size[2];
	double* tmp = (double*)calloc(size0, sizeof(double));
	double* tmp2 = (double*)calloc(size0, sizeof(double));
	double* data2 = (double*)calloc(size0, sizeof(double));

	OpNabla::nablaX(in, size, data);
	for (long add = 0; add < size0; add++) data[add] *= (double)dipxy[add] / 1000.0;
	OpNabla::nablaY(in, size, tmp2);
	for (long add = 0; add < size0; add++) data[add] += tmp2[add];

	for (int i = 0; i < size0; i++) data2[i] = data[i] * (double)dipxy[i]/1000.0;

	OpNabla::nablatX(data2, size, tmp);
	OpNabla::nablatY(data, size, tmp2);
	for (long add = 0; add < size0; add++) tmp[add] += tmp2[add];

	for (long add = 0; add < size0; add++) data[add] = tmp[add];


	if (size[2] > 1 && dipxz != nullptr )
	{
		OpNabla::nablaX(in, size, tmp);
		for (long add = 0; add < size0; add++) tmp[add] *= (double)dipxz[add] / 1000.0;
		OpNabla::nablaZ(in, size, tmp2);
		for (long add = 0; add < size0; add++) tmp[add] += tmp2[add];

		for (int i = 0; i < size0; i++) data2[i] = tmp[i] * (double)dipxz[i] / 1000.0;

		OpNabla::nablatX(data2, size, tmp2);
		OpNabla::nablatZ(tmp, size, data2);
		for (long add = 0; add < size0; add++) tmp2[add] += data2[add];

		for (long add = 0; add < size0; add++) data[add] += tmp2[add];
	}
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
		InverseLaplacian2* lap2;
	};
public:
	myCallBack2();
	int* size, size0;
	short* dipxy, *dipxz;
	double epsilon;
	std::vector<XData> Xdipxy, Xdipxz;
	std::vector<CONSTRAINT> constraints;
	void setSize(int* size);
	void CallBack(void* in, void* out);
	void Preconditionner(void* in, void* out);
	void CallBack(std::vector<XData> Xin, std::vector<XData> Xout);
	void Preconditionner(std::vector<XData> Xin, std::vector<XData> Xout);
	void setDipxy(short* dipxy);
	void setDipxy(std::vector<XData> Xdipxy);
	void setDipxz(short* dipxy);
	void setDipxz(std::vector<XData> Xdipxy);
	void setEpsilon(double epsilon);
	void setConstraints(std::vector<CONSTRAINT> constraints);
private:
	PARAM* param;
	std::vector<PARAM*> Xparam;
	void preconditionnerConstraintsApply(std::vector<XData> data, std::vector<CONSTRAINT> constraints);
};

myCallBack2::myCallBack2()
{
	param = nullptr;
}

void myCallBack2::setSize(int* size)
{
	this->size = size;
}

void myCallBack2::setDipxy(short* dipxy)
{
	this->dipxy = dipxy;
}

void myCallBack2::setDipxy(std::vector<XData> data)
{
	this->Xdipxy = data;
}

void myCallBack2::setDipxz(short* dipxz)
{
	this->dipxz = dipxz;
}

void myCallBack2::setDipxz(std::vector<XData> data)
{
	this->Xdipxz = data;
}

void myCallBack2::setEpsilon(double epsilon)
{
	this->epsilon = epsilon;
}

void myCallBack2::setConstraints(std::vector<CONSTRAINT> constraints)
{
	this->constraints = constraints;
}


void myCallBack2::CallBack(void* in, void* out)
{
	opDirect((double*)in, dipxy, dipxz, size, epsilon, (double*)out);
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

		param->lap2 = new InverseLaplacian2();
		param->lap2->setdataIn(in);
		param->lap2->setDataOut(out);
		param->lap2->setSize(size);

	}
	// param->lap->run();
	param->lap2->run();

	// memcpy(out, in, (size_t)size[0] * size[1] * size[2] * sizeof(double));
}

void myCallBack2::CallBack(std::vector<XData> Xin, std::vector<XData> Xout)
{
	CHRONOS_TIC(CHRONOS, CHRONOS_OPDIRECT);
	int N = Xin.size();
	for (int n = 0; n < N; n++)
	{
		double* in = (double*)Xin[n].getData();
		int* size = Xin[n].getSize();
		double* out = (double*)Xout[n].getData();
		short* dxy = (short*)Xdipxy[n].getData();
		short* dxz = (short*)Xdipxz[n].getData();
		opDirect((double*)in, dxy, dxz, size, epsilon, (double*)out);
	}
	CHRONOS_TOC(CHRONOS, CHRONOS_OPDIRECT);
}

void myCallBack2::preconditionnerConstraintsApply(std::vector<XData> data, std::vector<CONSTRAINT> constraints)
{
	int N = constraints.size();
	for (int n = 0; n < N; n++)
	{
		int* no = this->constraints[n].no;
		int* x0 = this->constraints[n].x0;
		int* y0 = this->constraints[n].y0;
		int* z0 = this->constraints[n].z0;
		int nbre = this->constraints[n].nbre;
		double sum = 0.0;
		for (int i = 0; i < nbre; i++)
		{
			double* in = (double*)data[no[i]].getData();
			sum += in[this->size[0] * y0[i] + x0[i]];
		}
		sum /= (double)nbre;
		for (int i = 0; i < nbre; i++)
		{
			double* in = (double*)data[no[i]].getData();
			in[this->size[0] * y0[i] + x0[i]] = sum;
		}
	}
}

void myCallBack2::Preconditionner(std::vector<XData> Xin, std::vector<XData> Xout)
{
	CHRONOS_TIC(CHRONOS, CHRONOS_PRECOND);
	int N = Xin.size();
	if (Xparam.size() == 0)
	{
		Xparam.resize(N);
		for (int n = 0; n < N; n++)
		{
			void* in = Xin[n].getData();
			int* size = Xin[n].getSize();
			void* out = Xout[n].getData();
			PARAM* param0 = new myCallBack2::PARAM();
			param0->lap = new InverseLaplacian();
			param0->lap->setdataIn(in);
			param0->lap->setDataOut(out);
			param0->lap->setSize(size[0], size[1], size[2]);

			param0->lap2 = new InverseLaplacian2();
			param0->lap2->setdataIn(in);
			param0->lap2->setDataOut(out);
			param0->lap2->setSize(size);

			Xparam[n] = param0;
		}
	}

	preconditionnerConstraintsApply(Xin, this->constraints);
	for (int n = 0; n < N; n++)
	{
		Xparam[n]->lap2->run();
		// debug_data_save((double*)Xout[n].getData(), Xin[n].getSize(), "d:\\text.raw");
		// int* size = Xin[n].getSize(); memcpy(Xout[n].getData(), Xin[n].getData(), size[0] * size[1] * size[2] * sizeof(double));
	}
	preconditionnerConstraintsApply(Xout, this->constraints);
	CHRONOS_TOC(CHRONOS, CHRONOS_PRECOND);
}



#define YSWAP(x) ((x & 0xff) << 8) | ((x & 0xff00) >> 8)


int main_geotimeConstraint(int argc, char** argv)
{
	CHRONOS = CHRONOS_INIT;

	// int size[] = { 350, 2000, 1 };
	int size[] = { 200, 201, 202 };

	char* filename = "D:\\TOTAL\\PLI\\seismic.xt";
	XT_FILE* pf = new XT_FILE();
	pf->openForRead(filename);
	size[0] = pf->get_dimx();
	size[1] = pf->get_dimy();
	size[2] = 1;

	size[0] = 200; size[1] = 1001; size[2] = 202;

	long size0 = (long)size[0] * size[1] * size[2];
	short* seismic = new short[size0];
	short* dipxy = (short*)calloc(size0, sizeof(short));
	short* dipxz = (short*)calloc(size0, sizeof(short));;
	// pf->inlineRead(200, seismic);
	delete pf;
	
		
	int Nc = 3;
	CONSTRAINT constr;
	constr.no = (int*)calloc(Nc, sizeof(int));
	constr.x0 = (int*)calloc(Nc, sizeof(int));
	constr.y0 = (int*)calloc(Nc, sizeof(int));
	constr.z0 = (int*)calloc(Nc, sizeof(int));
	constr.nbre = Nc;
	constr.no[0] = 0; constr.x0[0] = 137; constr.y0[0] = 467;
	constr.no[1] = 0; constr.x0[1] = 137; constr.y0[1] = 531;
	constr.no[2] = 0; constr.x0[2] = 131; constr.y0[2] = 707;

	std::vector<CONSTRAINT> constraints;
	constraints.resize(1);
	constraints[0] = constr;


	// dipEstimation2D(seismic, size, 1.0, 1.5, dipxy);
	dipEstimationSynthetic3D(size, dipxy, dipxz);


	double* tau = (double*)calloc(size0, sizeof(double));	
	double* rhs = (double*)calloc(size0, sizeof(double));
	// rhsInit(rhs, size, dipxy, constraints);
	rhsInit(rhs, size, dipxy, dipxz);


	debug_data_save(rhs, size, "d:\\text.raw");
	debug_data_save(tau, size, "d:\\text.raw");


	std::vector<XData> Xdipxy;
	Xdipxy.resize(1);
	XData mdipxy;
	mdipxy.setData(dipxy, size, XData::TYPE::SHORT);
	Xdipxy[0] = mdipxy;

	std::vector<XData> Xdipxz;
	Xdipxz.resize(1);
	XData mdipxz;
	mdipxz.setData(dipxz, size, XData::TYPE::SHORT);
	Xdipxz[0] = mdipxz;

	std::vector<XData> Xrhs;
	Xrhs.resize(1);
	XData mrhs;
	mrhs.setData(rhs, size, XData::TYPE::DOUBLE);
	Xrhs[0] = mrhs;

	std::vector<XData> Xtau;
	Xtau.resize(1);
	XData mtau;
	mtau.setData(tau, size, XData::TYPE::DOUBLE);
	Xtau[0] = mtau;


	myCallBack2* c = new myCallBack2();
	c->setSize(size);
	c->setDipxy(dipxy);
	c->setDipxy(Xdipxy);
	c->setDipxz(Xdipxz);
	c->setEpsilon(0.005);
	// c->setConstraints(constraints);

	ConjugateGradient* p = new ConjugateGradient();
	p->setCallback(c, nullptr);
	p->setSize(size0);
	p->setRhs(rhs);
	p->setRhs(Xrhs);
	p->setX(tau);
	p->setX(Xtau);
	p->setNbiter(50);


	auto start = std::chrono::steady_clock::now();
	p->run2();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	double t_direct = CHRONOS_GET(CHRONOS, CHRONOS_OPDIRECT);
	double t_precond = CHRONOS_GET(CHRONOS, CHRONOS_PRECOND);

	fprintf(stderr, "tdirect:  %f\ntprecond: %f\n", t_direct, t_precond);


	delete p;

	for (long z=0; z<size[2]; z++)
	{
		for (long y = 0; y < size[1]; y++)
		{
			for (long x = 0; x < size[0]; x++)
			{
				tau[size[0]*size[1]*z + size[0] * y + x] += x;
			}
		}
	}
	debug_data_save(tau, size, "d:\\text.raw");

	fprintf(stderr, "%s ok\n", __FILE__);
	return 0;
}


