

#include <stdio.h>
#include <math.h>
#include "ConjugateGradient.h"

ConjugateGradient::PARAM::PARAM(long size)
{
	this->z = new double[size];
	this->tmp = new double[size];
	this->r = new double[size];
	this->d = new double[size];
}

ConjugateGradient::PARAM::~PARAM()
{
	delete this->z;
	delete this->tmp;
	delete this->r;
	delete this->d;
}



ConjugateGradient::ConjugateGradient()
{
	nbIter = 100;
	size0 = 0;
	rhs = nullptr;
	x = nullptr;
	param = nullptr;
	this->callback = nullptr;
}

ConjugateGradient::~ConjugateGradient()
{
}

void ConjugateGradient::setCallback(void* callback, void* data)
{
	this->callback = (CALLBACK*)callback;
}

void ConjugateGradient::setPreconditionner(void* f, void* data)
{
}

void ConjugateGradient::setSize(int size)
{
	this->size0 = size;
}

void ConjugateGradient::setRhs(void* rhs)
{
	this->rhs = (double*)rhs;
}

void ConjugateGradient::setX(void* x)
{
	this->x = (double*)x;
}

void ConjugateGradient::setNbiter(int nbiter)
{
	this->nbIter = nbiter;
}


double ConjugateGradient::dot(double* v1, long size, double* v2)
{
	double d = 0.0;
	for (long add = 0; add < size; add++)
	{
		d += v1[add] * v2[add];
	}
	return d;
}

double ConjugateGradient::dot(float* v1, long size, float* v2)
{
	double d = 0.0;
	for (long add = 0; add < size; add++)
	{
		d += (double)v1[add] * v2[add];
	}
	return d;
}

int ConjugateGradient::run()
{
	double db = 0.0, rho0 = 0.0, rho1 = 0.0, beta = 0.0, denom = 1.0, alphax = 1.0;
	int cont = 1;
	double err = 0.0;

	const double minus_one = -1.0;
	const double plus_one = 1.0;

	if (param == nullptr)
	{
		param = new ConjugateGradient::PARAM(this->size0);
	}

	db = dot(rhs, size0, rhs);
	for (long i = 0; i < size0; i++) param->r[i] = rhs[i];
	callback->CallBack(x, param->tmp);
	for (long i = 0; i < size0; i++) param->r[i] -= param->tmp[i];
	callback->Preconditionner(param->r, param->z);
	rho0 = dot(param->r, size0, param->z);
	int iter = 0;

	while (cont)
	{
		if (iter % 10) fprintf(stderr, "%d \ %d [%G - %G]\n", iter, nbIter, alphax, denom);
		if (iter == 0)
		{
			for (long i = 0; i < size0; i++)
			{
				param->d[i] = param->z[i];
			}
		}
		else
		{
			for (long i = 0; i < size0; i++)
			{
				beta = rho0 / rho1;
				param->d[i] = param->z[i] + beta * param->d[i];
			}
		}
		callback->CallBack(param->d, param->tmp);
		denom = dot(param->d, size0, param->tmp);
		alphax = rho0 / denom;
		if ( fabs(denom) > 1e-6 && fabs(alphax) < 1e6 && fabs(alphax) > 1e-6)
		{			
			for (long i = 0; i < size0; i++)
			{
				x[i] += alphax * param->d[i];
			}
			for (long i = 0; i < size0; i++)
			{
				param->r[i] -= alphax * param->tmp[i];
			}
			callback->Preconditionner(param->r, param->z);
			rho1 = rho0;
			rho0 = dot(param->r, size0, param->z);
		}
		else
		{
			cont = 0;
		}
		
		iter++;
		if (iter >= nbIter)
		{
			cont = 0;
		}
	}

	return 0;
}

