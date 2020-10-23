// #pragma once

#ifndef __CONJUGATEGRADIEN__
#define __CONJUGATEGRADIENT__


class ConjugateGradient
{
public:
	class CALLBACK
	{
		public:
			virtual void CallBack(void* in, void* out) {};
			virtual void Preconditionner(void* in, void* out) {};
	};

	class PARAM
	{
	public:
		double* tmp, * z, * r, * d;
		PARAM(long size);
		~PARAM();
	};


private:
	long nbIter, size0;
	double* rhs, * x;
	CALLBACK* callback;
	double dot(double* v1, long size, double* v2);
	double dot(float* v1, long size, float* v2);
public:
	ConjugateGradient();
	~ConjugateGradient();

	void setCallback(void* callback, void* data);
	void setPreconditionner(void* f, void* data);
	void setSize(int size);
	void setRhs(void* rhs);
	void setX(void* x);
	void setNbiter(int nbiter);
	int run();
	PARAM* param;

};


#endif
