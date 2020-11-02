// #pragma once

#ifndef __CONJUGATEGRADIENT__
#define __CONJUGATEGRADIENT__

#include<vector>
#include <XData.h>


class ConjugateGradient
{
public:
	class CALLBACK
	{
		public:
			virtual void CallBack(void* in, void* out) {};
			virtual void Preconditionner(void* in, void* out) {};

			virtual void CallBack(std::vector<XData> xin, std::vector<XData> xout) {}
			virtual void Preconditionner(std::vector<XData> xin, std::vector<XData> xout) {}
	};

	class PARAM
	{
	public:
		double* tmp, * z, * r, * d;
		std::vector<XData> Xtmp, Xz, Xr, Xd;
		PARAM(long size);
		PARAM(std::vector<XData> data);
		~PARAM();
	};


private:
	long nbIter, size0;
	double* rhs, * x;

	std::vector<XData> Xrhs;
	std::vector<XData> Xx;

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
	int run2();

	double dot(std::vector<XData> v1, std::vector<XData> v2);
	void setRhs(std::vector<XData> rhs);
	void setX(std::vector<XData> x);

	PARAM* param;

};


#endif
