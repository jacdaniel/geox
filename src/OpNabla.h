// #pragma once


#ifndef __OPNABLA__
#define __OPNABLA__


class OpNabla
{
private:
	template <typename T> static int nablaX(T* f, int* size, T* g);
	template <typename T> static int nablaY(T* f, int* size, T* g);
	template <typename T> static int nablaZ(T* f, int* size, T* g);
	template <typename T> static int nablatX(T* f, int* size, T* g);
	template <typename T> static int nablatY(T* f, int* size, T* g);
	template <typename T> static int nablatZ(T* f, int* size, T* g);

	template <typename T> static int nablaX(short* f, double norm, int* size, T* g);
	template <typename T> static int nablaY(short* f, double norm, int* size, T* g);
	template <typename T> static int nablaZ(short* f, double norm, int* size, T* g);
	template <typename T> static int nablatX(short* f, double norm, int* size, T* g);
	template <typename T> static int nablatY(short* f, double norm, int* size, T* g);
	template <typename T> static int nablatZ(short* f, double norm, int* size, T* g);

public:
	static int nablaX(float* f, int* size, float* g);
	static int nablaY(float* f, int* size, float* g);
	static int nablaZ(float* f, int* size, float* g);
	static int nablaX(double* f, int* size, double* g);
	static int nablaY(double* f, int* size, double* g);
	static int nablaZ(double* f, int* size, double* g);

	static int nablatX(float* f, int* size, float* g);
	static int nablatY(float* f, int* size, float* g);
	static int nablatZ(float* f, int* size, float* g);
	static int nablatX(double* f, int* size, double* g);
	static int nablatY(double* f, int* size, double* g);
	static int nablatZ(double* f, int* size, double* g);

	static int nablaX(short* f, double norm, int* size, double* g);
	static int nablaY(short* f, double norm, int* size, double* g);
	static int nablaZ(short* f, double norm, int* size, double* g);
	static int nablatX(short* f, double norm, int* size, double* g);
	static int nablatY(short* f, double norm, int* size, double* g);
	static int nablatZ(short* f, double norm, int* size, double* g);
};


#endif