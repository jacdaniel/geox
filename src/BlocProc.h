// #pragma once
#ifndef __BLOCPROC__
#define __BLOCPROC__

class BlocProc
{
public:
	enum COMPUTE_TYPE { SAME, UTIL };
	class CALLBACK
	{
	public:
		virtual void f(int xread1, int xread2, int yread1, int yread2, int zread1, int zread2, int *x1, int *x2, int *xe1, int *xe2) {}
	};

	BlocProc();
	void set_bloc_size(int* size);
	void set_size(int* size);
	void set_border_size(int* size);
	void set_compute_type(int type);
	void setCallBack(CALLBACK* callback);
	void run();

private:
	typedef struct _PARAM
	{
		int x1, x2, xread1, xread2, xoffset, dx;
	}PARAM;

	void init_param();
	void block_dim_cut_param2(long size, long x0, long Nu, long border, PARAM* p);

	CALLBACK* pcallback;
	int bloc_size[3];
	int size[3];
	int border_size[3];
	int compute_type;

	int Nblocs[3], ** blocks;
};



#endif
