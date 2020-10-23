
#include <malloc.h>
#include <math.h>
#include <BlocProc.h>
#include <Gaussian.h>
#include <Convolution.h>
#include <DataToNormal.h>
#include "DipProcessing.h"


void debug_data_float_save(float* data, int* size, char* filename);
void debug_data_short_save(short* data, int* size, char* filename);



void DipProcessing::setDataIn(void* in) 
{
	this->dataIn = in;
}

void DipProcessing::setSize(int* size)
{
	this->sizeIn = size;
}

void DipProcessing::setBlocSize(int* size)
{
	this->blockSize = size;
}

void DipProcessing::setSigmaGradient(double sigma)
{
	this->sigmaGrad = sigma;
}

void DipProcessing::setSigmaTensor(double sigma)
{
	this->sigmaTens = sigma;
}

void DipProcessing::setDipxXY(void* data)
{
	this->dipXY = data;
}

void DipProcessing::setDipXZ(void* data)
{
	this->dipXZ = data;
}




class mycallback2 : public BlocProc::CALLBACK
{
private:
	class PARAM
	{
	public:
		int sizeMaskD, sizeMask2;
		float* smooth, *deriv, *smooth2, * gx, * gy, * gz;
		short* res;
		DataToNormal* pDataToNormal;
	};

	int* x1, * x2, * xe1, * xe2;
	int* border, *size0, *sizeC, *blockSize;
	short* data0 = nullptr, * data_e;
	float *data_c;
	short* dipXY;
	PARAM* param = nullptr;

public:
	mycallback2();
	void setData0(short* data);
	void setDataE(short* data);
	void setDataC(float* data);
	void setSize0(int* size0);
	void setSizeC(int* size);
	void setBorder(int *size);
	void setDipXY(short* data);
	void setBlockSize(int* size);

	void f(int xread1, int xread2, int yread1, int yread2, int zread1, int zread2, int* x1, int* x2, int* xe1, int* xe2);
};

mycallback2::mycallback2()
{
	this->param = nullptr;
}

void mycallback2::setData0(short* data)
{
	this->data0 = data;
}

void mycallback2::setDataE(short* data)
{
	this->data_e = data;
}

void mycallback2::setDataC(float* data)
{
	this->data_c = data;
}

void mycallback2::setSize0(int* size0)
{
	this->size0 = size0;
}

void mycallback2::setSizeC(int* size)
{
	this->sizeC = size;
}

void mycallback2::setBorder(int* size)
{
	this->border = size;
}

void mycallback2::setDipXY(short* data)
{
	this->dipXY = data;
}

void mycallback2::setBlockSize(int* size)
{
	this->blockSize = size;
}


void blocExtract(short* data0, int* size0, int xread1, int xread2, int yread1, int yread2, int zread1, int zread2, short* data_e, int dx, int dy, int dz)
{
	for (int z = 0; z < dz; z++)
	{
		for (int y = 0; y < dy; y++)
		{
			for (int x = 0; x < dx; x++)
			{
				// if (z + xe1[2] < size0[2] && y < size0[0] && x + xe1[0] < size0[0])
				{
					long add_e = dx * dy * z + dx * y + x;
					long add0 = size0[0] * size0[1] * (z + zread1) + size0[0] * (y + yread1) + (x + xread1);
					data_e[add_e] = data0[add0];
				}
			}
		}
	}

}

#define MODULOX2(x, m, res) {\
 res = x; \
 if ( res < 0 ) res = -res; \
 res = (res) % (2*(m)); \
 if ( res >= m ) res = 2 * (m) - res - 1; }


void data_mirror_short_float(short* in, int* size_in, int* offset, float* out, int* size_out)
{
	int xx, yy, zz;

	for (int z = 0; z < size_out[2]; z++)
	{
		MODULOX2(z - offset[2], size_in[2], zz);
		for (int y = 0; y < size_out[1]; y++)
		{
			MODULOX2(y - offset[1], size_in[1], yy);
			for (int x = 0; x < size_out[0]; x++)
			{
				MODULOX2(x - offset[0], size_in[0], xx);
				long add = size_out[0] * size_out[1] * z + size_out[0] * y + x;
				long add0 = size_in[0] * size_in[1] * zz + size_in[0] * yy + xx;
				out[add] = in[add0];
			}
		}
	}
}


// void mycallback2:f(void* _param, int* x1, int* x2, int* xe1, int* xe2)
void mycallback2::f(int xread1, int xread2, int yread1, int yread2, int zread1, int zread2, int* x1, int* x2, int* xe1, int* xe2)
{
	int size_in[3], offset[3];
	int sizeB = sizeC[0] * sizeC[1] * sizeC[2];

	if (this->param == nullptr)
	{
		param = new mycallback2::PARAM();
		param->res = new short[sizeB];	
		param->pDataToNormal = new DataToNormal();
		param->pDataToNormal->setSigmaGradient(1.0);
		param->pDataToNormal->setSigmaTensor(1.5);
		param->pDataToNormal->setSize(sizeC);
	}
	
	for (int i = 0; i < 3; i++)
		size_in[i] = xe2[i] - xe1[i] + 1;

	for (int i = 0; i < 3; i++)
		offset[i] = border[i] - (x1[i] - xe1[i]);
	
	int dx = xread2 - xread1 + 1;
	int dy = yread2 - yread1 + 1;
	int dz = zread2 - zread1 + 1;

	int sizex[] = { dx, dy, dz };

	blocExtract(data0, size0, xread1, xread2, yread1, yread2, zread1, zread2, data_e, dx, dy, dz);
	debug_data_short_save(data_e, sizex, "d:\\data_e.raw");
	data_mirror_short_float(data_e, sizex, offset, data_c, sizeC);
	debug_data_float_save(data_c, sizeC, "d:\\data_c.raw");

	param->pDataToNormal->dataToDip(data_c, sizeC, param->res);


	// dataToDipx

	for (int y=0; y< this->blockSize[1]; y++)
		for (int x = 0; x < this->blockSize[0]; x++)
		{
			int add0 = sizeC[0] * (y+ border[1]) + x + border[0];
			int add = size0[0] * (y + x1[1]) + x + x1[0];
			if (x + x1[0] < size0[0] && y + x1[1] < size0[1])
				this->dipXY[add] = param->res[add0];
		}

	int i = 0;
	debug_data_short_save(dipXY, size0, "d:\\dipXY.raw");



	/*
	for (int z = 0; z < size_in[2]; z++)
	{
		for (int y = 0; y < size_in[1]; y++)
		{
			for (int x = 0; x < size_in[0]; x++)
			{
				if (z + xe1[2] < size0[2] && y < size0[0] && x + xe1[0] < size0[0])
				{
					// param->data0[(long)size_in[0] * size_in[1] * z + (long)size_in[0] * y + x] = param->cache[(long)param->cache_size[0] * param->cache_size[1] * z + (long)param->cache_size[0] * y + x + xe1[0]];
				}
			}
		}
	}
	*/



	/*


	// to improve
	if (param->cache)
	{
		if (xe1[0] == 0)
		{
			fprintf(stderr, "cache\n");
			fileio_bloc_read_short(param->file, 0, xe1[1], xe1[2], param->native_size[0], param->size[1], param->size[2], param->cache);
		}
		for (int z = 0; z < size_in[2]; z++)
		{
			for (int y = 0; y < size_in[1]; y++)
			{
				for (int x = 0; x < size_in[0]; x++)
				{
					if (z < param->cache_size[2] && y < param->cache_size[1] && x + xe1[0] < param->cache_size[0])
					{
						param->data0[(long)size_in[0] * size_in[1] * z + (long)size_in[0] * y + x] = param->cache[(long)param->cache_size[0] * param->cache_size[1] * z + (long)param->cache_size[0] * y + x + xe1[0]];
					}
				}
			}
		}
	}
	else
	{
		fileio_bloc_read_short(param->file, xe1[0], xe1[1], xe1[2], xe2[0] - xe1[0] + 1, xe2[1] - xe1[1] + 1, xe2[2] - xe1[2] + 1, param->data0);
	}
	data_mirror_short_float(param->data0, size_in, offset, param->data, param->size);
	acp_gradient_run(param->p, param->data, param->size);
	int bloc_size[3];
	int pos[3];
	for (int n = 0; n < param->nbscales; n++)
	{
		int N = pow2(n);
		for (int i = 0; i < 3; i++) { bloc_size[i] = param->bloc_size[i] / N; }
		for (int i = 0; i < 3; i++) { pos[i] = x1[i] / N; }


		if (param->ni[n] != NULL && param->vector[n] != NULL)
		{
			for (int i = 0; i < 3; i++)
				data_bloc_insert_float(param->ni[n][i], param->out_size[n], param->vector[n][i], bloc_size, pos, NULL);
		}
	}
	*/
}




void DipProcessing::run()
{
	int b = (int)ceil(3.0 * this->sigmaGrad + 3.0 * this->sigmaTens);
	int borderSize[] = { b, b, 1 };
	int size_e[3];
	for (int i = 0; i < 3; i++) if (blockSize[i] > 1)  size_e[i] = this->blockSize[i] + 2 * borderSize[i]; else size_e[i] = 1;
	long s0 = size_e[0] * size_e[1] * size_e[2];
	short* data_e = (short*)calloc(s0, sizeof(short));
	float* data_c = (float*)calloc(s0, sizeof(float));
	
	mycallback2* c = new mycallback2();
	c->setData0((short*)this->dataIn);
	c->setDataE((short*)data_e);
	c->setDataC((float*)data_c);
	c->setSize0(this->sizeIn);
	c->setSizeC(size_e);
	c->setBorder(borderSize);
	c->setDipXY((short*)this->dipXY);
	c->setBlockSize(blockSize);

	BlocProc * p = new BlocProc();
	p->set_bloc_size(this->blockSize);
	p->set_size(this->sizeIn);
	p->set_border_size(borderSize);
	p->set_compute_type(BlocProc::UTIL);
	p->setCallBack(c);
	p->run();

	delete c;
	delete p;
}