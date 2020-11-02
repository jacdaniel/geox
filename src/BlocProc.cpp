
#include <algorithm>
#include <BlocProc.h>

#include <stdio.h>

BlocProc::BlocProc()
{
	for (int i = 0; i < 3; i++) this->bloc_size[i] = 32;
	for (int i = 0; i < 3; i++) this->size[i] = 1;
	for (int i = 0; i < 3; i++) this->border_size[i] = 0;
	for (int i = 0; i < 3; i++) this->Nblocs[i] = 0;
	this->blocks = nullptr;
	this->compute_type = COMPUTE_TYPE::SAME;
	pcallback = nullptr;
}

void BlocProc::set_bloc_size(int* size)
{
	if (size == nullptr) return;
	for (int i = 0; i < 3; i++) this->bloc_size[i] = size[i];
}

void BlocProc::set_size(int* size)
{
	if (size == nullptr) return;
	for (int i = 0; i < 3; i++) this->size[i] = size[i];
}

void BlocProc::set_border_size(int* size)
{
	if (size == nullptr) return;
	for (int i = 0; i < 3; i++) this->border_size[i] = size[i];
}

void BlocProc::set_compute_type(int type)
{
	this->compute_type = type;
}

void BlocProc::setCallBack(CALLBACK* callback)
{
	this->pcallback = callback;
}

void BlocProc::init_param()
{
	for (int i = 0; i < 3; i++)
	{
		this->Nblocs[i] = (this->size[i] - 1) / this->bloc_size[i] + 1;

	}
}

//void block_dim_cut_param2(long size, long x0, long Nu, long border, int* x1, int* x2, int* xe1, int* xe2)
//{
//	*x1 = std::min(x0, size - 1);
//	*x2 = std::min(x0 + Nu - 1, size - 1);
//	*xe1 = std::max((const int)(*x1 - border), 0);
//	*xe2 = std::min(*x2 + border, size - 1);
//}

void BlocProc::block_dim_cut_param2(long size, long x0, long Nu, long border, PARAM* p)
{
	p->x1 = std::min(x0, size - 1);
	p->x2 = std::min(x0 + Nu - 1, size - 1);
	p->xread1 = std::max((const int)(p->x1 - border), 0);
	p->xread2 = std::min(p->x2 + border, size - 1);
	p->dx = p->xread2 - p->xread1 + 1;
}

void BlocProc::run()
{
	// if (this->pcallback == nullptr) return;
	init_param();
	PARAM px, py, pz;
	int x1[3], x2[3], xe1[3], xe2[3];

	for (int z = 0; z < this->Nblocs[2]; z++)
	{
		int z0 = z * this->bloc_size[2];
		block_dim_cut_param2(this->size[2], z0, this->bloc_size[2], this->border_size[2], &pz);

		for (int y = 0; y < this->Nblocs[1]; y++)
		{
			int y0 = y * this->bloc_size[1];
			block_dim_cut_param2(this->size[1], y0, this->bloc_size[1], this->border_size[1], &py);

			for (int x = 0; x < this->Nblocs[0]; x++)
			{
				int x0 = x * this->bloc_size[0];
				block_dim_cut_param2(this->size[0], x0, this->bloc_size[0], this->border_size[0], &px);

				// fprintf(stderr, "xread1 %d %d %d\n", px.xread1, py.xread1, pz.xread1);
				// fprintf(stderr, "xread2 %d %d %d\n", px.xread2, py.xread2, pz.xread2);
				// fprintf(stderr, "dx     %d %d %d\n", px.dx, py.dx, pz.dx);


				x1[0] = std::min(x * this->bloc_size[0], size[0] - 1);
				x1[1] = std::min(y * this->bloc_size[1], size[1] - 1);
				x1[2] = std::min(z * this->bloc_size[2], size[2] - 1);
				for (int i = 0; i < 3; i++)
				{
					x2[i] = std::min(x1[i] + bloc_size[i] - 1, size[i] - 1);
					xe1[i] = std::max(x1[i] - border_size[i], 0);
					xe2[i] = std::min(x2[i] + border_size[i], size[i] - 1);
				}

				pcallback->f(px.xread1, px.xread2, py.xread1, py.xread2, pz.xread1, pz.xread2, x1, x2, xe1, xe2);



			}

		}
	}
}