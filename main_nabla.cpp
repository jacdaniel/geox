
#include <malloc.h>
#include <stdio.h>

#include <OpNabla.h>



int main_nabla(int argc, char** argv)
{
	int size[] = { 10, 11, 1 };
	int size0 = size[0] * size[1] * size[2];

	double* in = (double*)calloc(size0, sizeof(double));
	double* out = (double*)calloc(size0, sizeof(double));
	double* out2 = (double*)calloc(size0, sizeof(double));

	in[11*5+5] = 1;

	OpNabla::nablaY(in, size, out);
	OpNabla::nablatY(out, size, out2);


	return 0;
}