
#include <stdio.h>
#include <ConjugateGradient.h>

class myCallBack : public ConjugateGradient::CALLBACK
{
public:
	int size;
	void setSize(int size);
	void CallBack(void* in, void* out);
	void Preconditionner(void* in, void* out);
};

void myCallBack::setSize(int size)
{
	this->size = size;
}

void myCallBack::CallBack(void* in, void* out)
{
	double* pin = (double*)in, * pout = (double*)out;
	pout[0] = 283.0 * pin[0] + 368.0 * pin[1] + 318.0 * pin[2];
	pout[1] = 368.0 * pin[0] + 480.0 * pin[1] + 412.0 * pin[2];
	pout[2] = 318.0 * pin[0] + 412.0 * pin[1] + 362.0 * pin[2];
}


void myCallBack::Preconditionner(void* in, void* out)
{
	for (long add = 0; add < size; add++)
	{
		((double*)out)[add] = ((double*)in)[add];
	}
}


int main_conjugateGradient(int argc, char** argv)
{
	int size[] = { 3, 3, 1 };
	double rhs[] = { 4880.0, 6344.0, 5504.0 };
	double x[] = { 0.0, 0.0, 0.0 };
	long size0 = (long)size[0] * size[1] * size[2];
	myCallBack* c = new myCallBack();
	c->setSize(3);
	ConjugateGradient* p = new ConjugateGradient();

	p->setCallback(c, nullptr);
	p->setSize(3);
	p->setRhs(rhs);
	p->setX(x);
	p->setNbiter(4);
	p->run();

	delete p;
	fprintf(stderr, "%s ok\n", __FILE__);
	return 0;
}