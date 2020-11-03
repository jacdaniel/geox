
#include <xchronos.h>


Xchronos::Xchronos()
{
	nbre = 10;

	start = new std::chrono::steady_clock::time_point[nbre];
	end = new std::chrono::steady_clock::time_point[nbre];
	elapsed = new double[nbre];
	for (int n = 0; n < nbre; n++)
		reset(n);
}


Xchronos::~Xchronos()
{

}


void Xchronos::reset(int no)
{
	elapsed[no] = 0.0;
}

void Xchronos::tic(int no)
{
	start[no] = std::chrono::steady_clock::now();
}


void Xchronos::toc(int no)
{
	end[no] = std::chrono::steady_clock::now();
}

double Xchronos::get(int no)
{
	//std::chrono::duration<double> elapsed_seconds += hours(end[no] - start[no]);

	//return (double)elapsed_seconds.count();


	//std::chrono::duration<double> elapsed_seconds;
	//
	//std::timepoint_t t0;

	return 0.0;
}



