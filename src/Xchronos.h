// #pragma once


#include <chrono>

class Xchronos
{
private:
	int nbre;
	std::chrono::steady_clock::time_point* start, * end;
	double* elapsed;
public:
	Xchronos();
	~Xchronos();
	void reset(int no);
	void tic(int no);
	void toc(int no);
	double get(int no);
};
