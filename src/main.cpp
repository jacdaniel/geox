
#include <stdio.h>


void main_constraint_geotime_run(int argc, char **argv);
int main_conjugateGradient(int argc, char** argv);
int main_geotimeConstraint(int argc, char** argv);
int main_nabla(int argc, char** argv);


int main(int argc, char** argv)
{

	// main_constraint_geotime_run(argc, argv);
	// main_conjugateGradient(argc, argv);
	main_geotimeConstraint(argc, argv);
	// main_nabla(argc, argv);

	fprintf(stderr, "finish ok\n");
	return 1;
}