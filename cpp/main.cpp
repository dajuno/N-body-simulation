#include "NBS.h"
using namespace N_body_simulation;

int main() {
	try
	{
		
	int N = 100;
	double radius = 1e18;
	double time_step_size = 1e4;
	
	std::vector<Body> B = BuildBodies(N,radius);
	NBS nbs {B};
	nbs.solve(time_step_size);
		
	return 0;
	}
	catch (std::exception& e)
	{
		std::cerr<<e.what()<<"\n";
		return 1;
	}
	catch (...)
	{
		std::cerr<<"Some exception\n";
		return 2;
	}
}
