#include "NBS.h"
using namespace N_body_simulation;

int main() {
	try
	{
		
	int N = 1e2;
	double num_steps = 1e3;
	
	std::vector<Body> B = BuildBodies(N);
	NBS nbs {B};
	nbs.solve(num_steps);
		
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
