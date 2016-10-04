#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>

namespace N_body_simulation {

const double gravitational_constant=6.674e-11; 
const double solar_mass=1.98892e30;

struct Body {
	double x_pos;
	double y_pos;
	
	double x_vel;
	double y_vel;
	
	double x_acc;
	double y_acc;
	
	double mass; 
};

class NBS {
	using id=std::vector<Body>::size_type;
public:
	explicit NBS (std::vector<Body> BB)
		: B {BB} {}
	NBS (std::vector<Body> BB, double tss, double sc)
		: B {BB}, time_step_size {tss},softening_constant {sc} {}
	
	void solve(unsigned int num_steps);
	void solve_next();

private:
	void solve_acc(id); 
	void solve_vel(id); 
	void solve_pos(id);

	std::vector<Body> B;
	double time_step_size = 1e11;
	double softening_constant = 3e4;
};

std::vector<Body> BuildBodies(std::vector<Body>::size_type N, double r);

}
