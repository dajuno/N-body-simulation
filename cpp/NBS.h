#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>

namespace N_body_simulation {

const double gravitational_constant=6.674e-11; 
const double solar_mass=1.98892e30;

struct Body {
	double x;
	double y;
	
	double xv;
	double yv;
	
	double xa;
	double ya;
	
	double m; 
};

struct Force {
    double x;
    double y;
};

class NBS {
	using id=std::vector<Body>::size_type;
public:
	explicit NBS (std::vector<Body> BB)
		: B {BB} {}
	NBS (std::vector<Body> BB, double tss, double sc)
		: B {BB}, dt {tss},softening_constant {sc} {}
	
	void solve(unsigned int num_steps);

private:
	void solve_next();
	void solve_acc(id);
	void euler_improved();
	void forces(std::vector<Force>& F);

	std::vector<Body> B;
	double dt = 1e11; // Time step size
	double softening_constant = 3e4;
};

// N = number of bodies
std::vector<Body> BuildBodies(std::vector<Body>::size_type N);

}
