#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include <memory>

namespace N_body_simulation {

constexpr double gravitational_constant=6.674e-11; 
constexpr double solar_mass=1.98892e30;

struct Body {	
	double x;
	double y;
	
	double v_x;
	double v_y;
	
	double a_x;
	double a_y;
	
	double m; 
};

class BS {
	//~ BS (Body Solver) is an abstract class. The method "Advance" has 
	//~ to be implimented by subclasses (i.e. Euler, RK4 etc.).
public:
	using id=std::vector<Body>::size_type;
	
	explicit BS (std::vector<Body> BB)
		: B {BB} {}
	BS (std::vector<Body> BB, double tss, double sc)
		: B {BB}, dt {tss}, softening_constant {sc} {}
	
	void solve_for(unsigned int num_steps);

protected:
	virtual void Advance() = 0;
	void ComputeAcceleration(std::vector<Body> *R, std::vector<Body> *W);
	void WritePositionToFile(std::fstream &x_fst, std::fstream &y_fst);

	std::vector<Body> B;
	double dt = 1e11; // Time step size
	double softening_constant = 3e4;
};

class Euler_std : public BS {
public:
	explicit Euler_std(std::vector<Body> BB)
		: BS (BB) {}
	Euler_std (std::vector<Body> BB, double tss, double sc)
		: BS (BB, tss, sc) {}

	void Advance();
};

}
