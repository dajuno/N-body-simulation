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

class BodySolver {
	//~ BodySolver is an abstract class. The method "Advance" has 
	//~ to be implimented by subclasses (i.e. Euler, RK4 etc.).
public:
	using id=std::vector<Body>::size_type;
	
	explicit BodySolver (std::vector<Body> BB)
		: B {BB} {}
	BodySolver (std::vector<Body> BB, double tss, double sc)
		: B {BB}, dt {tss}, softening_constant {sc} {}
	
	void SolveTimeEvolution(unsigned int num_steps);

protected:
	virtual void Advance() = 0;
	void ComputeAcceleration(std::vector<Body> *B);
    double ComputeTotalEnergy();
	void WritePositionToFile(std::fstream &x_fst, std::fstream &y_fst);
    void WriteEnergyToFile(std::fstream &e_fst, double E);

	std::vector<Body> B;
	double dt = 1e11; // Time step size
	double softening_constant = 3e4;
};


class Euler : public BodySolver {
public:
	explicit Euler(std::vector<Body> BB)
		: BodySolver (BB) {}
	Euler(std::vector<Body> BB, double tss, double sc)
		: BodySolver (BB, tss, sc) {}

	void Advance();
};

class EulerAcc : public BodySolver {
public:
	explicit EulerAcc(std::vector<Body> BB)
		: BodySolver (BB) {}
	EulerAcc(std::vector<Body> BB, double tss, double sc)
		: BodySolver (BB, tss, sc) {}

	void Advance();
};

class EulerStabilized : public BodySolver {
public:
	explicit EulerStabilized(std::vector<Body> BB)
		: BodySolver (BB) {}
	EulerStabilized(std::vector<Body> BB, double tss, double sc)
		: BodySolver (BB, tss, sc) {}

	void Advance();
};

class EulerAccImproved : public BodySolver {
public:
	explicit EulerAccImproved(std::vector<Body> BB)
		: BodySolver (BB) {}
	EulerAccImproved(std::vector<Body> BB, double tss, double sc)
		: BodySolver (BB, tss, sc) {}

	void Advance();

private:
    std::vector<Body> B_prec;
};

class EulerImproved : public BodySolver {
public:
	explicit EulerImproved(std::vector<Body> BB)
		: BodySolver (BB) {}
	EulerImproved(std::vector<Body> BB, double tss, double sc)
		: BodySolver (BB, tss, sc) {}

	void Advance();

private:
    std::vector<Body> B_prec;
};

}
