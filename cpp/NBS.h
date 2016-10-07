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
	double x=0;
	double y=0;
	
	double v_x=0;
	double v_y=0;
	
	double a_x=0;
	double a_y=0;
	
	double m=0; 
};

bool body_in_quadrant(Body b, Body q, double l);

struct Node {
public:
	Node(Body qq, double ll) : q {qq}, l {ll} {}
	Node();
	
	Body q; // Quadrant
	double l; // Length of quadrant
	
	Body b;	// Body
	bool empty = true;

	std::shared_ptr<Node> NE = nullptr;
	std::shared_ptr<Node> NW = nullptr;
	std::shared_ptr<Node> SE = nullptr;
	std::shared_ptr<Node> SW = nullptr;
};

void insert_in_tree(std::shared_ptr<Node>, Body);
std::shared_ptr<Node> BuildTree(std::vector<Body> B, double length);
void print_tree(std::shared_ptr<Node> node);

void GetBodiesFromTree(std::shared_ptr<Node> node, Body b, std::vector<Body>* TB);

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
	void ComputeAcceleration(const std::vector<Body> &R, std::vector<Body> *W);
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
