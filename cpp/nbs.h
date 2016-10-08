#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include <memory>

namespace N_body_simulation {


//~ tpl
constexpr double gravitational_constant=6.674e-11; 
constexpr double solar_mass=1.98892e30;
struct tpl {
	double x=0;
	double y=0;
	
tpl& operator+=(const tpl& rhs) {
	x += rhs.x;
	y += rhs.y;
	return *this;
}
tpl& operator-=(const tpl& rhs) {
	x -= rhs.x;
	y -= rhs.y;
	return *this;
}
tpl& operator*=(const double& rhs) {
	x *= rhs;
	y *= rhs;
	return *this;
}
tpl& operator/=(const double& rhs) {
	x /= rhs;
	y /= rhs;
	return *this;
}
bool operator==(const tpl& rhs) {
	return x == rhs.x && y == rhs.y;
}
bool operator<=(const tpl& rhs) {
	return x <= rhs.x && y <= rhs.y;
}
bool operator<(const tpl& rhs) {
	return x < rhs.x && y < rhs.y;
}
bool operator>(const tpl& rhs) {
	return x > rhs.x && y > rhs.y;
}
bool operator>=(const tpl& rhs) {
	return x >= rhs.x && y >= rhs.y;
}
};
std::ostream& operator<<(std::ostream& ost, const tpl& tuple);
tpl operator+(const tpl& lhs, const tpl& rhs);
tpl operator-(const tpl& lhs, const tpl& rhs);
tpl operator*(const double& lhs, const tpl& rhs);
double operator*(const tpl lhs, const tpl& rhs);

//~ Body
struct Body {	
	tpl r;
	tpl v;
	tpl a;
	double m;
};
std::ostream& operator<<(std::ostream& ost, const Body& b);

//~ QuadTree
bool BodyIsInQuadrant(Body b, Body q, double l);
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
std::shared_ptr<Node> BuildTree(std::vector<Body> *B, double length);
void print_tree(std::shared_ptr<Node> node);
void GetBodiesFromTree(std::shared_ptr<Node> node, Body b, std::vector<Body> *TB);


//~ BodySolver

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
