#include "nbs.h"

namespace N_body_simulation {

std::vector<Body> BuildBodies(std::vector<Body>::size_type N) {
	//~ Builds a vector of starting Bodies which are sampled within in
	//~ circle of radius r and have apropriate starting velocities.
	//~ A large mass is put into the origin, to simulate dark matter :-)
	
	std::vector<Body> B (N+1);
		
	// Center with large mass
	Body center {0,0,0,0,0,0,1e6*solar_mass};
	B[0] = center;
	
	// Random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
	double radius = 1e18;
    std::uniform_real_distribution<> dis_radius(0, radius);
    std::uniform_real_distribution<> dis_angle(0, 2*M_PI);

	// Sample body positions and calculate inial velocities (physics?)
	for (std::vector<Body>::size_type i=1; i<N+1; ++i) {
		double r = dis_radius(gen);
		double angle = dis_angle(gen);
		double num = gravitational_constant*1e6*solar_mass;
		
		B[i].x = r * cos(angle);
		B[i].y = r * sin(angle);
		
		B[i].v_x = -sqrt(num/r) * sin(angle);
		B[i].v_y = sqrt(num/r) * cos(angle);
		
		B[i].a_x = 0;
		B[i].a_y = 0;
		
		B[i].m = 10*solar_mass;
	}
	return B;
}

void add_pos_bodies(std::vector<Body> &B, double dx, double dy) {
	for (Body& b : B) {
		b.x += dx;
		b.y += dy;
	}
}

void add_vel_bodies(std::vector<Body> &B, double dxv, double dyv) {
	for (Body& b : B) {
		b.v_x += dxv;
		b.v_y += dyv;
	}
}

void flip_vel_bodies(std::vector<Body> &B) {
	for (Body& b : B) {
		b.v_x *= -1;
		b.v_y *= -1;
	}
}

std::vector<Body> GenerateClash(int N) {
	std::vector<Body> B1 = BuildBodies(N);
	std::vector<Body> B2 = BuildBodies(N);

	add_pos_bodies(B1, -5e17, -5e17);
	add_pos_bodies(B2, 5e17, 5e17);

	//~ flip_vel_bodies(B1);
	flip_vel_bodies(B2);

	add_vel_bodies(B1, 5e3, 0);
	add_vel_bodies(B2, -4e3, 0);



	for (Body b : B2) B1.push_back(b);

	return B1;
} 
}

using namespace N_body_simulation;

int main() {
	try
	{
	
	
	int N = 1e2;
	double num_steps = 5e3;
	std::vector<Body> B = GenerateClash(N);
	
    EulerImproved elstd {B};
	// EulerLegacy elstd {B};
    // Euler elstd {B};
	elstd.SolveTimeEvolution(num_steps);

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

