#include "NBS.h"

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
	//~ double length = 4;
	//~ Body b0 {-0.5, -7.5, 0, 0, 0, 0, 0};
	//~ Body b1 {-1, -3, 0, 0, 0, 0, 1};
	//~ Body b2 {1, -1, 0, 0, 0, 0, 2};
	//~ Body b3 {1, 1, 0, 0, 0, 0, 3};
	//~ Body b4 {0.5, 0.5, 0, 0, 0, 0, 4};
	//~ Body b5 {2, 2, 0, 0, 0, 0, 5};
	//~ Body b6 {-3, -2, 0, 0, 0, 0, 6};
	//~ Body b7 {-3.5, -2.5, 0, 0, 0, 0, 7};
	//~ Body b8 {-3, -3, 0, 0, 0, 0, 8};
	
	//~ std::vector<Body> B {b1,b2,b3,b4,b5,b6,b7,b8};
	
	//~ std::shared_ptr<Node> root = BuildTree(B,length);
	//~ print_tree(root);
	
	//~ std::vector<Body> TB;
	
	
	//~ GetBodiesFromTree(root, b4, &TB);
	//~ for (Body b : TB) std::cout<<b.m<<std::endl;
	
	int N = 1e3;
	double num_steps = 5e3;
	std::vector<Body> B = GenerateClash(N);
	
	Euler_std elstd {B};
	elstd.solve_for(num_steps);


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

