#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include "NBS.h"

namespace N_body_simulation {

void NBS::solve_acc(id n) {
	//~ Calculates and overwrites acceleration using data 
	//~ from current positions
	B[n].x_acc=0;
	B[n].y_acc=0;
	
	for (id i=0; i<B.size(); ++i) {
		double dx = B[i].x_pos-B[n].x_pos;
		double dy = B[i].y_pos-B[n].y_pos;
		double d = sqrt(dx*dx + dy*dy)+softening_constant;
		double m = B[i].mass;
		
		B[n].x_acc += dx*m / (d*d*d);
		B[n].y_acc += dy*m / (d*d*d);
	}
	
	B[n].x_acc *= gravitational_constant;
	B[n].y_acc *= gravitational_constant;
}

void NBS::solve_vel(id n) {
	//~ Calculates and overwrites velocities using data from current 
	//~ accelerations
	B[n].x_vel += time_step_size*B[n].x_acc;
	B[n].y_vel += time_step_size*B[n].y_acc;
}

void NBS::solve_pos(id n) {
	//~ Calculates and overwrites positions using data from current  
	//~ accelerations and velocities
	B[n].x_pos += B[n].x_acc *time_step_size*time_step_size/2 
				+ B[n].x_vel *time_step_size;
	B[n].y_pos += B[n].y_acc *time_step_size*time_step_size/2 
				+ B[n].y_vel *time_step_size;
}

void NBS::solve_next() {
	// Solves the next iterations by putting the three operations 
	//~ together (in the right order!)
	for (id i=0; i<B.size(); ++i) solve_acc(i);
	for (id i=0; i<B.size(); ++i) solve_vel(i);
	for (id i=0; i<B.size(); ++i) solve_pos(i);
}

void NBS::solve(unsigned int num_steps) {
	//~ Advances the system num_steps times and, after each step,
	//~ writes a line with the	current positions into the files "x" 
	//~ and "y"
	
	// Open files
	std::fstream x_fst("x", std::fstream::out | std::fstream::trunc);
	std::fstream y_fst("y", std::fstream::out | std::fstream::trunc);
	
	// Write current positions to files
	for (unsigned int s=0; s<num_steps; ++s) {
		for (Body b : B) {
			x_fst<<b.x_pos<<" ";
			y_fst<<b.y_pos<<" ";
		}
		x_fst<<std::endl;
		y_fst<<std::endl;
	
	// Solve the next iteration
	solve_next();
	}
}

std::vector<Body> BuildBodies(std::vector<Body>::size_type N, double r) {
	//~ Builds a vector of starting Bodies which are sampled within in
	//~ circle of radius r and have apropriate starting velocities.
	//~ A large mass is put into the origin, to simulate dark matter :-)
	
	std::vector<Body> B (N);
	
	// Center with large mass
	Body center {0,0,0,0,0,0,1e6*solar_mass};
	B.push_back(center);
	
	// Random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_radius(0, r);
    std::uniform_real_distribution<> dis_angle(0, 2*M_PI);

	// Sample body positions and calculate inial velocities (physics?)
	for (std::vector<Body>::size_type i=0; i<N; ++i) {
		double radius = dis_radius(gen);
		double angle = dis_angle(gen);
		double num = gravitational_constant*1e6*solar_mass;
		
		Body b;
		b.x_pos = radius * cos(angle);
		b.y_pos = radius * sin(angle);
		
		b.x_vel = -sqrt(num/radius) * sin(angle);
		b.y_vel = sqrt(num/radius) * cos(angle);
		
		b.x_acc = 0;
		b.y_acc = 0;
		
		b.mass = 10*solar_mass;
	
		B.push_back(b);
	}
	

	
	return B;
}

}
