#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include "NBS.h"

namespace N_body_simulation {

inline void NBS::solve_acc(id n) {
	//~ Calculates and overwrites acceleration using data 
	//~ from current positions
	B[n].xa=0;
	B[n].ya=0;
	
	for (id i=0; i<B.size(); ++i) {
		double dx = B[i].x-B[n].x;
		double dy = B[i].y-B[n].y;
		double d = sqrt(dx*dx + dy*dy + softening_constant);
		double m = B[i].m;
		
		B[n].xa += dx*m / (d*d*d);
		B[n].ya += dy*m / (d*d*d);
	}
	
	B[n].xa *= gravitational_constant;
	B[n].ya *= gravitational_constant;
}

void NBS::solve_next() {
	// Estimates the positions and velocities after dt
	
	#pragma omp parallel for
	for (id i=0; i<B.size(); ++i) solve_acc(i);

	//~ Calculate and overwrite velocities and positions using data  
	//~ from current accelerations.
	#pragma omp parallel for
	for (id n=0; n<B.size(); ++n) {
		B[n].xv += dt*B[n].xa;
		B[n].yv += dt*B[n].ya;		
		
		//~ B[n].x += B[n].xa *dt*dt/2 + B[n].xv *dt;
		//~ B[n].y += B[n].ya *dt*dt/2 + B[n].yv *dt;
		B[n].x += B[n].xv *dt;
		B[n].y += B[n].yv *dt;
	}
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
			x_fst<<b.x<<" ";
			y_fst<<b.y<<" ";
		}
		x_fst<<std::endl;
		y_fst<<std::endl;
	
	// Solve the next iteration
	solve_next();
	}
}

std::vector<Body> BuildBodies(std::vector<Body>::size_type N) {
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
	double radius = 1e18;
    std::uniform_real_distribution<> dis_radius(0, radius);
    std::uniform_real_distribution<> dis_angle(0, 2*M_PI);

	// Sample body positions and calculate inial velocities (physics?)
	for (std::vector<Body>::size_type i=0; i<N; ++i) {
		double r = dis_radius(gen);
		double angle = dis_angle(gen);
		double num = gravitational_constant*1e6*solar_mass;
		
		Body b;
		b.x = r * cos(angle);
		b.y = r * sin(angle);
		
		b.xv = -sqrt(num/r) * sin(angle);
		b.yv = sqrt(num/r) * cos(angle);
		
		b.xa = 0;
		b.ya = 0;
		
		b.m = 10*solar_mass;
	
		B.push_back(b);
	}
	

	
	return B;
}

}
