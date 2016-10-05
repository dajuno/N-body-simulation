#include "NBS.h"

namespace N_body_simulation {

inline void BS::ComputeAcceleration(std::vector<Body> *R, std::vector<Body> *W) {
	//~ Calculates acceleration using data of R and writes it to W
	#pragma omp parallel for
	for (id n=0; n<W->size(); ++n) {
	
		Body &b = (*W)[n];
		b.a_x=0;
		b.a_y=0;
	
		for (id i=0; i<R->size(); ++i) {
			Body &a = (*R)[i];
			
			double dx = a.x - b.x;
			double dy = a.y - b.y;
			double d = sqrt(dx*dx + dy*dy + softening_constant);
			
			b.a_x += dx*a.m / (d*d*d);
			b.a_y += dy*a.m / (d*d*d);
		}
		
		b.a_x *= gravitational_constant;
		b.a_y *= gravitational_constant;
	}
}

void BS::WritePositionToFile(std::fstream &x_fst, std::fstream &y_fst) {
	for (Body b : B) {
		x_fst<<b.x<<" ";
		y_fst<<b.y<<" ";
	}
	x_fst<<std::endl;
	y_fst<<std::endl;
}

void BS::solve_for(unsigned int num_steps) {
	//~ Advances the system num_steps times and, after each step,
	//~ writes a line with the	current positions into the files "x" 
	//~ and "y"
	
	// Empty and open files
	std::fstream x_fst("x", std::fstream::out | std::fstream::trunc);
	std::fstream y_fst("y", std::fstream::out | std::fstream::trunc);
	
	for (unsigned int s=0; s<num_steps; ++s) {
		WritePositionToFile(x_fst, y_fst);
		Advance();
	}
}

void Euler_std::Advance() {
	// Estimates the positions and velocities after dt using a standard
	// Euler approach
	
	ComputeAcceleration(&B, &B);
	
	//~ Calculate and overwrite velocities
	for (Body &b : B) {
		b.x += b.v_x * dt;
		b.y += b.v_y * dt;
		
		b.v_x += dt * b.a_x;
		b.v_y += dt * b.a_y;
	}
}

}
