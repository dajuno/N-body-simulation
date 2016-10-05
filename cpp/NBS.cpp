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
		double dx = B[i].x - B[n].x;
		double dy = B[i].y - B[n].y;
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
    // standard first order Euler method
	
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

void NBS::forces(std::vector<Force>* F) {
	// Calculates forces without overwriting

    double dx, dy, d, m;

#pragma omp parallel for
    for (id n=0; n<B.size(); ++n) {
        (*F)[n].x = 0;
        (*F)[n].y = 0;
        for (id i=0; i<B.size(); ++i) {
            dx = B[i].x - B[n].x;
            dy = B[i].y - B[n].y;
            d = sqrt(dx*dx + dy*dy + softening_constant);
            m = B[i].m;
            
            (*F)[n].x += dx*m / (d*d*d);
            (*F)[n].y += dy*m / (d*d*d);
        }
        
        (*F)[n].x *= gravitational_constant;
        (*F)[n].y *= gravitational_constant;
    }

}

void NBS::euler_improved() {
    // compute one time step with the Euler corrected integration scheme
    // cf. Aarseth (2003): Gravitational N-Body-Simulations 
    //
    // better not to do these declarations at every time step!?
    int N = B.size();
    std::vector<Force> F1(N);
    std::vector<Force> F2(N);
    std::vector<double> x0(N);
    std::vector<double> y0(N);
    // std::vector<double> vx0(N);
    // std::vector<double> vy0(N);

    double kt = 0.5*dt*dt;
    double Fmean_x, Fmean_y;

    // 1. first order step
    forces(&F1);

	for (id n=0; n<B.size(); ++n) {
        // save positions and velocities at beginning of time interval
        x0[n] = B[n].x;
        y0[n] = B[n].y;
        // vx0[n] = B[n].xv;
        // vy0[n] = B[n].yv;

		// B[n].xv += dt*F1[n].x;
		// B[n].yv += dt*F1[n].y;
		
        B[n].x += F1[n].x * kt + B[n].xv*dt;
        B[n].y += F1[n].y * kt + B[n].yv*dt;
	}

    // 2. correction with mean of forces F2=F(t1), F1=F(t0)
    forces(&F2);

	for (id n=0; n < B.size(); ++n) {
        Fmean_x = 0.5*(F1[n].x + F2[n].x);
        Fmean_y = 0.5*(F1[n].y + F2[n].y);
		
        B[n].x = Fmean_x * kt + B[n].xv *dt + x0[n];
        B[n].y = Fmean_y * kt + B[n].yv *dt + y0[n];

		B[n].xv += dt*Fmean_x;
		B[n].yv += dt*Fmean_y;
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
    // solve_next();
    euler_improved();
    std::cout << s << std::endl;
	}
}

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
		
		B[i].xv = -sqrt(num/r) * sin(angle);
		B[i].yv = sqrt(num/r) * cos(angle);
		
		B[i].xa = 0;
		B[i].ya = 0;
		
		B[i].m = 10*solar_mass;
	}
	return B;
}

}
