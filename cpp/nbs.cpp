#include "nbs.h"

namespace N_body_simulation {

inline void BodySolver::ComputeAcceleration(std::vector<Body> *W) {
	//~ Calculates acceleration due to gravity forces between bodies W
#pragma omp parallel for
	for (id n=0; n<W->size(); ++n) {
	
		Body &b = (*W)[n];
		b.a_x = 0;
		b.a_y = 0;
	
		for (id i=0; i<W->size(); ++i) {
			const Body &a = (*W)[i];
			
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

void BodySolver::WritePositionToFile(std::fstream &x_fst, std::fstream &y_fst) {
	for (Body b : B) {
		x_fst<<b.x<<" ";
		y_fst<<b.y<<" ";
	}
	x_fst<<std::endl;
	y_fst<<std::endl;
}

void BodySolver::SolveTimeEvolution(unsigned int num_steps) {
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

void Euler::Advance() {
	// Estimates the positions and velocities after dt using a standard
	// Euler approach.
	
	ComputeAcceleration(&B);
	
	//~ Calculate and overwrite velocities
	for (Body &b : B) {
		b.x += dt * b.v_x;
		b.y += dt * b.v_y;
		
		b.v_x += dt * b.a_x;
		b.v_y += dt * b.a_y;
	}
}

void EulerAcc::Advance() {
	// Estimates the positions and velocities after dt using a standard
	// Euler approach. Version including acceleration.
	double kt = 0.5*dt*dt;
	
	ComputeAcceleration(&B);
	
	//~ Calculate and overwrite velocities
	for (Body &b : B) {
        b.x += b.a_x*kt + b.v_x*dt;
        b.y += b.a_y*kt + b.v_y*dt;

		b.v_x += dt * b.a_x;
		b.v_y += dt * b.a_y;
	}
}

void EulerStabilized::Advance() {
	// Estimates the positions and velocities after dt using a stabilized
	// Euler method, via a "semi"-implicit position update.
	
	ComputeAcceleration(&B);
	
	for (Body &b : B) {
		b.v_x += dt * b.a_x;
		b.v_y += dt * b.a_y;

        b.x += b.v_x*dt;
        b.y += b.v_y*dt;
	}
}

void EulerImproved::Advance() {
    // Perform one time step with the improved Euler method, also referred to
    // as Heun's method and 2-stage Runge-Kutta method.

    // 0.) update acceleration for current positions at time t0
	ComputeAcceleration(&B);
    // copy the set of bodies
    B_prec = B;
    // 1.a) perform 1 Euler step to get a position prediction at time t1
    for (Body &b : B_prec) {
        b.x += dt*b.v_x;
        b.y += dt*b.v_y;
        b.v_x += dt*b.a_x;
        b.v_y += dt*b.a_y;
    }
    // 1.b) acceleration prediction at time t1
    ComputeAcceleration(&B_prec);
    // 2.a) compute mean acceleration and
    //   b) use it to perform correction time step 
	for (unsigned int i=0; i<B.size(); ++i) {
        B[i].x += 0.5*dt*(B[i].v_x + B_prec[i].v_x);
        B[i].y += 0.5*dt*(B[i].v_y + B_prec[i].v_y);
		
        B[i].v_x += 0.5*dt*(B[i].a_x + B_prec[i].a_x);
        B[i].v_y += 0.5*dt*(B[i].a_y + B_prec[i].a_y);
	}

}
void EulerAccImproved::Advance() {
    // Perform one time step with the improved Euler method, also referred to
    // as Heun's method and 2-stage Runge-Kutta method.

	double kt = 0.5*dt*dt;
    // 0.) update acceleration for current positions at time t0
	ComputeAcceleration(&B);
    // copy the set of bodies
    B_prec = B;
    // 1.a) perform 1 Euler step to get a position prediction at time t1
    for (Body &b : B_prec) {
        b.x += b.a_x*kt + b.v_x*dt;
        b.y += b.a_y*kt + b.v_y*dt;
    }
    // 1.b) acceleration prediction at time t1
    ComputeAcceleration(&B_prec);
    // 2.a) compute mean acceleration and
    //   b) use it to perform correction time step 
	for (unsigned int i=0; i<B.size(); ++i) {
        // can I use (Body a, Body b:B, B_prec) construct?
        double acc_mean_x = 0.5*(B[i].a_x + B_prec[i].a_x);
        double acc_mean_y = 0.5*(B[i].a_y + B_prec[i].a_y);
        // or overwrite B[i].a_x with acc_mean_x ... does not matter
		B[i].x += acc_mean_x*kt + B[i].v_x*dt;
		B[i].y += acc_mean_y*kt + B[i].v_y*dt;
		
        B[i].v_x += acc_mean_x*dt;
        B[i].v_y += acc_mean_y*dt;
	}

}

}
