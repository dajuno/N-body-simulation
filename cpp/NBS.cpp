#include "NBS.h"

namespace N_body_simulation {

inline void BS::ComputeAcceleration(std::vector<Body> *R, std::vector<Body> *W)
{
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

/* THIS FUNCTION SHOULD BE OBSOLETE ! CHECK XXX
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
*/

/* INTEGRATE THIS FUNCTION -> CLASS INHERITANCE
void NBS::euler_improved() {
    // compute one time step with the Euler corrected integration scheme
    // cf. Aarseth (2003): Gravitational N-Body-Simulations 
    //
    // better not to do these declarations at every time step!?
    std::vector<Force> F1(B.size());
    std::vector<Force> F2(B.size());
    std::vector<double> x0(B.size());
    std::vector<double> y0(B.size());
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
		
        B[n].x = Fmean_x*kt + B[n].xv *dt + x0[n];
        B[n].y = Fmean_y*kt + B[n].yv *dt + y0[n];

		B[n].xv += dt*Fmean_x;
		B[n].yv += dt*Fmean_y;
	}


}
*/

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
