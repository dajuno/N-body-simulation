#include "nbs.h"

namespace N_body_simulation {

inline double eucl_dist(Body a, Body b) {
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	return sqrt(dx*dx + dy*dy);
}

bool BodyIsInQuadrant(Body b, Body q, double l) {
	return ( b.x<q.x+l && b.x>=q.x-l ) 
		&& ( b.y<q.y+l && b.y>=q.y-l );
}

void insert_in_tree(std::shared_ptr<Node> node, Body b) {	
	// First add weight to the current quadrant
	node->q.m += b.m;
	
	// If current node is empty
	if (node->empty) {
		node->b = b;
		node->empty = false;
	} else {
		// Build possible sub quadrants
		double x = node->q.x;
		double y = node->q.y;
		double hl = node->l/2;
		
		Body NE_sq {x+hl, y+hl, 0, 0, 0, 0, 0};
		Body NW_sq {x-hl, y+hl, 0, 0, 0, 0, 0};
		Body SE_sq {x+hl, y-hl, 0, 0, 0, 0, 0};
		Body SW_sq {x-hl, y-hl, 0, 0, 0, 0, 0};
		
		//~ Check in which sub quadrant b should go.
		//~ If there is no node, create one and place the quadrant 
		//~ inside.
		//~ Put b insided by recursing.
	
		if ( BodyIsInQuadrant(b,NE_sq,hl) ) {
			if (!node->NE) node->NE = std::shared_ptr<Node> ( new Node(NE_sq,hl) );
			insert_in_tree(node->NE,b);
		}
		
		if ( BodyIsInQuadrant(b,NW_sq,hl) ) {
			if (!node->NW) node->NW = std::shared_ptr<Node> ( new Node(NW_sq,hl) );
			insert_in_tree(node->NW,b);
		}
		
		if ( BodyIsInQuadrant(b,SE_sq,hl) ) {
			if (!node->SE) node->SE = std::shared_ptr<Node> ( new Node(SE_sq,hl) );
			insert_in_tree(node->SE,b);
		}
		
		if ( BodyIsInQuadrant(b,SW_sq,hl) ) {
			if (!node->SW) node->SW = std::shared_ptr<Node> ( new Node(SW_sq,hl) );
			insert_in_tree(node->SW,b);
		}
	}
}

std::shared_ptr<Node> BuildTree(std::vector<Body> *B, double universe_length) {
	// Builds a QuadTree for bodies inside (!) the Quadrant centered
	//~ at the origin and of side length 2l.
	Body origin;
	
	std::shared_ptr<Node> root = 
		std::shared_ptr<Node> ( new Node {origin, universe_length} );
	
	for (Body b : (*B)) {
		if ( BodyIsInQuadrant(b,origin,universe_length) ) 
			insert_in_tree(root,b);
	}
		
	return root;
}

void print_tree(std::shared_ptr<Node> node) {
	if (node) {
		std::cout<<"Quadrant pos: ("<<node->q.x<<","<<node->q.y<<")";
		std::cout<<", Length pos: "<<node->l;
		std::cout<<", Quadrant mass: "<<node->q.m;
		std::cout<<", Body mass: "<<node->b.m;
		std::cout<<std::endl;
		
		print_tree(node->NE);
		print_tree(node->NW);
		print_tree(node->SE);
		print_tree(node->SW);
	}
}

void GetBodiesFromTree(std::shared_ptr<Node> node, Body b, std::vector<Body> *TB) {
		if (node) {
			// If node is a leaf, push back its body
			bool is_leaf = !node->NE && !node->NW && !node->SE && !node->SW;
			
			if (is_leaf) {
				(*TB).push_back(node->b);
				//~ std::cout<<"Adding body: "<<"("<<node->b.x<<","<<node->b.y<<") of mass "<<node->b.m<<std::endl;
			}
			// If node is not a leaf but far away, push back its center
			else if (eucl_dist(b,node->q) > 2*(node->l)) {
				//~ std::cout<<"Adding quadrant: "<<"("<<node->q.x<<","<<node->q.y<<") of mass "<<node->q.m<<std::endl;
				(*TB).push_back(node->q);
			}
			// Otherwise recurse
			else { 
				(*TB).push_back(node->b);
				//~ std::cout<<"Adding body: "<<"("<<node->b.x<<","<<node->b.y<<") of mass "<<node->b.m<<std::endl;
				GetBodiesFromTree(node->NE, b, TB);
				GetBodiesFromTree(node->NW, b, TB);
				GetBodiesFromTree(node->SE, b, TB);
				GetBodiesFromTree(node->SW, b, TB);
			}
		}
	}

//~ inline void BodySolver::ComputeAcceleration(std::vector<Body> *B) {
	//~ // BRUTE FORCE
	//~ Calculates acceleration due to gravity forces between bodies B
//~ #pragma omp parallel for
	//~ for (id n=0; n<B->size(); ++n) {
	
		//~ Body &b = (*B)[n];
		//~ b.a_x = 0;
		//~ b.a_y = 0;
	
		//~ for (id i=0; i<B->size(); ++i) {
			//~ const Body &a = (*B)[i];
			
			//~ double dx = a.x - b.x;
			//~ double dy = a.y - b.y;
			//~ double d = sqrt(dx*dx + dy*dy + softening_constant);
			
			//~ b.a_x += dx*a.m / (d*d*d);
			//~ b.a_y += dy*a.m / (d*d*d);
		//~ }
		
		//~ b.a_x *= gravitational_constant;
		//~ b.a_y *= gravitational_constant;
	//~ }
//~ }

inline void BodySolver::ComputeAcceleration(std::vector<Body> *B) {
	// TREE
	//~ Calculates acceleration due to gravity forces between bodies B
	std::shared_ptr<Node> root = BuildTree(B,5e18);
#pragma omp parallel for
	for (id n=0; n<B->size(); ++n) {
	
		Body &b = (*B)[n];
		b.a_x = 0;
		b.a_y = 0;

		std::vector<Body> TB;
		GetBodiesFromTree(root, b, &TB);
	
		for (id i=0; i<TB.size(); ++i) {
			const Body &a = TB[i];
			
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

void BodySolver::WriteEnergyToFile(std::fstream &e_fst, const double E) {
    e_fst << E << std::endl;
}


void BodySolver::SolveTimeEvolution(unsigned int num_steps) {
	//~ Advances the system num_steps times and, after each step,
	//~ writes a line with the	current positions into the files "x" 
	//~ and "y"
	
	// Empty and open files
	std::fstream x_fst("x", std::fstream::out | std::fstream::trunc);
	std::fstream y_fst("y", std::fstream::out | std::fstream::trunc);
	std::fstream e_fst("energy", std::fstream::out | std::fstream::trunc);
	
	for (unsigned int s=0; s<num_steps; ++s) {
		WritePositionToFile(x_fst, y_fst);
        WriteEnergyToFile(e_fst, ComputeTotalEnergy());

		Advance();
	}
}

double BodySolver::ComputeTotalEnergy() {
    // compute total energy of the N-body system at the current time
    // this could be merged with ComputeAcceleration and enabled with an option
    // switch
    double E_tot = 0;
    double E_kin = 0;
    double E_pot = 0;

#pragma omp parallel for    
    for (unsigned int i=0; i<B.size(); ++i) {
        E_kin += 0.5*B[i].m*(B[i].v_x*B[i].v_x + B[i].v_y*B[i].v_y);

        double E_pot_loc = 0;
        for (unsigned int k=0; k<B.size(); ++k) {
            if (i != k) {
                double dx = B[i].x - B[k].x;
                double dy = B[i].y - B[k].y;
                double d = sqrt(dx*dx + dy*dy + softening_constant);
                E_pot_loc += B[k].m/d;
            }
        }
        E_pot += 0.5*B[i].m*gravitational_constant*E_pot_loc;
    }

    E_tot = E_kin + E_pot;
    return E_tot;
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
