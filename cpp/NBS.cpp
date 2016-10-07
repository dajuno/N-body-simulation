#include "NBS.h"

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

std::shared_ptr<Node> BuildTree(std::vector<Body> B, double universe_length) {
	// Builds a QuadTree for bodies inside (!) the Quadrant centered
	//~ at the origin and of side length 2l.
	Body origin;
	
	std::shared_ptr<Node> root = 
		std::shared_ptr<Node> ( new Node {origin, universe_length} );
	
	for (Body b : B) {
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

void GetBodiesFromTree(std::shared_ptr<Node> node, Body b, std::vector<Body>* TB) {
		if (node) {
			// If node is a leaf, push back its body
			bool is_leaf = !node->NE && !node->NW && !node->SE && !node->SW;
			
			if (is_leaf) {
				(*TB).push_back(node->b);
				std::cout<<"Adding body: "<<"("<<node->b.x<<","<<node->b.y<<") of mass "<<node->b.m<<std::endl;
			}
			// If node is not a leaf but far away, push back its center
			else if (eucl_dist(b,node->q) > 2*(node->l)) {
				std::cout<<"Adding quadrant: "<<"("<<node->q.x<<","<<node->q.y<<") of mass "<<node->q.m<<std::endl;
				(*TB).push_back(node->q);
			}
			// Otherwise recurse
			else { 
				(*TB).push_back(node->b);
				std::cout<<"Adding body: "<<"("<<node->b.x<<","<<node->b.y<<") of mass "<<node->b.m<<std::endl;
				GetBodiesFromTree(node->NE, b, TB);
				GetBodiesFromTree(node->NW, b, TB);
				GetBodiesFromTree(node->SE, b, TB);
				GetBodiesFromTree(node->SW, b, TB);
			}
		}
	}



inline void BS::ComputeAcceleration(const std::vector<Body> &R, std::vector<Body> *W) {
	//~ Calculates acceleration using data of R and writes it to W
	#pragma omp parallel for
	for (id n=0; n<W->size(); ++n) {
	
		Body &b = (*W)[n];
		b.a_x=0;
		b.a_y=0;
	
		for (const Body &a : R) {			
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
	
	ComputeAcceleration(B, &B);
	
	//~ Calculate and overwrite velocities
	for (Body &b : B) {
		b.x += b.v_x * dt;
		b.y += b.v_y * dt;
		
		b.v_x += dt * b.a_x;
		b.v_y += dt * b.a_y;
	}
}

}
