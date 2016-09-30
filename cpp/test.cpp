#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>

//~ const double G=6.674e-11; //Gravitational constant
const double G=1; //Gravitational constant

struct Body {
	double x_pos;
	double y_pos;
	
	double x_vel;
	double y_vel;
	
	double x_acc;
	double y_acc;
	
	double mass; 
};
double cubed_dist(Body i, Body j) {
	return pow(
				pow(i.x_pos - j.x_pos,1.5)
				+pow(i.y_pos - j.y_pos,1.5),
				2);
}

class NBS {
	using id=std::vector<Body>::size_type;
public:
	NBS (std::vector<Body> BB, double time_step_sz, double softening_const);
	void solve(unsigned int num_steps);
	
private:
	void solve_next(); // Solves for the next step
	void solve_acc(id);
	void solve_vel(id);
	void solve_pos(id);

	std::vector<Body> B;
	double time_step_size = 1;
	double softening_constant = 0.001;
};

NBS::NBS (std::vector<Body> BB, double tt, double ss) 
	: B {BB}, time_step_size {tt}, softening_constant {ss} {
		}

void NBS::solve_acc(id n) {
	B[n].x_acc=0;
	B[n].y_acc=0;
	
	for (id i=0; i<B.size(); ++i) {
		if (i==n) {continue;}
		double cube_dist = std::max(softening_constant,cubed_dist(B[i],B[n]));
		B[n].x_acc += G*B[i].x_pos*B[i].mass / cube_dist;
		B[n].y_acc += G*B[i].y_pos*B[i].mass / cube_dist;
	}
}
void NBS::solve_vel(id n) {
	B[n].x_vel += time_step_size*B[n].x_acc;
	B[n].y_vel += time_step_size*B[n].y_acc;
}
void NBS::solve_pos(id n) {
	B[n].x_pos += B[n].x_acc *time_step_size*time_step_size/2 + B[n].x_vel *time_step_size;
	B[n].y_pos += B[n].y_acc *time_step_size*time_step_size/2 + B[n].y_vel *time_step_size;
}
void NBS::solve_next() {
	for (id i=0; i<B.size(); ++i) solve_acc(i);
	for (id i=0; i<B.size(); ++i) {
		solve_vel(i);
		solve_pos(i);
	}
}

void NBS::solve(unsigned int num_steps) {
	for (unsigned int s=0; s<num_steps; ++s) {
		std::cout<<"Step "<<s<<std::endl;
		std::cout<<"Pos: ";
		for (Body b : B) {
			std::cout<<"("<<b.x_pos<<","<<b.y_pos<<") ";
		}
		std::cout<<std::endl;
		std::cout<<"Vel: ";
		for (Body b : B) {
			std::cout<<"("<<b.x_vel<<","<<b.y_vel<<") ";
		}
		std::cout<<std::endl;
		std::cout<<"Acc: ";
		for (Body b : B) {
			std::cout<<"("<<b.x_acc<<","<<b.y_acc<<") ";
		}
		std::cout<<std::endl;
		
		solve_next();
	}
	
	std::fstream x_fst("x");
	std::fstream y_fst("y");
		
	for (unsigned int s=0; s<num_steps; ++s) {
		for (Body b : B) {
			x_fst<<b.x_pos<<" ";
			y_fst<<b.y_pos<<" ";
		}
		x_fst<<std::endl;
		y_fst<<std::endl;
		
		solve_next();
	}
}

std::vector<Body> BuildBodies(std::vector<Body>::size_type N, double r) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_radius(0, r);
    std::uniform_real_distribution<> dis_angle(0, M_PI);

	std::vector<Body> B (N);
	for (std::vector<Body>::size_type i=0; i<N; ++i) {
		//~ B[i].x_pos = dis_radius(gen) * cos(dis_angle(gen));
		//~ B[i].y_pos = dis_radius(gen) * sin(dis_angle(gen));
		B[i].x_pos = (i+1)*r/N * cos((i+1)*3.1/N);
		B[i].y_pos = (i+1)*r/N * sin((i+1)*3.1/N);
		
		B[i].x_vel = 0;
		B[i].y_vel = 0;
		
		B[i].x_acc = 0;
		B[i].y_acc = 0;
		
		B[i].mass = 100;
	}
	
	return B;
}


int main() {
	try
	{
	std::vector<Body> B = BuildBodies(5,1);
	NBS nbs {B,10,1e-10};
	nbs.solve(4);
		
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
