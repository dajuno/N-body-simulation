#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <queue>
#include <unordered_map>
#include <map>
#include <boost/functional/hash.hpp>
#include <iomanip>      // std::setw
#include <limits>
#include <bitset>         // std::bitset
#include <algorithm>    // std::fi

struct Body {
	double x_loc;
	double y_loc;
	
	double x_vel;
	double y_vel;
	
	double x_acc;
	double y_acc;
};

class NBS {
public:
	NBS (vector<Body> BB0, double hh) : B0 {BB0}, h {hh} {}
	
	void solve();
	
private
	double step_size = 1; // time step size
	vector<Body> B0; // Current
	vector<Body> B1; // Next step
};
