#ifndef WELLFLOW_HPP_
#define WELLFLOW_HPP_

#include "src/Well.hpp"

#include <cmath>

#define BAR 1.E+5

struct Parameters
{	
	// Dimensions
	double x_dim, t_dim, p_dim;
	
	// Spacial params
	Point sizes;
	Point r1;	// at z = 0
	double alpha;
	Point r2;
	
	double rw;
	
	// Other params
	double visc;
	double perm;
	double rate;
	
	// Numbers
	int K;	
	int M, N, L;	
	int I;
	
	// Grid sizes
	int nx, ny, nz;
	
	// Integral division limit
	double xi_c;
};

class WellFlow
{
	protected:
		Parameters props;
		Well* well;
		
		void findRateDistribution();
		
		double highIntegral2D(const Point& r, double xi_c);
		double lowIntegral2D(const Point& r, double xi_c);
		double calc3D(const Point& r);
		
	public:
		WellFlow(const Parameters& _props);
		~WellFlow();
		
		virtual double calcPressure(const Point& r, bool debug);
};

#endif /* WELLFLOW_HPP_ */
