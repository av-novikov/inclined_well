#ifndef INCLINEDSUM_HPP_
#define INCLINEDSUM_HPP_

#include <cmath>

#include "src/Well.hpp"

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

class InclinedSum
{
protected:
	const Parameters* props;
	const Well* well;
public:
	InclinedSum(const Parameters* _props, const Well* _well);
	~InclinedSum();
		
	virtual double getPres(const Point& r);
	virtual double get2D(const Point& r);
	virtual double get3D(const Point& r);
};

#endif /* INCLINEDSUM_HPP_ */
