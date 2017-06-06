#ifndef PROPERTIES_HPP_
#define PROPERTIES_HPP_

#include <string>
#include "src/utils/Point.hpp"

struct MainProperties
{
	// Dimensions
	double x_dim, t_dim, p_dim;
	// Other params
	double visc;
	double kx, kz;
	double perm;
	double rate;

	Point sizes;
	// Other props
	static const double porosity;
	static const double compressibility;
	// Boundary fluxes
	double fx1, fx2, fy1, fy2, fz1, fz2;
};
struct WellGeomProperties
{
	// Spacial params
	Point rc;	// middle point of well
	double length;
	double alpha;
	Point r1, r2;
	double rw;
	// Number of segments
	int seg_num;
};
struct FracGeomProperties
{
};
struct SummatorProperties
{
	// Numbers
	int K;
	int M, N, L;
	int I;
	// Integral division limit
	double xi_c;
	// Observation point
	Point r_obs;
};

#endif /* PROPERTIES_HPP_ */