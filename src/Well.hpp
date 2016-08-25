#ifndef WELL_HPP_
#define WELL_HPP_

#include <iostream>
#include <vector>
#include <algorithm>
#include <new>
#include <cmath>

#define EQUALITY_TOLERANCE 1.E-8
#define BAR 1.E+5

struct Point
{
	union
	{
		double coords [3];
		struct
		{
			double x;	double y;	double z;
		};
	};
	
	Point() {};
	Point(const double _x, const double _y, const double _z) : x(_x), y(_y), z(_z) { };
	
	Point(const Point& a)
	{
		(*this) = a;
	};
	Point& operator=(const Point& rhs)
	{
		x = rhs.x,	y = rhs.y,	z = rhs.z;
		return *this;
	};
	Point& operator/=(const double k)
	{
		x /= k;	y /= k;	z /= k;
		return *this;
	};
};

inline std::ostream& operator<<(std::ostream& os, const Point& a)
{
	os << a.x << " " << a.y << " " << a.z << std::endl;
	return os;
}

inline bool operator==(const Point& a1, const Point& a2)
{
	if( (fabs(a2.x - a1.x) > EQUALITY_TOLERANCE) || 
		(fabs(a2.y - a1.y) > EQUALITY_TOLERANCE) || 
		(fabs(a2.z - a1.z) > EQUALITY_TOLERANCE) )
		return false;
	else
		return true;
};	
inline Point operator-(const Point& rhs)
{
	return Point(-rhs.x, -rhs.y, -rhs.z);
};
inline Point operator-(const Point& a1, const Point& a2)
{
	return Point(a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
};
inline Point operator+(const Point& rhs)
{
	return Point(rhs.x, rhs.y, rhs.z);
};
inline Point operator+(const Point& a1, const Point& a2)
{
	return Point(a1.x + a2.x, a1.y + a2.y, a1.z + a2.z);
};
inline Point operator*(const Point& a1, double k)
{
	return Point(a1.x * k, a1.y * k, a1.z * k);
};
inline Point operator*(double k, const Point& a1)
{
	return a1 * k;
};
inline Point operator/(const Point& a1, double k)
{
	return Point(a1.x / k, a1.y / k, a1.z / k);
};
inline double operator*(const Point& a1, const Point& a2)
{
	return a1.x * a2.x + a1.y * a2.y + a1.z * a2.z;
};

struct WellSegment
{
	const Point r1;
	const Point r2;
	const Point r_bhp;
	
	double tau1, tau2;
	double length;
	double pres;
	double pres2D;
	double pres3D;
	double rate;
	
	WellSegment(const Point& _r1, const Point& _r2, const Point _r_bhp, double _tau1, double _tau2) : r1(_r1), r2(_r2), r_bhp(_r_bhp), tau1(_tau1), tau2(_tau2)
	{
		length = sqrt((r2 - r1) * (r2 -r1));
		pres = pres2D = pres3D = rate = 0.0;
	};
};

struct Parameters
{	
	// Dimensions
	double x_dim, t_dim, p_dim;
	
	// Spacial params
	Point sizes;
	Point rc;	// middle point of well
	double length;
	double alpha;
	Point r1, r2;
	
	double rw;
	
	// Other params
	double visc;
	double kx, kz;
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
	
	// Observation point
	Point r_obs;
};

class Well
{
	protected:
		
		const Point r1;
		const Point r2;
		const int num;
		const double r_w;
	
		double alpha;
		double length;
		double rate;
		
	public:
		Well(const Point& _r1, const Point& _r2, const int _num, const double _r_w);
		~Well();
		
		void setRate(double _rate);
		void setUniformRate();
		
		std::vector<WellSegment> segs;
		
		double pres_av;
		double pres_dev;
		
		void printRates(const Parameters* props);
		void writeRates(const Parameters* props);

};

#endif /* WELL_HPP_ */

