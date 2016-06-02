#ifndef WELL_HPP_
#define WELL_HPP_

#include <iostream>
#include <vector>
#include <algorithm>
#include <new>
#include <cmath>

#define EQUALITY_TOLERANCE 1.E-6

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
	
	double length;
	double pres;
	double rate;
	
	WellSegment(const Point& _r1, const Point& _r2) : r1(_r1), r2(_r2)
	{
		length = sqrt((r2 - r1) * (r2 -r1));
		pres = rate = 0.0;
	};
};

class Well
{
	protected:
		
		const Point r1;
		const Point r2;
	
		double length;
		double pres_av;
		double rate;
		
		const int num;
		
		void uniformRate(WellSegment& seg);
		
	public:
		Well(const Point& _r1, const Point& _r2, const int _num);
		~Well();
		
		void setRate(double _rate);
		void setUniformRate();
		
		std::vector<WellSegment> segs;
};

#endif /* WELL_HPP_ */

