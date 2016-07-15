#ifndef WELL_CUH_
#define WELL_CUH_

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#define EQUALITY_TOLERANCE 1.E-8
#define BAR 1.E+5

struct Point
{
	union
	{
		double coords[3];
		struct
		{
			double x;	double y;	double z;
		};
	};

	__host__ __device__ Point() {};
	__host__ __device__ Point(const double _x, const double _y, const double _z) : x(_x), y(_y), z(_z) { };

	__host__ __device__ Point(const Point& a)
	{
		(*this) = a;
	};
	__host__ __device__ Point& operator=(const Point& rhs)
	{
		x = rhs.x, y = rhs.y, z = rhs.z;
		return *this;
	};
	__host__ __device__ Point& operator/=(const double k)
	{
		x /= k;	y /= k;	z /= k;
		return *this;
	};
};

__host__ __device__ inline bool operator==(const Point& a1, const Point& a2)
{
	if ((fabs(a2.x - a1.x) > EQUALITY_TOLERANCE) ||
		(fabs(a2.y - a1.y) > EQUALITY_TOLERANCE) ||
		(fabs(a2.z - a1.z) > EQUALITY_TOLERANCE))
		return false;
	else
		return true;
};
__host__ __device__ inline Point operator-(const Point& rhs)
{
	return Point(-rhs.x, -rhs.y, -rhs.z);
};
__host__ __device__ inline Point operator-(const Point& a1, const Point& a2)
{
	return Point(a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
};
__host__ __device__ inline Point operator+(const Point& rhs)
{
	return Point(rhs.x, rhs.y, rhs.z);
};
__host__ __device__ inline Point operator+(const Point& a1, const Point& a2)
{
	return Point(a1.x + a2.x, a1.y + a2.y, a1.z + a2.z);
};
__host__ __device__ inline Point operator*(const Point& a1, double k)
{
	return Point(a1.x * k, a1.y * k, a1.z * k);
};
__host__ __device__ inline Point operator*(double k, const Point& a1)
{
	return a1 * k;
};
__host__ __device__ inline Point operator/(const Point& a1, double k)
{
	return Point(a1.x / k, a1.y / k, a1.z / k);
};
__host__ __device__ inline double operator*(const Point& a1, const Point& a2)
{
	return a1.x * a2.x + a1.y * a2.y + a1.z * a2.z;
};

struct WellSegment
{
	Point r1;
	Point r2;
	Point r_bhp;

	double length;
	double pres;
	double pres2D;
	double pres3D;
	double rate;

	__host__ __device__ WellSegment(const Point& _r1, const Point& _r2, const Point _r_bhp) : r1(_r1), r2(_r2), r_bhp(_r_bhp)
	{
		length = sqrt((r2 - r1) * (r2 - r1));
		pres = pres2D = pres3D = rate = 0.0;
	};

	__host__ __device__ WellSegment& operator=(const WellSegment& rhs) = default;
	/*{
		r1 = rhs.r1;
	};*/
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

	thrust::host_vector<WellSegment> segs;

	double pres_av;
	double pres_dev;

	void printRates(const Parameters* props);
	void writeRates(const Parameters* props);

};

class Well_Device
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
	__host__ __device__ Well_Device(const Point& _r1, const Point& _r2, const int _num, const double _r_w);
	__host__ __device__ ~Well_Device();

	thrust::device_vector<WellSegment> segs;

	double pres_av;
	double pres_dev;
};

#endif /* WELL_CUH_ */