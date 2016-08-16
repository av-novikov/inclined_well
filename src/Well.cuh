#ifndef WELL_CUH_
#define WELL_CUH_

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#define EQUALITY_TOLERANCE 1.E-8
#define BAR 1.E+5

template <class T>
struct Point
{
	union
	{
		T coords[3];
		struct
		{
			T x;	T y;	T z;
		};
	};

	__host__ __device__ Point() {};
	__host__ __device__ Point(const T _x, const T _y, const T _z) : x(_x), y(_y), z(_z) { };

	__host__ __device__ Point(const Point& a)
	{
		(*this) = a;
	};
	__host__ __device__ Point& operator=(const Point& rhs)
	{
		x = rhs.x, y = rhs.y, z = rhs.z;
		return *this;
	};
	__host__ __device__ Point& operator/=(const T k)
	{
		x /= k;	y /= k;	z /= k;
		return *this;
	};
};

template <class T>
__host__ __device__ inline bool operator==(const Point<T>& a1, const Point<T>& a2)
{
	if ((fabs(a2.x - a1.x) > EQUALITY_TOLERANCE) ||
		(fabs(a2.y - a1.y) > EQUALITY_TOLERANCE) ||
		(fabs(a2.z - a1.z) > EQUALITY_TOLERANCE))
		return false;
	else
		return true;
};

template <class T>
__host__ __device__ inline Point<T> operator-(const Point<T>& rhs)
{
	return Point<T>(-rhs.x, -rhs.y, -rhs.z);
};

template <class T>
__host__ __device__ inline Point<T> operator-(const Point<T>& a1, const Point<T>& a2)
{
	return Point<T>(a1.x - a2.x, a1.y - a2.y, a1.z - a2.z);
};

template <class T>
__host__ __device__ inline Point<T> operator+(const Point<T>& rhs)
{
	return Point<T>(rhs.x, rhs.y, rhs.z);
};

template <class T>
__host__ __device__ inline Point<T> operator+(const Point<T>& a1, const Point<T>& a2)
{
	return Point<T>(a1.x + a2.x, a1.y + a2.y, a1.z + a2.z);
};

template <class T>
__host__ __device__ inline Point<T> operator*(const Point<T>& a1, T k)
{
	return Point<T>(a1.x * k, a1.y * k, a1.z * k);
};

template <class T>
__host__ __device__ inline Point<T> operator*(T k, const Point<T>& a1)
{
	return a1 * k;
};

template <class T>
__host__ __device__ inline Point<T> operator/(const Point<T>& a1, T k)
{
	return Point<T>(a1.x / k, a1.y / k, a1.z / k);
};

template <class T>
__host__ __device__ inline T operator*(const Point<T>& a1, const Point<T>& a2)
{
	return a1.x * a2.x + a1.y * a2.y + a1.z * a2.z;
};

template <class T>
struct WellSegment
{
	Point<T> r1;
	Point<T> r2;
	Point<T> r_bhp;

	T length;
	T pres;
	T pres2D;
	T pres3D;
	T rate;

	__host__ __device__ WellSegment()
	{
	};

	__host__ __device__ WellSegment(const Point<T>& _r1, const Point<T>& _r2, const Point<T> _r_bhp) : r1(_r1), r2(_r2), r_bhp(_r_bhp)
	{
		length = sqrt((r2 - r1) * (r2 - r1));
		pres = pres2D = pres3D = rate = 0.0;
	};

	__host__ __device__ WellSegment& operator=(const WellSegment& rhs) = default;
};

template <class T>
struct Parameters
{
	// Dimensions
	T x_dim, t_dim, p_dim;

	// Spacial params
	Point<T> sizes;
	Point<T> rc;	// middle point of well
	T length;
	T alpha;
	Point<T> r1, r2;

	T rw;

	// Other params
	T visc;
	T kx, kz;
	T perm;
	T rate;

	// Numbers
	int K;
	int M, N, L;
	int I;

	// Grid sizes
	int nx, ny, nz;

	// Integral division limit
	T xi_c;

	// Observation point
	Point<T> r_obs;
};

template <class T>
class Well
{
protected:

	const Point<T> r1;
	const Point<T> r2;
	const int num;
	const T r_w;

	T alpha;
	T length;
	T rate;

public:
	Well(const Point<T>& _r1, const Point<T>& _r2, const int _num, const T _r_w);
	~Well();

	void setRate(T _rate);
	void setUniformRate();

	WellSegment<T>* segs;

	T pres_av;
	T pres_dev;

	void printRates(const Parameters<T>* props);
	void writeRates(const Parameters<T>* props);

	Well& operator=(const Well& well);
};

#endif /* WELL_CUH_ */