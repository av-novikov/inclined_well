#include "boost/math/special_functions/expint.hpp"
#include "src\inclined_sum\VerticalNeumann.h"

using boost::math::expint;

VerticalNeumann::VerticalNeumann(const Parameters* _props, const Well* _well) : BaseSum(_props, _well)
{
}

VerticalNeumann::~VerticalNeumann()
{
}

void VerticalNeumann::prepare()
{
	F2d[0] = directSum() + fourierSum();
}

double VerticalNeumann::get2D(int seg_idx)
{
	return F2d[0];
}

double VerticalNeumann::directSum()
{
	double sum = 0.0, sum_m = 0.0, sum_n = 0.0;
	double buf1, buf2;
	const Point& r = well->segs[0].r_bhp;

	for (int m = 1; m < props->M; m++)
	{
		// 1-D N-series summation
		buf1 = cos(M_PI * (double)(m)* r.y / props->sizes.y) * cos(M_PI * (double)(m)* props->rc.y / props->sizes.y);
		buf2 = M_PI * M_PI * (double)(m * m) / props->sizes.y / props->sizes.y;
		sum_n += buf1 / buf2 * exp(-buf2 * props->xi_c);

		// 1-D M-series summation
		buf1 = cos(M_PI * (double)(m)* props->rc.x / props->sizes.x);
		buf1 *= buf1;
		buf2 = M_PI * M_PI * (double)(m * m) / props->sizes.x / props->sizes.x;
		sum_m += buf1 / buf2 * exp(-buf2 * props->xi_c);

		// 2-D M-N-series summation
		for (int n = 1; n < props->M; n++)
		{
			buf2 = M_PI * M_PI * ((double)(m * m) / props->sizes.x / props->sizes.x + (double)(n * n) / props->sizes.y / props->sizes.y);

			sum += exp(-buf2 * props->xi_c) * buf1 / buf2 *
				cos(M_PI * (double)(n)* r.y / props->sizes.y) * cos(M_PI * (double)(n)* props->rc.y / props->sizes.y);
		}


	}

	sum *= (4.0 * props->visc * props->rate / props->sizes.x / props->sizes.y / props->sizes.z / props->kx);
	sum_m *= (2.0 * props->visc * props->rate / props->sizes.x / props->sizes.y / props->sizes.z / props->kx);
	sum_n *= (2.0 * props->visc * props->rate / props->sizes.x / props->sizes.y / props->sizes.z / props->kx);

	return sum + sum_m + sum_n;
}

double VerticalNeumann::fourierSum()
{
	double sum = 0.0;
	double buf11, buf12, buf21, buf22;
	const Point& r = well->segs[0].r_bhp;

	for (int p = -props->I; p < props->I; p++)
	{
		for (int q = -props->I; q < props->I; q++)
		{
			buf11 = (r.x - props->rc.x + 2.0 * p * props->sizes.x) / 2.0;
			buf12 = (r.x + props->rc.x + 2.0 * p * props->sizes.x) / 2.0;
			buf21 = (r.y - props->rc.y + 2.0 * q * props->sizes.y) / 2.0;
			buf22 = (r.y + props->rc.y + 2.0 * q * props->sizes.y) / 2.0;

			sum += expint(1, (buf11 * buf11 + buf21 * buf21) / props->xi_c) + expint(1, (buf11 * buf11 + buf22 * buf22) / props->xi_c) +
				expint(1, (buf12 * buf12 + buf21 * buf21) / props->xi_c) + expint(1, (buf12 * buf12 + buf22 * buf22) / props->xi_c);
		}
	}
	
	return props->visc * props->rate / props->kx * (sum / props->sizes.z / 4.0 / M_PI - props->xi_c / props->sizes.x / props->sizes.y / props->sizes.z);
}