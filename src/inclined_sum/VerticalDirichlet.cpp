#include <iostream>

#include "boost/math/special_functions/expint.hpp"
#include "src/inclined_sum/VerticalDirichlet.h"

using boost::math::expint;

VerticalDirichlet::VerticalDirichlet(const Parameters* _props, const Well* _well) : BaseSum(_props, _well)
{
}

VerticalDirichlet::~VerticalDirichlet()
{
}

void VerticalDirichlet::prepare()
{
	F2d[0] = directSum() + fourierSum();
}

double VerticalDirichlet::get2D(int seg_idx)
{
	return F2d[0];
}

double VerticalDirichlet::getPres(const Point& p)
{
	return 0.0;
}

double VerticalDirichlet::directSum()
{
	double sum = 0.0;
	double buf1, buf2;
	const Point& r = well->segs[0].r_bhp;

	for (int m = 1; m < props->M; m++)
	{
		buf1 = sin(M_PI * (double)(m)* props->rc.x / props->sizes.x) * sin(M_PI * (double)(m)* props->rc.x / props->sizes.x);
		for (int n = 1; n < props->N; n++)
		{
			buf2 = M_PI * M_PI * ((double)(m * m) / props->sizes.x / props->sizes.x + (double)(n * n) / props->sizes.y / props->sizes.y);

			sum += exp(-buf2 * props->xi_c) * buf1 / buf2 *
				sin(M_PI * (double)(n)* r.y / props->sizes.y) * sin(M_PI * (double)(n)* props->rc.y / props->sizes.y);
		}
	}

	return 4.0 * props->visc * props->rate / props->sizes.x / props->sizes.y / props->sizes.z / props->kx * sum;
}

double VerticalDirichlet::fourierSum()
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

			sum += expint(1, (buf11 * buf11 + buf21 * buf21) / props->xi_c) - expint(1, (buf11 * buf11 + buf22 * buf22) / props->xi_c) -
				expint(1, (buf12 * buf12 + buf21 * buf21) / props->xi_c) + expint(1, (buf12 * buf12 + buf22 * buf22) / props->xi_c);
		}
	}

	return props->visc * props->rate / props->sizes.z / props->kx / 4.0 / M_PI * sum;
}

double VerticalDirichlet::getAnalyticalPres()
{
	const Point& r = well->segs[0].r_bhp;
	double qwe = r.y - props->rc.y;
	return -props->visc * props->rate / props->sizes.z / props->kx / 4.0 / M_PI *
			(log((1.0 - 2.0 * exp(-M_PI / props->sizes.x * (r.y - props->rc.y)) * cos(M_PI * (r.x - props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI / props->sizes.x * (r.y - props->rc.y))) /
			((1.0 - 2.0 * exp(-M_PI / props->sizes.x * (r.y - props->rc.y)) * cos(M_PI * (r.x + props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI / props->sizes.x * (r.y - props->rc.y))))) +
			log((1.0 - 2.0 * exp(-M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y - props->rc.y) / props->sizes.y)) * cos(M_PI * (r.x - props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y - props->rc.y) / props->sizes.y))) /
			((1.0 - 2.0 * exp(-M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y - props->rc.y) / props->sizes.y)) * cos(M_PI * (r.x + props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y - props->rc.y) / props->sizes.y))))) -
			log((1.0 - 2.0 * exp(-M_PI / props->sizes.x * (r.y + props->rc.y)) * cos(M_PI * (r.x - props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI / props->sizes.x * (r.y + props->rc.y))) /
			((1.0 - 2.0 * exp(-M_PI / props->sizes.x * (r.y + props->rc.y)) * cos(M_PI * (r.x + props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI / props->sizes.x * (r.y + props->rc.y))))) -
			log((1.0 - 2.0 * exp(-M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y + props->rc.y) / props->sizes.y)) * cos(M_PI * (r.x - props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y + props->rc.y) / props->sizes.y))) /
			((1.0 - 2.0 * exp(-M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y + props->rc.y) / props->sizes.y)) * cos(M_PI * (r.x + props->rc.x) / props->sizes.x) +
				exp(-2.0 * M_PI * props->sizes.y / props->sizes.x * (2.0 - (r.y + props->rc.y) / props->sizes.y))))));
}