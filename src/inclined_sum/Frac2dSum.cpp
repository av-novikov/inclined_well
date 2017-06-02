#include "src/inclined_sum/Frac2dSum.hpp"
#include "boost/math/special_functions/expint.hpp"

#include <valarray>

using boost::math::expint;

Frac2dSum::Frac2dSum(const Parameters* _props, const Well* _well) : BaseSum(_props, _well)
{
}
Frac2dSum::~Frac2dSum()
{
}
double Frac2dSum::getPres(const Point& point)
{
	return 0.0;
}
double Frac2dSum::get2D(int seg_idx)
{
	double sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F2d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (4.0 * props->visc / props->sizes.x / props->sizes.y / props->sizes.z / props->kx / cos(props->alpha));
	return sum;
}
double Frac2dSum::get3D(int seg_idx)
{
	double sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F3d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (props->visc / 4.0 / M_PI / props->sizes.z / props->kx / cos(props->alpha));
	return sum;
}
void Frac2dSum::prepare()
{
	prepareDirect();
	prepareFourier();
}
void Frac2dSum::prepareDirect()
{
	double sum = 0.0;
	double sum_prev;
	double buf, F;
	int break_idx = 0;

	for (int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[arr_idx % props->K];
		const Point& r = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;

		F2d[arr_idx] = sum_prev = 0.0;

		break_idx = 0;
		for (int m = 1; m <= props->M; m++)
		{
			for (int n = 1; n <= props->M; n++)
			{
				F = ((sin(M_PI * (double)m * seg.r2.x / props->sizes.x - M_PI * (double)m * seg.r2.y / props->sizes.y) -
					sin(M_PI * (double)m * seg.r1.x / props->sizes.x - M_PI * (double)m * seg.r1.y / props->sizes.y)) /
					(M_PI * (double)m / props->sizes.x - M_PI * (double)n / props->sizes.y * tan(props->alpha)) -
					(sin(M_PI * (double)m * seg.r2.x / props->sizes.x + M_PI * (double)m * seg.r2.y / props->sizes.y) -
						sin(M_PI * (double)m * seg.r1.x / props->sizes.x + M_PI * (double)m * seg.r1.y / props->sizes.y)) /
						(M_PI * (double)m / props->sizes.x + M_PI * (double)n / props->sizes.y * tan(props->alpha))) / 2.0;
				buf = M_PI * M_PI * ((double)m * (double)m / props->sizes.x / props->sizes.x + 
										(double)n * (double)n / props->sizes.y / props->sizes.y);
				F2d[arr_idx] += F * sin(M_PI * (double)m * r.x / props->sizes.x) * 
									sin(M_PI * (double)n * r.y / props->sizes.y) *
									exp(-props->xi_c * buf) / buf;
			}

			if (fabs(F2d[arr_idx] - sum_prev) > F2d[arr_idx] * EQUALITY_TOLERANCE)
			{
				sum_prev = F2d[arr_idx];
				break_idx = 0;
			}
			else
				break_idx++;

			if (break_idx > 1)
			{
				//std::cout << m << std::endl;
				break;
			}
		}
	}
}
void Frac2dSum::prepareFourier()
{
	double buf1[3], buf2[3], buf3[3], buf4[3], res;
	Point v1, v2, v3, v4, v5, v6;
	int break_idx = 0;

	for (int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[arr_idx % props->K];
		const Point& r = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;

		F3d[arr_idx] = 0.0;

		for (int p = -props->I; p <= props->I; p++)
		{
			for (int q = -props->I; q <= props->I; q++)
			{
				Point pt((double)p, (double)q, 0.0);
				const Point av_r = (seg.r1 + seg.r2) / 2.0;
				v1 = (r - seg.r1 + 2.0 * product(pt, props->sizes));	v1 = product(v1, v1) / 4.0 / props->xi_c;
				v2 = (r + seg.r1 + 2.0 * product(pt, props->sizes));	v2 = product(v2, v2) / 4.0 / props->xi_c;
				v3 = (r - seg.r2 + 2.0 * product(pt, props->sizes));	v3 = product(v3, v3) / 4.0 / props->xi_c;
				v4 = (r + seg.r2 + 2.0 * product(pt, props->sizes));	v4 = product(v4, v4) / 4.0 / props->xi_c;
				v5 = (r - av_r + 2.0 * product(pt, props->sizes));		v5 = product(v5, v5) / 4.0 / props->xi_c;
				v6 = (r + av_r + 2.0 * product(pt, props->sizes));		v6 = product(v6, v6) / 4.0 / props->xi_c;

				buf1[0] = v1.x + v1.y;		buf2[0] = v1.x + v2.y;
				buf3[0] = v2.x + v1.y;		buf4[0] = v2.x + v2.y;
				buf1[1] = v3.x + v3.y;		buf2[1] = v3.x + v4.y;
				buf3[1] = v4.x + v3.y;		buf4[1] = v4.x + v4.y;
				buf1[2] = v5.x + v5.y;		buf2[2] = v5.x + v6.y;
				buf3[2] = v6.x + v5.y;		buf4[2] = v6.x + v6.y;

				res = (expint(-buf1[0]) + expint(-buf1[1]) + 4.0 * expint(-buf1[2]));
				res -= (expint(-buf2[0]) + expint(-buf2[1]) + 4.0 * expint(-buf2[2]));
				res -= (expint(-buf3[0]) + expint(-buf3[1]) + 4.0 * expint(-buf3[2]));
				res += (expint(-buf4[0]) + expint(-buf4[1]) + 4.0 * expint(-buf4[2]));

				F3d[arr_idx] -= res * (seg.r2.x - seg.r1.x) / 6.0;
			}
		}
	}
}
