#include "src/inclined_sum/Frac2dSum.hpp"
#include "boost/math/special_functions/expint.hpp"

#include <valarray>

using boost::math::expint;

Frac2dSum::Frac2dSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well) : BaseSum(_sprops, _props, _well)
{
}
Frac2dSum::~Frac2dSum()
{
}
double Frac2dSum::get2D(int seg_idx)
{
	double sum = 0.0;
	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F2d[seg_idx * sprops.K + k] * seg.rate / seg.length;
	}

	sum *= (4.0 * props->visc / props->sizes.x / props->sizes.y / props->sizes.z / props->kx / cos(gprops->alpha));
	return sum;
}
double Frac2dSum::get3D(int seg_idx)
{
	double sum = 0.0;

	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F3d[seg_idx * sprops.K + k] * seg.rate / seg.length;
	}

	sum *= (props->visc / 4.0 / M_PI / props->sizes.z / props->kx / cos(gprops->alpha));
	return sum;
}
void Frac2dSum::prepare()
{
	size = segs->size() * sprops.K;
	F2d = new double[size];		F3d = new double[size];
	prepareDirect();
	prepareFourier();
}
void Frac2dSum::prepareDirect()
{
	double sum_prev;
	double buf, F;
	int break_idx = 0;

	for (int arr_idx = 0; arr_idx < size; arr_idx++)
	{
		const WellSegment& seg = well->segs[arr_idx % sprops.K];
		const Point& r = (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->r_bhp;

		F2d[arr_idx] = sum_prev = 0.0;

		break_idx = 0;
		for (int m = 1; m <= sprops.M; m++)
		{
			for (int n = 1; n <= sprops.M; n++)
			{
				F = ((sin(M_PI * (double)m * seg.r2.x / props->sizes.x - M_PI * (double)n * seg.r2.y / props->sizes.y) -
					sin(M_PI * (double)m * seg.r1.x / props->sizes.x - M_PI * (double)n * seg.r1.y / props->sizes.y)) /
					(M_PI * (double)m / props->sizes.x - M_PI * (double)n / props->sizes.y * tan(gprops->alpha)) -
					(sin(M_PI * (double)m * seg.r2.x / props->sizes.x + M_PI * (double)n * seg.r2.y / props->sizes.y) -
						sin(M_PI * (double)m * seg.r1.x / props->sizes.x + M_PI * (double)n * seg.r1.y / props->sizes.y)) /
						(M_PI * (double)m / props->sizes.x + M_PI * (double)n / props->sizes.y * tan(gprops->alpha))) / 2.0;
				buf = M_PI * M_PI * ((double)m * (double)m / props->sizes.x / props->sizes.x + 
										(double)n * (double)n / props->sizes.y / props->sizes.y);
				F2d[arr_idx] += F * sin(M_PI * (double)m * r.x / props->sizes.x) * 
									sin(M_PI * (double)n * r.y / props->sizes.y) *
									exp(-sprops.xi_c * buf) / buf;
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
				//break;
			}
		}
	}
}
void Frac2dSum::prepareFourier()
{
	int break_idx = 0;

	for (int arr_idx = 0; arr_idx < size; arr_idx++)
	{
		const WellSegment& seg = well->segs[arr_idx % sprops.K];
		const Point& r = (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->r_bhp;

		F3d[arr_idx] = 0.0;

		for (int p = -sprops.I; p <= sprops.I; p++)
		{
			for (int q = -sprops.I; q <= sprops.I; q++)
			{
				Point pt((double)p, (double)q, 0.0);
				auto getFoo = [=, this](const Point& point) -> double
				{
					Point vv1, vv2;
					vv1 = (r - point + 2.0 * product(pt, props->sizes));	vv1 = product(vv1, vv1) / 4.0 / sprops.xi_c;
					vv2 = (r + point + 2.0 * product(pt, props->sizes));	vv2 = product(vv2, vv2) / 4.0 / sprops.xi_c;
					double bbuf1, bbuf2, bbuf3, bbuf4;
					bbuf1 = vv1.x + vv1.y;		bbuf2 = vv1.x + vv2.y;
					bbuf3 = vv2.x + vv1.y;		bbuf4 = vv2.x + vv2.y;
					return expint(-bbuf1) - expint(-bbuf2) - expint(-bbuf3) + expint(-bbuf4);
				};

				Double2d* integrand = new Double2d[PART_SIZE + 1];
				for (int i = 0; i < PART_SIZE + 1; i++)
				{
					Point point = seg.r1 + (double)i * (seg.r2 - seg.r1) / (double)PART_SIZE;
					integrand[i].x = point.x;
					integrand[i].y = getFoo(point);
				}
				integr = new Integral(integrand, PART_SIZE + 1);
				F3d[arr_idx] -= integr->Calculate(seg.r1.x, seg.r2.x);

				delete integrand;
				delete integr;
			}
		}
	}
}

double Frac2dSum::getPressure(const Point& r)
{
	double sum1 = 0.0, sum2 = 0.0;
	double sum_prev = 0.0;
	double buf, F;
	int break_idx = 0;

	for (size_t seg_idx = 0; seg_idx < sprops.K; seg_idx++)
	{
		const WellSegment& seg = well->segs[seg_idx];

		break_idx = 0;
		for (int m = 1; m <= sprops.M; m++)
		{
			for (int n = 1; n <= sprops.M; n++)
			{
				F = ((sin(M_PI * (double)m * seg.r2.x / props->sizes.x - M_PI * (double)n * seg.r2.y / props->sizes.y) -
					sin(M_PI * (double)m * seg.r1.x / props->sizes.x - M_PI * (double)n * seg.r1.y / props->sizes.y)) /
					(M_PI * (double)m / props->sizes.x - M_PI * (double)n / props->sizes.y * tan(gprops->alpha)) -
					(sin(M_PI * (double)m * seg.r2.x / props->sizes.x + M_PI * (double)n * seg.r2.y / props->sizes.y) -
						sin(M_PI * (double)m * seg.r1.x / props->sizes.x + M_PI * (double)n * seg.r1.y / props->sizes.y)) /
						(M_PI * (double)m / props->sizes.x + M_PI * (double)n / props->sizes.y * tan(gprops->alpha))) / 2.0;
				buf = M_PI * M_PI * ((double)m * (double)m / props->sizes.x / props->sizes.x +
					(double)n * (double)n / props->sizes.y / props->sizes.y);
				sum1 += seg.rate / seg.length * F * sin(M_PI * (double)m * r.x / props->sizes.x) *
					sin(M_PI * (double)n * r.y / props->sizes.y) *
					exp(-sprops.xi_c * buf) / buf;
			}

			if (fabs(sum1 - sum_prev) > sum1 * EQUALITY_TOLERANCE)
			{
				sum_prev = sum1;
				break_idx = 0;
			}
			else
				break_idx++;

			if (break_idx > 1)
			{
				//std::cout << m << std::endl;
				//break;
			}
		}
	}
	sum1 *= (4.0 * props->visc / props->sizes.x / props->sizes.y / props->sizes.z / props->kx / cos(gprops->alpha));

	break_idx = 0;
	for (size_t seg_idx = 0; seg_idx < sprops.K; seg_idx++)
	{
		const WellSegment& seg = well->segs[seg_idx];

		for (int p = -sprops.I; p <= sprops.I; p++)
		{
			for (int q = -sprops.I; q <= sprops.I; q++)
			{
				Point pt((double)p, (double)q, 0.0);
				auto getFoo = [=, this](const Point& point) -> double
				{
					Point vv1, vv2;
					vv1 = (r - point + 2.0 * product(pt, props->sizes));	vv1 = product(vv1, vv1) / 4.0 / sprops.xi_c;
					vv2 = (r + point + 2.0 * product(pt, props->sizes));	vv2 = product(vv2, vv2) / 4.0 / sprops.xi_c;
					double bbuf1, bbuf2, bbuf3, bbuf4;
					bbuf1 = vv1.x + vv1.y;		bbuf2 = vv1.x + vv2.y;
					bbuf3 = vv2.x + vv1.y;		bbuf4 = vv2.x + vv2.y;
					return expint(-bbuf1) - expint(-bbuf2) - expint(-bbuf3) + expint(-bbuf4);
				};

				Double2d* integrand = new Double2d[PART_SIZE + 1];
				for (int i = 0; i < PART_SIZE + 1; i++)
				{
					Point point = seg.r1 + (double)i * (seg.r2 - seg.r1) / (double)PART_SIZE;
					integrand[i].x = point.x;
					integrand[i].y = getFoo(point);
				}
				integr = new Integral(integrand, PART_SIZE + 1);
				sum2 -= seg.rate / seg.length * integr->Calculate(seg.r1.x, seg.r2.x);

				delete integrand;
				delete integr;
			}
		}
	}
	sum2 *= (props->visc / 4.0 / M_PI / props->sizes.z / props->kx / cos(gprops->alpha));

	return sum1 + sum2;
}