#include "src/inclined_sum/Frac2dSum.hpp"

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

	sum *= (4.0 * props->visc * props->rate / props->sizes.x / props->sizes.y / props->sizes.z / props->kx / cos(props->alpha));
	return sum;
}
void Frac2dSum::prepare()
{
	double sum = 0.0;
	double sum_prev;
	double buf;
	int break_idx = 0;

	for (int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[arr_idx % props->K];
		const Point& r = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;

		F2d[arr_idx] = sum_prev = 0.0;

		auto getInt = [=, this](const int m, const int n) -> double
		{
			return		((sin(M_PI * (double)m * seg.r2.x / props->sizes.x - M_PI * (double)m * seg.r2.y / props->sizes.y) -
				sin(M_PI * (double)m * seg.r1.x / props->sizes.x - M_PI * (double)m * seg.r1.y / props->sizes.y)) /
				(M_PI * (double)m / props->sizes.x - M_PI * (double)n / props->sizes.y * tan(props->alpha)) -
				(sin(M_PI * (double)m * seg.r2.x / props->sizes.x + M_PI * (double)m * seg.r2.y / props->sizes.y) -
					sin(M_PI * (double)m * seg.r1.x / props->sizes.x + M_PI * (double)m * seg.r1.y / props->sizes.y)) /
					(M_PI * (double)m / props->sizes.x + M_PI * (double)n / props->sizes.y * tan(props->alpha))) / 2.0;
		};

		break_idx = 0;
		for (int m = 1; m <= props->M; m++)
		{
			buf = sin(M_PI * (double)m * r.x / props->sizes.x) *
					exp(-props->xi_c * M_PI * M_PI * (double)m * (double)m / props->sizes.x / props->sizes.x);
			for (int n = 1; n <= props->M; n++)
			{
				double qwe = getInt(m, n);
				F2d[arr_idx] += buf * getInt(m, n) * sin(M_PI * (double)n * r.y / props->sizes.y) *
					exp(-props->xi_c * M_PI * M_PI * (double)n * (double)n / props->sizes.y / props->sizes.y) / 
					M_PI / M_PI / ((double)m * (double)m / props->sizes.x / props->sizes.x + (double)n * (double)n / props->sizes.y / props->sizes.y);
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
double Frac2dSum::fourierSum(const Point& r)
{
	/*double sum = 0.0;
	double buf11, buf12, buf21, buf22;

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
	*/
	return 0.0;
}