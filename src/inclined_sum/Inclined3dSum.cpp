#include <new>

#include "src/inclined_sum/Inclined3dSum.h"

Inclined3dSum::Inclined3dSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well) : BaseSum(_sprops, _props, _well)
{
}
Inclined3dSum::~Inclined3dSum()
{
}
double Inclined3dSum::getPres(const Point& point)
{
	return 0.0;
}
double Inclined3dSum::get2D(int seg_idx)
{
	double sum = 0.0;
	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F2d[seg_idx * sprops.K + k] * seg.rate / seg.length;
	}

	sum *= (8.0 * props->visc * gprops->length / props->sizes.x / props->sizes.y / props->sizes.z / props->kx);
	return sum;
}
double Inclined3dSum::get3D(int seg_idx)
{
	double sum = 0.0;
	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F3d[seg_idx * sprops.K + k] * seg.rate / seg.length;
	}

	sum *= (props->visc * gprops->length / 8.0 / M_PI / props->kx);
	return sum;
}
void Inclined3dSum::prepare()
{
	size = segs->size() * sprops.K;
	F2d = new double[size];		F3d = new double[size];
	prepareDirect();	prepareFourier();
}
void Inclined3dSum::prepareDirect()
{
	double sum = 0.0;
	double sum_prev_m, sum_prev_n;
	double buf, F1, F2;
	int break_idx_m = 0, break_idx_n = 0;

	for (int arr_idx = 0; arr_idx < size; arr_idx++)
	{
		const WellSegment& seg = well->segs[arr_idx % sprops.K];
		const Point& r = (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->r_bhp;
		const Point rad = gprops->r2 - gprops->r1;

		F2d[arr_idx] = sum_prev_m = 0.0;

		break_idx_m = 0;
		for (int m = 1; m <= sprops.M; m++)
		{
			sum_prev_n = 0;	 break_idx_n = 0;
			for (int n = 1; n <= sprops.M; n++)
			{
				F1 = ((sin(M_PI * (double)m * seg.r2.x / props->sizes.x - M_PI * (double)n * seg.r2.y / props->sizes.y) -
						sin(M_PI * (double)m * seg.r1.x / props->sizes.x - M_PI * (double)n * seg.r1.y / props->sizes.y)) /
						(M_PI * (double)m / props->sizes.x * rad.x - M_PI * (double)n / props->sizes.y * rad.y) -
					(sin(M_PI * (double)m * seg.r2.x / props->sizes.x + M_PI * (double)n * seg.r2.y / props->sizes.y) -
						sin(M_PI * (double)m * seg.r1.x / props->sizes.x + M_PI * (double)n * seg.r1.y / props->sizes.y)) /
						(M_PI * (double)m / props->sizes.x * rad.x + M_PI * (double)n / props->sizes.y * rad.y)) / 2.0;
				buf = M_PI * M_PI * ((double)m * (double)m / props->sizes.x / props->sizes.x + (double)n * (double)n / props->sizes.y / props->sizes.y);
				F2d[arr_idx] += 1.0 / 2.0 * F1 * sin(M_PI * (double)m * r.x / props->sizes.x) *	sin(M_PI * (double)n * r.y / props->sizes.y) * exp(-sprops.xi_c * buf) / buf;

				for (int l = 1; l <= sprops.L; l++)
				{
					auto getFoo = [=, this](const Point& pt) -> double
					{
						const Point point = M_PI * pt / props->sizes;
						const Point p1((double)m, (double)n, -(double)l);
						const Point p2((double)m, -(double)n, (double)l);
						const Point p3((double)m, -(double)n, -(double)l);
						const Point p4((double)m, (double)n, (double)l);
						return 1.0 / 4.0 * (
							-sin(p1 * point) / (M_PI * p1 * (rad / props->sizes)) +
							sin(p2 * point) / (M_PI * p2 * (rad / props->sizes)) +
							sin(p3 * point) / (M_PI * p3 * (rad / props->sizes)) -
							sin(p4 * point) / (M_PI * p4 * (rad / props->sizes)));
					};

					F2 = getFoo(seg.r2) - getFoo(seg.r1);
					buf = M_PI * M_PI * ((double)m * (double)m / props->sizes.x / props->sizes.x +
						(double)n * (double)n / props->sizes.y / props->sizes.y + (double)l * (double)l / props->sizes.z / props->sizes.z);
					F2d[arr_idx] += F2 * sin(M_PI * (double)m * r.x / props->sizes.x) *	sin(M_PI * (double)n * r.y / props->sizes.y) *
						cos(M_PI * (double)l * r.z / props->sizes.z) * exp(-sprops.xi_c * buf) / buf;
				}

				if (fabs(F2d[arr_idx] - sum_prev_n) > F2d[arr_idx] * EQUALITY_TOLERANCE)
				{
					sum_prev_n = F2d[arr_idx];
					break_idx_n = 0;
				}
				else
					break_idx_n++;

				if (break_idx_n > 1)
				{
					//std::cout << m << std::endl;
					//break;
				}
			}

			if (fabs(F2d[arr_idx] - sum_prev_m) > F2d[arr_idx] * EQUALITY_TOLERANCE)
			{
				sum_prev_m = F2d[arr_idx];
				break_idx_m = 0;
			}
			else
				break_idx_m++;

			if (break_idx_m > 1)
			{
				//std::cout << m << std::endl;
				//break;
			}
		}
	}
}
void Inclined3dSum::prepareFourier()
{
	for (int arr_idx = 0; arr_idx < size; arr_idx++)
	{
		const WellSegment& seg = well->segs[arr_idx % sprops.K];
		const Point& r = (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->r_bhp;

		F3d[arr_idx] = 0.0;

		for (int p = -sprops.I; p <= sprops.I; p++)
		{
			for (int q = -sprops.I; q <= sprops.I; q++)
			{
				for (int s = -100 * sprops.I; s <= 100 * sprops.I; s++)
				{
					Point pt((double)p, (double)q, (double)s);
					auto getFoo = [=, this](const Point& point) -> double
					{
						Point vv1, vv2;
						vv1 = (r - point + 2.0 * product(pt, props->sizes));	vv1 = product(vv1, vv1) / 4.0;
						vv2 = (r + point + 2.0 * product(pt, props->sizes));	vv2 = product(vv2, vv2) / 4.0;
						double bbuf1, bbuf2, bbuf3, bbuf4, bbuf5, bbuf6, bbuf7, bbuf8;
						bbuf1 = vv1.x + vv1.y + vv1.z;		bbuf2 = vv1.x + vv1.y + vv2.z;
						bbuf3 = vv1.x + vv2.y + vv1.z;		bbuf4 = vv1.x + vv2.y + vv2.z;
						bbuf5 = vv2.x + vv1.y + vv1.z;		bbuf6 = vv2.x + vv1.y + vv2.z;
						bbuf7 = vv2.x + vv2.y + vv1.z;		bbuf8 = vv2.x + vv2.y + vv2.z;
						return	erfc(sqrt(bbuf1 / sprops.xi_c)) / sqrt(bbuf1)
							+ erfc(sqrt(bbuf2 / sprops.xi_c)) / sqrt(bbuf2)
							- erfc(sqrt(bbuf3 / sprops.xi_c)) / sqrt(bbuf3)
							- erfc(sqrt(bbuf4 / sprops.xi_c)) / sqrt(bbuf4)
							- erfc(sqrt(bbuf5 / sprops.xi_c)) / sqrt(bbuf5)
							- erfc(sqrt(bbuf6 / sprops.xi_c)) / sqrt(bbuf6)
							+ erfc(sqrt(bbuf7 / sprops.xi_c)) / sqrt(bbuf7)
							+ erfc(sqrt(bbuf8 / sprops.xi_c)) / sqrt(bbuf8);
					};

					Double2d* integrand = new Double2d[PART_SIZE + 1];
					for (int i = 0; i < PART_SIZE + 1; i++)
					{
						Point point = seg.r1 + (double)i * (seg.r2 - seg.r1) / (double)PART_SIZE;
						double tau = seg.tau1 + (double)i * (seg.tau2 - seg.tau1) / (double)PART_SIZE;
						integrand[i].x = tau;
						double qwe = getFoo(point);
						integrand[i].y = getFoo(point);
					}
					integr = new Integral(integrand, PART_SIZE + 1);
					F3d[arr_idx] += integr->Calculate(seg.tau1, seg.tau2);

					delete integrand;
					delete integr;
				}
			}
		}
	}
}