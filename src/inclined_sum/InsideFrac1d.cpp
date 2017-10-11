#include "src/inclined_sum/InsideFrac1D.hpp"

InsideFrac1d::InsideFrac1d(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well) : BaseSum(_sprops, _props, _well)
{
}
InsideFrac1d::~InsideFrac1d()
{
}
double InsideFrac1d::get2D(int seg_idx)
{
	double sum = 0.0;
	for (int k = mid_idx; k < mid_idx + seg_idx; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += seg.rate;
	}

	return props->visc / well->getGeomProps()->rw / props->kf / props->sizes.z * (props->rate / 2.0 - sum);
	
	/*const double length = well->getGeomProps()->length;
	double sum = 0.0;
	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum -= F2d[seg_idx * sprops.K + k] * seg.rate / seg.length;
	}

	return 2.0 * props->visc / 10000.0 / props->kx / props->sizes.z / 2.0 / well->getGeomProps()->rw / length * sum + 
			p1 + (p2 - p1) * ((*segs)[seg_idx]->tau1 + (*segs)[seg_idx]->tau2) / 2.0;*/
}
double InsideFrac1d::get3D(int seg_idx)
{
	return 0.0;
}
void InsideFrac1d::setBounds(double _p1, double _p2)
{
	//p1 = _p1;		p2 = _p2;
}
void InsideFrac1d::prepare()
{
	mid_idx = sprops.K / 2;
	/*size = segs->size() * sprops.K;
	F2d = new double[size];		F3d = new double[size];
	
	double sum = 0.0;
	double sum_prev;
	double buf, F;
	int break_idx = 0;

	for (int arr_idx = 0; arr_idx < size; arr_idx++)
	{
		const WellSegment& seg = well->segs[arr_idx % sprops.K];
		const double length = well->getGeomProps()->length;
		const double x1 = length * (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->tau1;
		const double x2 = length * (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->tau2;
		const double x = (x1 + x2) / 2.0;

		F2d[arr_idx] = sum_prev = 0.0;
		break_idx = 0;
		for (int m = 1; m < sprops.M; m++)
		{
			F2d[arr_idx] += sin(M_PI * (double)m * x / length) *
				((cos(M_PI * (double)m * x1 / length) - cos(M_PI * (double)m * x2 / length)) / 
				(M_PI * (double)m / length) - props->rate * seg.length / seg.rate * sin(M_PI * (double)m / 2.0)) /
				(M_PI * M_PI * (double)m * (double)m / length / length);

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
	}*/
}
