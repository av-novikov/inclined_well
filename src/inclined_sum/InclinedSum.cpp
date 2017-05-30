#include <iomanip>
#include <new>

#include "src/inclined_sum/InclinedSum.hpp"

InclinedSum::InclinedSum(const Parameters* _props, const Well* _well) : BaseSum(_props, _well)
{
}
InclinedSum::~InclinedSum()
{	
}
double InclinedSum::get2D(int seg_idx)
{
	double sum = 0.0;
	
	for(int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F2d[seg_idx * props->K + k] * seg.rate / seg.length;
	}
				
	sum *= (props->visc * props->sizes.x / M_PI / M_PI / props->sizes.z / props->kx / sin(props->alpha));
				
	return sum;	
}
double InclinedSum::get3D(int seg_idx)
{
	double sum = 0.0;

	for(int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F3d[seg_idx * props->K + k] * seg.rate / seg.length;
	}
				
	sum *= (2.0 * props->visc / M_PI / props->sizes.x / props->sizes.z / props->kx / cos(props->alpha));
				
	return sum;	
}

void InclinedSum::prepare()
{
	prepare2D();	prepare3D();
}
void InclinedSum::prepare2D()
{
	double buf;
	double sum_prev = 0.0;
	int break_idx = 0;

	for(int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[arr_idx % props->K];
		const Point& r = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;
		
		F2d[arr_idx] = sum_prev = 0.0;
		
		break_idx = 0;
		
		for(int m = 1; m <= props->M; m++)
		{
			buf = sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
					(cos(M_PI * (double)(m) * seg.r1.x / props->sizes.x) -
					cos(M_PI * (double)(m) * seg.r2.x / props->sizes.x));
			
			F2d[arr_idx] += buf / (double)(m) / (double)(m) *
				(exp(-M_PI * (double)(m)* fabs(r.y - props->r1.y) / props->sizes.x) -
					exp(-M_PI * (double)(m)* fabs(r.y + props->r1.y) / props->sizes.x));
			for(int i = 1; i <= props->I; i++)
			{
				F2d[arr_idx] += buf / (double)(m) / (double)(m) * 
					( exp(-M_PI * (double)(m) * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) - 
					exp(-M_PI * (double)(m) * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) );
				F2d[arr_idx] += buf / (double)(m) / (double)(m) *
					(exp(-M_PI * (double)(m)* fabs(r.y - props->r1.y - 2.0 * (double)(i)* props->sizes.y) / props->sizes.x) -
						exp(-M_PI * (double)(m)* fabs(r.y + props->r1.y - 2.0 * (double)(i)* props->sizes.y) / props->sizes.x));
			}
			
			if( fabs(F2d[arr_idx] - sum_prev) > F2d[arr_idx] * EQUALITY_TOLERANCE )
			{
				sum_prev = F2d[arr_idx];
				break_idx = 0;
			} else
				break_idx++;
				
			if(break_idx > 1)
			{
				//std::cout << m << std::endl;
				break;
			}
			
		}
	}
}
void InclinedSum::prepare3D()
{
	double buf1, buf2;
	double sum_prev1 = 0.0, sum_prev2 = 0.0;
	int break_idx1, break_idx2;	
	int k;
	
	for(int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[ arr_idx % props->K ];
		const Point& r = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;
		
		F3d[arr_idx] = sum_prev1 = sum_prev2 = 0.0;
		
		break_idx2 = 0;
		
		for(int m = 1; m <= props->M; m++)
		{
			break_idx1 = 0;
			
			for(int l = 1; l <= props->L; l++)
			{
				buf1 = ((cos(M_PI * (double)(m) * seg.r1.x / props->sizes.x - M_PI * (double)(l) * seg.r1.z / props->sizes.z) -
							cos(M_PI * (double)(m) * seg.r2.x / props->sizes.x - M_PI * (double)(l) * seg.r2.z / props->sizes.z)) /
						( M_PI * (double)(m) * tan(props->alpha) / props->sizes.x + M_PI * (double)(l) / props->sizes.z ) + 
							(cos(M_PI * (double)(m) * seg.r1.x / props->sizes.x + M_PI * (double)(l) * seg.r1.z / props->sizes.z) -
							cos(M_PI * (double)(m) * seg.r2.x / props->sizes.x + M_PI * (double)(l) * seg.r2.z / props->sizes.z)) /
						( M_PI * (double)(m) * tan(props->alpha) / props->sizes.x - M_PI * (double)(l) / props->sizes.z )
							) / 2.0;
							
				buf2 = sqrt((double)(m) * (double)(m) / props->sizes.x / props->sizes.x + (double)(l) * (double)(l) / props->sizes.z / props->sizes.z);
			
				F3d[arr_idx] += buf1 / buf2 *
					(exp(-M_PI * buf2 * fabs(r.y - props->r1.y)) -
						exp(-M_PI * buf2 * fabs(r.y + props->r1.y))) *
					sin(M_PI * (double)(m)* r.x / props->sizes.x) *
					cos(M_PI * (double)(l)* r.z / props->sizes.z);
				for(int i = 1; i <= props->I; i++)
				{	
					F3d[arr_idx] += buf1 / buf2 *
						( exp(-M_PI * buf2 * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y)) - 
						exp(-M_PI * buf2 * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y)) ) *
						sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
						cos(M_PI * (double)(l) * r.z / props->sizes.z);
					F3d[arr_idx] += buf1 / buf2 *
						(exp(-M_PI * buf2 * fabs(r.y - props->r1.y - 2.0 * (double)(i)* props->sizes.y)) -
							exp(-M_PI * buf2 * fabs(r.y + props->r1.y - 2.0 * (double)(i)* props->sizes.y))) *
						sin(M_PI * (double)(m)* r.x / props->sizes.x) *
						cos(M_PI * (double)(l)* r.z / props->sizes.z);
				}
				
				if( fabs(F3d[arr_idx] - sum_prev1) > F3d[arr_idx] * EQUALITY_TOLERANCE )
				{
					sum_prev1 = F3d[arr_idx];
					break_idx1 = 0;
				} else
					break_idx1++;
				
				if(break_idx1 > 1)
				{
					//std::cout << "l=" << l << std::endl;
					break;
				}
			}
			
			if( fabs(F3d[arr_idx] - sum_prev2) > F3d[arr_idx] * EQUALITY_TOLERANCE )
			{
				sum_prev2 = F3d[arr_idx];
				break_idx2 = 0;
			} else
				break_idx2++;
				
			if(break_idx2 > 1)
			{
				std::cout << "m=" << m << std::endl;
				break;
			}
		}
	}
}
double InclinedSum::getPres(const Point& r)
{
	double sum2d = 0.0, sum3d = 0.0;
	double buf1, buf2, buf2d;
	double sum_prev1 = 0.0, sum_prev2 = 0.0, sum_prev_2d = 0.0;
	int break_idx1, break_idx2, break_idx_2d;	

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment seg = well->segs[k];
		break_idx2 = 0;

		for (int m = 1; m <= props->M; m++)
		{
			break_idx1 = 0;

			// 2D-summation
			buf2d = seg.rate / seg.length * sin(M_PI * (double)(m)* r.x / props->sizes.x) *
				(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x) -
					cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x));

			for (int i = -props->I; i <= props->I; i++)
			{
				sum2d += buf2d / (double)(m) / (double)(m) *
					(exp(-M_PI * (double)(m)* fabs(r.y - props->r1.y + 2.0 * (double)(i)* props->sizes.y) / props->sizes.x) -
						exp(-M_PI * (double)(m)* fabs(r.y + props->r1.y + 2.0 * (double)(i)* props->sizes.y) / props->sizes.x));
			}

			for (int l = 1; l <= props->L; l++)
			{
				buf1 = seg.rate / seg.length * ((cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x - M_PI * (double)(l)* seg.r1.z / props->sizes.z) -
					cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x - M_PI * (double)(l)* seg.r2.z / props->sizes.z)) /
					(M_PI * (double)(m)* tan(props->alpha) / props->sizes.x + M_PI * (double)(l) / props->sizes.z) +
					(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x + M_PI * (double)(l)* seg.r1.z / props->sizes.z) -
						cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x + M_PI * (double)(l)* seg.r2.z / props->sizes.z)) /
						(M_PI * (double)(m)* tan(props->alpha) / props->sizes.x - M_PI * (double)(l) / props->sizes.z)
					) / 2.0;

				buf2 = sqrt((double)(m) * (double)(m) / props->sizes.x / props->sizes.x + (double)(l) * (double)(l) / props->sizes.z / props->sizes.z);

				for (int i = -props->I; i <= props->I; i++)
				{
					sum3d += buf1 / buf2 *
						(exp(-M_PI * buf2 * fabs(r.y - props->r1.y + 2.0 * (double)(i)* props->sizes.y)) -
							exp(-M_PI * buf2 * fabs(r.y + props->r1.y + 2.0 * (double)(i)* props->sizes.y))) *
						sin(M_PI * (double)(m)* r.x / props->sizes.x) *
						cos(M_PI * (double)(l)* r.z / props->sizes.z);
				}

				if (fabs(sum3d - sum_prev1) > sum3d * EQUALITY_TOLERANCE)
				{
					sum_prev1 = sum3d;
					break_idx1 = 0;
				}
				else
					break_idx1++;

				if (break_idx1 > 1)
				{
					//std::cout << "l=" << l << std::endl;
					break;
				}
			}

			if (fabs(sum3d - sum_prev2) > sum3d * EQUALITY_TOLERANCE)
			{
				sum_prev2 = sum3d;
				break_idx2 = 0;
			}
			else
				break_idx2++;

			if (fabs(sum2d - sum_prev_2d) > sum2d * EQUALITY_TOLERANCE)
			{
				sum_prev_2d = sum2d;
				break_idx_2d = 0;
			}
			else
				break_idx_2d++;

			if (break_idx2 > 1 && break_idx_2d > 1)
			{
				std::cout << "m=" << m << std::endl;
				break;
			}
		}
	}

	sum2d *= (props->visc * props->sizes.x / M_PI / M_PI / props->sizes.z / props->kx / sin(props->alpha));
	sum3d *= (2.0 * props->visc / M_PI / props->sizes.x / props->sizes.z / props->kx / cos(props->alpha));

	return sum2d + sum3d;
}