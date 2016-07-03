#include <iomanip>
#include <new>

#include "src/inclined_sum/InclinedSum.hpp"

InclinedSum::InclinedSum(const Parameters* _props, const Well* _well) : props(_props), well(_well)
{
	F = new double [props->K * props->K];

	prepare3D();
}

InclinedSum::~InclinedSum()
{	
	delete [] F;
}

double InclinedSum::get2D(int seg_idx)
{
	double sum = 0.0;
	double sum_prev = 0.0;
	int break_idx = 0;
	const Point& r = well->segs[seg_idx].r_bhp;

	for(int k = 0; k < props->K; k++)
	{
		break_idx = 0;
		
		for(int m = 1; m <= props->M; m++)
		{
			for(int i = -props->I; i <= props->I; i++)
			{
				sum += well->segs[k].rate / well->segs[k].length / (double)(m) / (double)(m) * 
					( exp(-M_PI * (double)(m) * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) - 
					exp(-M_PI * (double)(m) * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) ) *
					sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props->alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props->sizes.x));
			}
			
			if( fabs(sum - sum_prev) > sum * EQUALITY_TOLERANCE )
			{
				sum_prev = sum;
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
	
	sum *= (props->visc * props->sizes.x / M_PI / M_PI / props->sizes.z / props->kx / sin(props->alpha));

	return sum;	
}

double InclinedSum::get3D(int seg_idx)
{
	double sum = 0.0;
	
	int k;
	
	for(k = 0; k < props->K; k++)
	{
		sum += F[seg_idx * props->K + k] * well->segs[k].rate;
	}
				
	sum *= (2.0 * props->visc / M_PI / props->sizes.x / props->sizes.z / props->kx / cos(props->alpha));
				
	return sum;	
}

void InclinedSum::prepare3D()
{
	double buf1, buf2;
	double sum_prev1 = 0.0, sum_prev2 = 0.0;
	int break_idx1, break_idx2;	
	int arr_idx;
	
	for(int seg_idx = 0; seg_idx < props->K; seg_idx++)
	{
		const Point& r = well->segs[seg_idx].r_bhp;
					
		for(int k = 0; k < props->K; k++)
		{
			arr_idx = seg_idx * props->K + k;
			F[arr_idx] = 0.0;
			
			break_idx2 = 0;
			
			for(int m = 1; m <= props->M; m++)
			{
				break_idx1 = 0;
				
				for(int l = 1; l <= props->L; l++)
				{
					buf1 = ((cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x - M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
								cos(M_PI * (double)(m) * well->segs[k].r2.x / props->sizes.x - M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							( M_PI * (double)(m) * tan(props->alpha) / props->sizes.x + M_PI * (double)(l) / props->sizes.z ) + 
								(cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x + M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
								cos(M_PI * (double)(m) * well->segs[k].r2.x / props->sizes.x + M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							( M_PI * (double)(m) * tan(props->alpha) / props->sizes.x - M_PI * (double)(l) / props->sizes.z )
								) / 2.0;
								
					buf2 = sqrt((double)(m) * (double)(m) / props->sizes.x / props->sizes.x + (double)(l) * (double)(l) / props->sizes.z / props->sizes.z);
				
					for(int i = -props->I; i <= props->I; i++)
					{	
						
						F[arr_idx] += buf1 / well->segs[k].length / buf2 *
							( exp(-M_PI * buf2 * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y)) - 
							exp(-M_PI * buf2 * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y)) ) *
							sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
							cos(M_PI * (double)(l) * r.z / props->sizes.z);
					}
					
					if( fabs(F[arr_idx] - sum_prev1) > F[arr_idx] * EQUALITY_TOLERANCE )
					{
						sum_prev1 = F[arr_idx];
						break_idx1 = 0;
					} else
						break_idx1++;
					
					if(break_idx1 > 1)
					{
						//std::cout << "l=" << l << std::endl;
						break;
					}
				}
				
				if( fabs(F[arr_idx] - sum_prev2) > F[arr_idx] * EQUALITY_TOLERANCE )
				{
					sum_prev2 = F[arr_idx];
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
}

double InclinedSum::getPres(int seg_idx)
{
	double s1, s2;
	
	s1 = get2D(seg_idx);
	s2 = get3D(seg_idx);
	
	return s1 + s2;
}
