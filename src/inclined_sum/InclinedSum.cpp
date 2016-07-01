#include <iomanip>
#include <omp.h>
#include <new>

#include "src/inclined_sum/InclinedSum.hpp"

InclinedSum::InclinedSum(const Parameters* _props, const Well* _well) : props(_props), well(_well)
{
	F = new double** [props->K];
	buf = new double** [props->K];
	for(int i = 0; i < props->K; i++)
	{
		F[i] = new double* [props->M+1];
		buf[i] = new double* [props->M+1];
		
		for(int j = 0; j < props->M+1; j++)
		{
			F[i][j] = new double [props->L+1];
			buf[i][j] = new double [props->L+1];
		}
	}
	
	prepare3D();
}

InclinedSum::~InclinedSum()
{
	for(int i = 0; i < props->K; i++)
	{
		for(int j = 0; j < props->M+1; j++)
		{
			delete [] F[i][j];
			delete [] buf[i][j];
		}
		
		delete [] F[i];
		delete [] buf[i];
	}
	
	delete [] F;
	delete [] buf;
}

double InclinedSum::get2D(const Point& r)
{
	double sum = 0.0;
	double sum_prev = 0.0;
	int break_idx = 0;

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

double InclinedSum::get3D(const Point& r)
{
	double sum = 0.0;
	double sum_prev1 = 0.0, sum_prev2 = 0.0;
	int break_idx1, break_idx2;
	
	int k;
	
	#pragma omp parallel for reduction(+: sum) private(sum_prev1, sum_prev2, break_idx1, break_idx2) schedule(static) 
	for(k = 0; k < props->K; k++)
	{
		break_idx2 = 0;
		
		for(int m = 1; m <= props->M; m++)
		{
			break_idx1 = 0;
			
			for(int l = 1; l <= props->L; l++)
			{
				for(int i = -props->I; i <= props->I; i++)
				{	
					sum += F[k][m][l] * well->segs[k].rate / well->segs[k].length / buf[k][m][l] * 
						( exp(-M_PI * buf[k][m][l] * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y)) - 
						exp(-M_PI * buf[k][m][l] * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y)) ) *
						sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
						cos(M_PI * (double)(l) * r.z / props->sizes.z);
				}
				
				if( fabs(sum - sum_prev1) > sum * EQUALITY_TOLERANCE )
				{
					sum_prev1 = sum;
					break_idx1 = 0;
				} else
					break_idx1++;
				
				if(break_idx1 > 1)
				{
					//std::cout << "l=" << l << std::endl;
					break;
				}
			}
			
			if( fabs(sum - sum_prev2) > sum * EQUALITY_TOLERANCE )
			{
				sum_prev2 = sum;
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
				
	sum *= (2.0 * props->visc / M_PI / props->sizes.x / props->sizes.z / props->kx / cos(props->alpha));
				
	return sum;	
}

void InclinedSum::prepare3D()
{
	for(int k = 0; k < props->K; k++)
		for(int m = 1; m <= props->M; m++)
		{			
			for(int l = 1; l <= props->L; l++)
			{
				buf[k][m][l] = sqrt((double)(m) * (double)(m) / props->sizes.x / props->sizes.x + (double)(l) * (double)(l) / props->sizes.z / props->sizes.z);
			
				F[k][m][l] = ((cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x - M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
							cos(M_PI * (double)(m) * well->segs[k].r2.x / props->sizes.x - M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
						( M_PI * (double)(m) * tan(props->alpha) / props->sizes.x + M_PI * (double)(l) / props->sizes.z ) + 
							(cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x + M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
							cos(M_PI * (double)(m) * well->segs[k].r2.x / props->sizes.x + M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
						( M_PI * (double)(m) * tan(props->alpha) / props->sizes.x - M_PI * (double)(l) / props->sizes.z )
							) / 2.0;
			}
		}
}

double InclinedSum::getPres(const Point& r)
{
	double s1, s2;
	
	s1 = get2D(r);
	s2 = get3D(r);
	
	std::cout << std::setprecision(10) << props->x_dim * r << "2d = " << s1 << "\t3d = " << s2 << std::endl;
	return s1 + s2;
}
