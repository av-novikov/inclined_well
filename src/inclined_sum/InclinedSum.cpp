#include "src/inclined_sum/InclinedSum.hpp"

InclinedSum::InclinedSum(const Parameters* _props, const Well* _well) : props(_props), well(_well)
{
}

InclinedSum::~InclinedSum()
{
}

double InclinedSum::get2D(const Point& r)
{
	double sum = 0.0;
	
	for(int k = 0; k < props->K; k++)
		for(int m = 1; m <= props->M; m++)
			for(int i = -props->I; i <= props->I; i++)
			{
				sum += /*well->segs[k].rate / well->segs[k].length*/ 1.0 / (double)(m) / (double)(m) * 
					( exp(-M_PI * (double)(m) * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) - 
					exp(-M_PI * (double)(m) * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) ) *
					sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props->alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props->sizes.x));
			}
	//sum *= (props->visc * props->sizes.x / M_PI / M_PI / props->sizes.z / props->perm / sin(props->alpha));
	//sum *= (well->segs[0].length);
		
	return sum;	
}

double InclinedSum::get3D(const Point& r)
{
	double sum = 0.0;
	double F, buf;
	
	for(int k = 0; k < props->K; k++)
		for(int m = 1; m <= props->M; m++)
			for(int l = 1; l <= props->L; l++)
			{
				buf = sqrt((double)(m) * (double)(m) / props->sizes.x / props->sizes.x + (double)(l) * (double)(l) / props->sizes.z / props->sizes.z);
				
				for(int i = -props->I; i <= props->I; i++)
				{
					F = ((cos(M_PI * (double)(m) * (well->segs[k].r1.x + 
							tan(props->alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props->sizes.x + 
							M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
						cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x +
							M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							(M_PI * (double)(l) / props->sizes.z - M_PI * (double)(m) * tan(props->alpha) / props->sizes.x) - 
						(cos(M_PI * (double)(m) * (well->segs[k].r1.x + 
							tan(props->alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props->sizes.x - 
							M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
						cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x -
							M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							(M_PI * (double)(l) / props->sizes.z + M_PI * (double)(m) * tan(props->alpha) / props->sizes.x)
							) / 2.0;
						
						sum += well->segs[k].rate / well->segs[k].length / buf * 
							( exp(-M_PI * buf * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y)) - 
							exp(-M_PI * buf * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y)) ) *
							sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
							cos(M_PI * (double)(l) * r.z / props->sizes.z);
				}
			}
				
	sum *= (2.0 * props->visc / M_PI / props->sizes.x / props->sizes.z / props->perm / cos(props->alpha));
				
	return sum;	
}

double InclinedSum::getPres(const Point& r)
{
	return get2D(r);
	//return ( get2D(r) + get3D(r) ) * props->p_dim;
}
