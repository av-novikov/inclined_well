#include <cassert>

#include "WellFlow.hpp"

WellFlow::WellFlow(const Parameters& _props) : props(_props)
{
	assert( fabs(props.alpha) > EQUALITY_TOLERANCE);
	
	well = new Well(props.r1, props.r2, props.K);
	well->setRate(props.rate);
	
	findRateDistribution();
}

WellFlow::~WellFlow()
{
}

void WellFlow::findRateDistribution()
{
	well->setUniformRate();
}

double WellFlow::highIntegral2D(const Point& r, double xi_c)
{
	double sum = 0.0;
	double buf;
	
	for(int k = 0; k < props.K; k++)
		for(int m = 1; m <= props.M; m++)
			for(int n = 1; n <= props.N; n++)
			{
				buf = ((double)(m) * (double)(m) / props.sizes.x / props.sizes.x + (double)(n) * (double)(n) / props.sizes.y / props.sizes.y);
				
				sum += well->segs[k].rate / well->segs[k].length / (double)(m) / buf * 
					sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x)) * 
					sin(M_PI * (double)(n) * r.y / props.sizes.y) * sin(M_PI * (double)(n) * props.r1.y / props.sizes.y) * 
					exp(-buf * M_PI * M_PI * xi_c);
			}
	sum *= (4.0 * props.visc / M_PI / M_PI / M_PI / props.sizes.y / props.sizes.z / props.perm / sin(props.alpha));
	
	return sum;
}

double WellFlow::lowIntegral2D(const Point& r, double xi_c)
{
	double sum = 0.0;
	double Iminus, Iplus;
	
	for(int k = 0; k < props.K; k++)
		for(int m = 1; m <= props.M; m++)
			for(int i = -props.I; i <= props.I; i++)
			{
				Iminus = exp(-M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
							(1.0 - erf((fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) - 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) + 
							exp(2.0 * M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
							(erf((fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) + 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) - 1.0));
				
				Iplus = exp(-M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
									(1.0 - erf((fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) - 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) + 
							exp(2.0 * M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
							(erf((fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) + 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) - 1.0));
				
				sum += well->segs[k].rate / well->segs[k].length / (double)(m) / (double)(m) * 
					( exp(-M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) - 
					  exp(-M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) ) *
					sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x)) * 
					(Iminus - Iplus);
					
				//std::cout << "Iminus = " << Iminus << "\t" << "Iplus = " << Iplus << std::endl;
				//std::cout << "sum = " << sum << std::endl;
			}
	sum *= (props.visc * props.sizes.x / 2.0 / M_PI / M_PI / props.sizes.z / props.perm / sin(props.alpha));
	
	return sum;
}

double WellFlow::calcPressure(const Point& r)
{	
	//double sum1 = 0.0;
	//double sum2 = 0.0;
	
	// l is equil to 0 case
	/*int k, m, n, i, l;
	
	for(k = 0; k < props.K; k++)
		for(m = 1; m <= props.M; m++)
			for(i = -props.I; i <= props.I; i++)
			{
				sum1 += well->segs[k].rate / well->segs[k].length / (double)(m) / (double)(m) * 
					( exp(-M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) - 
					  exp(-M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) ) *
					sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x));
			}
	sum1 *= (props.visc * props.sizes.x / M_PI / M_PI / props.sizes.z / props.perm / sin(props.alpha));
	
	for(k = 0; k < props.K; k++)
		for(m = 1; m <= props.M; m++)
			for(n = 1; n <= props.N; n++)
			{
				sum1 += well->segs[k].rate / well->segs[k].length / (double)(m) / 
				((double)(m) * (double)(m) / props.sizes.x / props.sizes.x + (double)(n) * (double)(n) / props.sizes.y / props.sizes.y) * 
				sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
				(cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x) -
				cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x)) * 
				sin(M_PI * (double)(n) * r.y / props.sizes.y) * sin(M_PI * (double)(n) * props.r1.y / props.sizes.y);
			}
	sum1 *= (4.0 * props.visc / M_PI / M_PI / M_PI / props.sizes.y / props.sizes.z / props.perm / sin(props.alpha));
	*/
	// l is not equil to 0 case
	/*double F;
	for(k = 0; k < props.K; k++)
		for(m = 1; m <= props.M; m++)
			for(l = 1; l <= props.L; l++)
				for(i = -props.I; i <= props.I; i++)
				{
					F = (cos(M_PI * (double)(m) / props->sizes.x * ()) / (M_PI * (double)(l) / props->sizes.z - M_PI * (double)(m) / props->sizes.x * tan(props->alpha))) / 2.0;
					sum2 += 
				}
	sum2 *= (8.0 * props.visc / props,sizes.x / props.sizes.y / props.sizes.z / props.perm / cos(props.alpha))*/
		
	return (lowIntegral2D(r, props.xi_c) + highIntegral2D(r, props.xi_c)) * props.p_dim;
}
