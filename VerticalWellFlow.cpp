#include "VerticalWellFlow.hpp"

VerticalWellFlow::VerticalWellFlow(const Parameters& _props) : WellFlow(_props)
{
}

VerticalWellFlow::~VerticalWellFlow()
{
}

double VerticalWellFlow::calcPressure(const Point& r)
{
	double sum = 0.0;
	
	int m, n;
	
	for(m = 1; m <= props.M; m++)
		for(n = 1; n <= props.N; n++)
		{
			sum += sin(M_PI * m * r.x / props.sizes.x) * sin(M_PI * m * props.r1.x / props.sizes.x) * 
					sin(M_PI * n * r.y / props.sizes.y) * sin(M_PI * n * props.r1.y / props.sizes.y) /
					(m * m / props.sizes.x / props.sizes.x + n * n / props.sizes.y / props.sizes.y);
		}
	sum *= (4.0 * props.rate * props.visc / M_PI / M_PI / props.sizes.x / props.sizes.y / props.sizes.z / props.perm);			

	return sum;
}
