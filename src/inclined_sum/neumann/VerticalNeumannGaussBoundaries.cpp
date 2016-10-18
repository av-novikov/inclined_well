#include "src/inclined_sum/neumann/VerticalNeumannGaussBoundaries.h"

VerticalNeumannGaussBoundaries::VerticalNeumannGaussBoundaries(const Parameters* _props, const Well* _well) : VerticalNeumann(_props, _well)
{
}

VerticalNeumannGaussBoundaries::~VerticalNeumannGaussBoundaries()
{
}

double VerticalNeumannGaussBoundaries::get2D(int seg_idx)
{
	const Point& r = well->segs[0].r_bhp;

	return F2d[0] + getBoundaries(r) - getPresAvg();
}

double VerticalNeumannGaussBoundaries::getPres(const Point& point)
{
	return F2d[0] + getBoundaries(point) - getPresAvg();
}

double VerticalNeumannGaussBoundaries::getBoundaries(const Point& point) const
{
	double sum = 0.0;
	
	sum += point.x * point.x / 2.0 / props->sizes.x * (props->fx2 - props->fx1) + props->fx1 * point.x +
			point.y * point.y / 2.0 / props->sizes.y * (props->fy2 - props->fy1) + props->fy1 * point.y;

	/*sum += sqrt(2.0 * M_PI * sigma * sigma) / props->sizes.x / props->sizes.y / props->sizes.z * 
												(erf(props->sizes.y / (2.0 * sqrt(2.0 * sigma * sigma))) * 
												(point.x * point.x * (Bx2 - Bx1) / 2.0 + props->sizes.x * point.x * Bx1) +
												erf(props->sizes.x / (2.0 * sqrt(2.0 * sigma * sigma))) *
												(point.y * point.y * (By2 - By1) / 2.0 + props->sizes.y * point.y * By1));

	double sumX = 0.0;
	double sumY = 0.0;
	for (int i = 1; i < props->M; i++)
	{
		sumX += cos(2.0 * M_PI * (double)(i)* point.x / props->sizes.x) / (double)(i * i * i) *
			(By2 * (exp(2.0 * M_PI * (double)(i) * props->sizes.y / props->sizes.x * (point.y / props->sizes.y - 1.0)) + 
					exp(2.0 * M_PI * (double)(i)* props->sizes.y / props->sizes.x * (-point.y / props->sizes.y - 1.0))) - 
			By1 * (exp(-2.0 * M_PI * (double)(i) / props->sizes.x * point.y) +
					exp(2.0 * M_PI * (double)(i)* props->sizes.y / props->sizes.x * (point.y / props->sizes.y - 2.0))));

		sumY += cos(2.0 * M_PI * (double)(i)* point.y / props->sizes.y) / (double)(i * i * i) *
			(Bx2 * (exp(2.0 * M_PI * (double)(i)* props->sizes.x / props->sizes.y * (point.x / props->sizes.x - 1.0)) +
					exp(2.0 * M_PI * (double)(i)* props->sizes.x / props->sizes.y * (-point.x / props->sizes.x - 1.0))) -
			Bx1 * (exp(-2.0 * M_PI * (double)(i) / props->sizes.y * point.x) +
					exp(2.0 * M_PI * (double)(i)* props->sizes.x / props->sizes.y * (point.x / props->sizes.x - 2.0))));
	}

	sumX *= (props->sizes.x * props->sizes.x * props->sizes.x / 4 / M_PI / M_PI / M_PI / sigma / sigma);
	sumY *= (props->sizes.y * props->sizes.y * props->sizes.y / 4 / M_PI / M_PI / M_PI / sigma / sigma);

	sum += sumX + sumY;*/

	return sum;
}

double VerticalNeumannGaussBoundaries::getPresAvg() const
{
	return 0.0;
}