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
	return 0.0;
	//return point.x * point.x / 2.0 / props->sizes.x * (props->fx2 - props->fx1) + props->fx1 * point.x +
	//	point.y * point.y / 2.0 / props->sizes.y * (props->fy2 - props->fy1) + props->fy1 * point.y;
}

double VerticalNeumannGaussBoundaries::getPresAvg() const
{
	return 0.0;
}