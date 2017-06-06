#include "src/inclined_sum/neumann/VerticalNeumannUniformBoundaries.h"

VerticalNeumannUniformBoundaries::VerticalNeumannUniformBoundaries(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well) : VerticalNeumann(_sprops, _props, _well)
{
}
VerticalNeumannUniformBoundaries::~VerticalNeumannUniformBoundaries()
{
}
double VerticalNeumannUniformBoundaries::get2D(int seg_idx)
{
	const Point& r = well->segs[0].r_bhp;

	return F2d[0] + getBoundaries(r) - getPresAvg();
}
double VerticalNeumannUniformBoundaries::getPres(const Point& point)
{
	return F2d[0] + getBoundaries(point) - getPresAvg();
}
double VerticalNeumannUniformBoundaries::getBoundaries(const Point& point) const
{
	return point.x * point.x / 2.0 / props->sizes.x * (props->fx2 - props->fx1)  + props->fx1 * point.x +
			point.y * point.y / 2.0 / props->sizes.y * (props->fy2 - props->fy1) + props->fy1 * point.y;
}
double VerticalNeumannUniformBoundaries::getPresAvg() const
{
	return (props->fx2 + 2.0 * props->fx1) * props->sizes.x / 6.0 + (props->fy2 + 2.0 * props->fy1) * props->sizes.y / 6.0;
}