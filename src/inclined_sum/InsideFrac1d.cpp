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
}
