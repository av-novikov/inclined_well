#include <string>

#include "src/inclined_sum/BaseSum.h"

using std::stoi;
using std::stod;

BaseSum::BaseSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well) : sprops(_sprops), props(_props), well(_well), gprops(_well->getGeomProps())
{
	size = sprops.K * sprops.K;
	F2d = new double[size];
	F3d = new double[size];
}
BaseSum::~BaseSum()
{
	delete[] F2d;
	delete[] F3d;
}
double BaseSum::getPres(int seg_idx)
{
	return get2D(seg_idx) + get3D(seg_idx);
}
double BaseSum::get2D(int seg_idx)
{
	return 0.0;
}
double BaseSum::get3D(int seg_idx)
{
	return 0.0;
}
const SummatorProperties* BaseSum::getSumProps() const
{
	return &sprops;
}
const Well* BaseSum::getWell() const
{
	return well;
}
Well* BaseSum::getWell()
{
	return const_cast<Well*>(well);
}