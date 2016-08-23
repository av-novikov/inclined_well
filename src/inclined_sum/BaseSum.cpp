#include "src/inclined_sum/BaseSum.h"

BaseSum::BaseSum(const Parameters* _props, const Well* _well) : props(_props), well(_well)
{
	size = props->K * props->K;
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