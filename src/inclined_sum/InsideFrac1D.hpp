#ifndef INSIDEFRAC1D_HPP_
#define INSIDEFRAC1D_HPP_

#include "src/inclined_sum/BaseSum.h"

class InsideFrac1d : public BaseSum
{
protected:
	double p1, p2;
	int mid_idx;
public:
	InsideFrac1d(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~InsideFrac1d();

	double get3D(int seg_idx);
	double get2D(int seg_idx);
	void prepare();
	void setBounds(double _p1, double _p2);
};

#endif /* INSIDEFRAC1D_HPP_ */
