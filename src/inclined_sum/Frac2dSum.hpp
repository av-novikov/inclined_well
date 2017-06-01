#ifndef FRAC2DSUM_HPP_
#define FRAC2DSUM_HPP_

#include "src/inclined_sum/BaseSum.h"

class Frac2dSum : public BaseSum
{
protected:
	void prepareDirect();
	void prepareFourier();
public:
	Frac2dSum(const Parameters* _props, const Well* _well);
	~Frac2dSum();

	double getPres(const Point& point);
	double get3D(int seg_idx);
	double get2D(int seg_idx);
	void prepare();
};

#endif /* FRAC2DSUM_HPP_ */
