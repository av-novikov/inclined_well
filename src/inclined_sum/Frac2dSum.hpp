#ifndef FRAC2DSUM_HPP_
#define FRAC2DSUM_HPP_

#include "src/inclined_sum/BaseSum.h"

class Frac2dSum : public BaseSum
{
protected:
	double fourierSum(const Point& r);

public:
	Frac2dSum(const Parameters* _props, const Well* _well);
	~Frac2dSum();

	double getPres(const Point& point);
	double get2D(int seg_idx);
	void prepare();
};

#endif /* FRAC2DSUM_HPP_ */
