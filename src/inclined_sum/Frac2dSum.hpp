#ifndef FRAC2DSUM_HPP_
#define FRAC2DSUM_HPP_

#include "src/inclined_sum/BaseSum.h"
#include "math/integral.h"

class Frac2dSum : public BaseSum
{
protected:
	void prepareDirect();
	void prepareFourier();

	// For RV3D numerical integration
	Integral* integr;
	Double2d* integrand;
	static const int PART_SIZE = 20;
public:
	Frac2dSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~Frac2dSum();

	double get3D(int seg_idx) const;
	double get2D(int seg_idx) const;
	void prepare();

	double getPressure(const Point& r);
};

#endif /* FRAC2DSUM_HPP_ */
