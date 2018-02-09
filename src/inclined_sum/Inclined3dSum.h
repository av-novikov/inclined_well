#ifndef INCLINED3DSUM_H_
#define INCLINED3DSUM_H_

#include "src/inclined_sum/BaseSum.h"
#include "math/integral.h"

class Inclined3dSum : public BaseSum
{
protected:
	void prepareDirect();
	void prepareFourier();

	// For RV3D numerical integration
	Integral* integr;
	Double2d* integrand;
	static const int PART_SIZE = 10;
public:
	Inclined3dSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~Inclined3dSum();

	void prepare();
	double getPressure(const Point& point);
	double get2D(int seg_idx) const;
	double get3D(int seg_idx) const;
};

#endif /* INCLINED3DSUM_H_ */

