#ifndef INCLINED3DSUM_H_
#define INCLINED3DSUM_H_

#include "src/inclined_sum/BaseSum.h"

class Inclined3dSum : public BaseSum
{
protected:
	void prepare2D();
	void prepare3D();

public:
	Inclined3dSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~Inclined3dSum();

	void prepare();
	double getPres(const Point& point);
	double get2D(int seg_idx);
	double get3D(int seg_idx);
};

#endif /* INCLINED3DSUM_H_ */

