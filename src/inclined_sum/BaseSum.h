#ifndef BASESUM_H_
#define BASESUM_H_

#define _USE_MATH_DEFINES
#include <math.h>

#include "src/Well.hpp"

class BaseSum
{
protected:
	const Parameters* props;
	const Well* well;

	int size;
	double* F2d;
	double* F3d;

public:
	BaseSum(const Parameters* _props, const Well* _well);
	virtual ~BaseSum();

	virtual void prepare() = 0;
	virtual double get2D(int seg_idx);
	virtual double get3D(int seg_idx);
	double getPres(int seg_idx);
};

#endif /* BASESUM_H_ */
