#ifndef VERTICALNEUMANN_H_
#define VERTICALNEUMANN_H_

#include "src/inclined_sum/BaseSum.h"

class VerticalNeumann : public BaseSum
{
protected:
	double directSum(const Point& r);
	double fourierSum(const Point& r);

public:
	VerticalNeumann(const Parameters* _props, const Well* _well);
	~VerticalNeumann();

	double getPres(const Point& point);
	double get2D(int seg_idx);
	void prepare();
};

#endif /* VERTICALNEUMANN_H_ */