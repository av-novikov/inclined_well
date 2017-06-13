#ifndef VERTICALNEUMANN_H_
#define VERTICALNEUMANN_H_

#include "src/inclined_sum/BaseSum.h"

class VerticalNeumann : public BaseSum
{
protected:
	double directSum(const Point& r);
	double fourierSum(const Point& r);

public:
	VerticalNeumann(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~VerticalNeumann();

	double get2D(int seg_idx);
	void prepare();
};

#endif /* VERTICALNEUMANN_H_ */