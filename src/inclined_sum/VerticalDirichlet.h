#ifndef VERTICALDIRICHLET_H_
#define VERTICALDIRICHLET_H_

#include "src/inclined_sum/BaseSum.h"

class VerticalDirichlet : public BaseSum
{
protected:
	double directSum();
	double fourierSum();

public:
	VerticalDirichlet(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~VerticalDirichlet();

	double get2D(int seg_idx);
	void prepare();

	double getPres(const Point& p);
	double getAnalyticalPres();
};

#endif /* VERTICALDIRICHLET_H_ */