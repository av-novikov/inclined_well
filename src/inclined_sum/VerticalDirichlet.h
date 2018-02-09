#ifndef VERTICALDIRICHLET_H_
#define VERTICALDIRICHLET_H_

#include "src/inclined_sum/BaseSum.h"

class VerticalDirichlet : public BaseSum
{
protected:
	void directSum();
	void fourierSum();

public:
	VerticalDirichlet(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~VerticalDirichlet();

	double get2D(int seg_idx) const;
	double get3D(int seg_idx) const;
	void prepare();

	double getPressure(const Point& point);
	double getAnalyticalPres(const int seg_idx) const;
};

#endif /* VERTICALDIRICHLET_H_ */