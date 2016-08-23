#ifndef VERTICALDIRICHLET_H_
#define VERTICALDIRICHLET_H_

#include "src\inclined_sum\BaseSum.h"

class VerticalDirichlet : public BaseSum
{
protected:
	double directSum();
	double fourierSum();

public:
	VerticalDirichlet(const Parameters* _props, const Well* _well);
	~VerticalDirichlet();

	double get2D(int seg_idx);
	void prepare();
};

#endif /* VERTICALDIRICHLET_H_ */