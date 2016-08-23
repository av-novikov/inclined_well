#ifndef VERTICALNEUMANN_H_
#define VERTICALNEUMANN_H_

#include "src\inclined_sum\BaseSum.h"

class VerticalNeumann : public BaseSum
{
protected:
	double directSum();
	double fourierSum();

public:
	VerticalNeumann(const Parameters* _props, const Well* _well);
	~VerticalNeumann();

	double get2D(int seg_idx);
	void prepare();
};

#endif /* VERTICALNEUMANN_H_ */