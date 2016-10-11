#ifndef VERTICALNEUMANNUNIFORMBOUNDARIES_H_
#define VERTICALNEUMANNUNIFORMBOUNDARIES_H_

#include "src/inclined_sum/neumann/VerticalNeumann.h"

class VerticalNeumannUniformBoundaries : public VerticalNeumann
{
protected:
	double getBoundaries(const Point& point) const;

public:
	VerticalNeumannUniformBoundaries(const Parameters* _props, const Well* _well);
	~VerticalNeumannUniformBoundaries();

	double getPresAvg() const;

	double getPres(const Point& point);
	double get2D(int seg_idx);
};

#endif /* VERTICALNEUMANNUNIFORMBOUNARIES_H_ */