#ifndef VERTICALNEUMANNBOUNDARIES_H_
#define VERTICALNEUMANNBOUNDARIES_H_

#include "src/inclined_sum/neumann/VerticalNeumann.h"

class VerticalNeumannBoundaries : public VerticalNeumann
{
protected:
	double getBoundaries(const Point& point) const;

public:
	VerticalNeumannBoundaries(const Parameters* _props, const Well* _well);
	~VerticalNeumannBoundaries();

	double getPresAvg() const;

	double getPres(const Point& point);
	double get2D(int seg_idx);
};

#endif /* VERTICALNEUMANNBOUNARIES_H_ */