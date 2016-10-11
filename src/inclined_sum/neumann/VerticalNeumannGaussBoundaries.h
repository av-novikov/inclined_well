#ifndef VERTICALNEUMANNGAUSSBOUNDARIES_H_
#define VERTICALNEUMANNGAUSSBOUNDARIES_H_

#include "src/inclined_sum/neumann/VerticalNeumann.h"

class VerticalNeumannGaussBoundaries : public VerticalNeumann
{
protected:
	double getBoundaries(const Point& point) const;

public:
	VerticalNeumannGaussBoundaries(const Parameters* _props, const Well* _well);
	~VerticalNeumannGaussBoundaries();

	double getPresAvg() const;

	double getPres(const Point& point);
	double get2D(int seg_idx);
};

#endif /* VERTICALNEUMANNGAUSSBOUNARIES_H_ */
