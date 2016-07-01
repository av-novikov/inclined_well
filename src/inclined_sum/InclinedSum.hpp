#ifndef INCLINEDSUM_HPP_
#define INCLINEDSUM_HPP_

#include <cmath>
#include <mpi.h>

#include "src/Well.hpp"

class InclinedSum
{
protected:
	const Parameters* props;
	const Well* well;
	
	double*** F;
	double*** buf;
	
	void prepare3D();
		
public:
	InclinedSum(const Parameters* _props, const Well* _well);
	~InclinedSum();
		
	virtual double getPres(const Point& r);
	virtual double get2D(const Point& r);
	virtual double get3D(const Point& r);
};

#endif /* INCLINEDSUM_HPP_ */
