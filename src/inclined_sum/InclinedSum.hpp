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
	
	double** F;
	//double*** buf;
	
	void prepare3D();
		
public:
	InclinedSum(const Parameters* _props, const Well* _well);
	~InclinedSum();
		
	virtual double getPres(int seg_idx);
	virtual double get2D(int seg_idx);
	virtual double get3D(int seg_idx);
};

#endif /* INCLINEDSUM_HPP_ */
