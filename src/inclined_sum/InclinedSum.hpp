#ifndef INCLINEDSUM_HPP_
#define INCLINEDSUM_HPP_

#include "src/inclined_sum/BaseSum.h"

class InclinedSum : public BaseSum
{
protected:	
	void prepare2D();
	void prepare3D();
		
public:
	InclinedSum(const Parameters* _props, const Well* _well);
	~InclinedSum();
	
	void prepare();
	double getPres(const Point& point);
	double get2D(int seg_idx);
	double get3D(int seg_idx);
};

#endif /* INCLINEDSUM_HPP_ */
