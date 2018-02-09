#ifndef INCLINEDSUM_HPP_
#define INCLINEDSUM_HPP_

#include "src/inclined_sum/BaseSum.h"

class InclinedSum : public BaseSum
{
protected:	
	void prepare2D();
	void prepare3D();
		
public:
	InclinedSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~InclinedSum();
	
	void prepare();
	double get2D(int seg_idx) const;
	double get3D(int seg_idx) const;
};

#endif /* INCLINEDSUM_HPP_ */
