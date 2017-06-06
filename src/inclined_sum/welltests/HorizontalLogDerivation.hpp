#ifndef HORIZONTALLOGDERIVATION_HPP_
#define HORIZONTALLOGDERIVATION_HPP_

#include "src/inclined_sum/BaseSum.h"

class HorizontalLogDerivation : public BaseSum
{
protected:
	double time;
public:
	HorizontalLogDerivation(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	~HorizontalLogDerivation();

	void prepare();
	double get2D(int seg_idx);
	double get3D(int seg_idx);
	double getPres(const Point& p);
	double getLogDerivative();

	void setTime(const double _time);
};

#endif /* HORIZONTALLOGDERIVATION_HPP_ */