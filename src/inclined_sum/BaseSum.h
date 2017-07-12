#ifndef BASESUM_H_
#define BASESUM_H_

#define _USE_MATH_DEFINES
#include <math.h>

#include "src/Well.hpp"
#include "src/Properties.hpp"

class BaseSum
{
protected:
	const SummatorProperties sprops;
	const MainProperties* props;
	const WellGeomProperties* gprops;
	const Well* well;
	std::vector<WellSegment*>* segs;

	int size;
	double* F2d;
	double* F3d;
public:
	BaseSum(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well);
	virtual ~BaseSum();

	const SummatorProperties* getSumProps() const;
	const Well* getWell() const;
	void setSegments(std::vector<WellSegment*>* _segs);
	Well* getWell();

	virtual void prepare() = 0;
	virtual double get2D(int seg_idx);
	virtual double get3D(int seg_idx);
	double getPres(int seg_idx);
	virtual double getAnalyticalPres() const;

	virtual double getPressure(const Point& r);
};

#endif /* BASESUM_H_ */
