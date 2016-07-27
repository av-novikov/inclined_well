#ifndef WELLFLOW_H_
#define WELLFLOW_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>

#include "src/Well.cuh"
#include "src/InclinedSum.cuh"

class WellFlow
{
protected:
	Parameters props;
	Well* well;

	void findRateDistribution();
	void loadTask(const std::string fileName);

public:
	WellFlow(const std::string fileName);
	~WellFlow();

	InclinedSum* inclSum;

	void setSummator(InclinedSum* _inclSum);
	const Parameters* getProps() const;
	const Well* getWell() const;

	double getP_bhp();
};

#endif /* WELLFLOW_H_ */