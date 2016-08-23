#ifndef WELLFLOW_H_
#define WELLFLOW_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>

#include "paralution.hpp"

#include "src/Well.hpp"
#include "src/inclined_sum/BaseSum.h"

class WellFlow
{
protected:
	Parameters props;
	Well* well;
	BaseSum* inclSum;

	int matSize;
	paralution::Inversion<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double> ls;

	double* q;
	double* dq;
	double* b;
	double* x;
	double* a;
	double** dpdq;

	int* ind_i;
	int* ind_j;

	paralution::LocalMatrix<double> A;
	paralution::LocalVector<double> X;
	paralution::LocalVector<double> B;

	void findRateDistribution();
	void loadTask(const std::string fileName);

public:
	WellFlow(const std::string fileName);
	~WellFlow();

	void setSummator(BaseSum* _inclSum);
	const Parameters* getProps() const;
	const Well* getWell() const;

	void calcPressure();
	double getP_bhp();
};

#endif /* WELLFLOW_H_ */
