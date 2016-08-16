#ifndef WELLFLOW_H_
#define WELLFLOW_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>

#include "paralution.hpp"

#include "src/Well.cuh"
#include "src/summators/InclinedSum.cuh"
#include "src/summators/Inclined3dSum.cuh"

template <class T>
class WellFlow
{
protected:
	Parameters<T> props;
	Well<T>* well;
	BaseSum<T>* inclSum;

	int matSize;
	paralution::Inversion<paralution::LocalMatrix<T>, paralution::LocalVector<T>, T> ls;

	T* q;
	T* dq;
	T* b;
	T* x;
	T* a;
	T** dpdq;

	int* ind_i;
	int* ind_j;

	paralution::LocalMatrix<T> A;
	paralution::LocalVector<T> X;
	paralution::LocalVector<T> B;

	void findRateDistribution();
	void loadTask(const std::string fileName);

public:
	WellFlow(const std::string fileName);
	~WellFlow();

	void setSummator(BaseSum<T>* _inclSum);
	const Parameters<T>* getProps() const;
	const Well<T>* getWell() const;

	void calcPressure();
	T getP_bhp();
};

#endif /* WELLFLOW_H_ */
