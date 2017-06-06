#ifndef WELLFLOW_H_
#define WELLFLOW_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <vector>
#include <cassert>

#include "paralution.hpp"
#include "tinyxml2.h"

#include "src/Well.hpp"
#include "src/Fracture.hpp"
#include "src/inclined_sum/BaseSum.h"

class WellFlow
{
protected:
	std::vector<Well> wells;
	std::vector<Fracture> fractures;
	std::vector<WellSegment*> segs;
	std::vector<BaseSum*> summators;
	MainProperties props;

	// For linear system treating
	int seg_num;
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
	void load(const std::string fileName);
	void loadGeomProps(tinyxml2::XMLElement* xml_well, WellGeomProperties* geom_props);
	void loadSumProps(tinyxml2::XMLElement* xml_summator, SummatorProperties* sum_props);

	bool firstTime;

	double pres_av, pres_dev;
public:
	WellFlow(const std::string fileName);
	~WellFlow();

	const MainProperties* getProps() const;
	const Well* getWell(const int i) const;
	const BaseSum* getSummator(const int i) const;

	void calcPressure();
	double getP_bhp();
};

#endif /* WELLFLOW_H_ */
