#include <iostream>
#include <string>
#include <iomanip>

#include "paralution.hpp"

#include "src/WellFlow.h"
#include "src/utils/perf-utils.hpp"

#include "src/inclined_sum/InclinedSum.hpp"
#include "src/inclined_sum/Inclined3dSum.h"
#include "src/inclined_sum/VerticalDirichlet.h"
#include "src/inclined_sum/neumann/VerticalNeumann.h"
#include "src/inclined_sum/neumann/VerticalNeumannBoundaries.h"

using namespace std;
using namespace paralution;

void testNeumann()
{
	WellFlow solver("task/config2d.xml");
	double p_bhp;

	auto t = measure_time(
		[&]() {
		VerticalNeumannBoundaries inclSum(solver.getProps(), solver.getWell());
		solver.setSummator(&inclSum);
		const Parameters* props = solver.getProps();
		p_bhp = solver.getP_bhp() * props->p_dim / BAR;
	}, 1);

	print_test_results("PERF_TEST", t);

	cout << "P_bhp = " << p_bhp << endl;
}

void testNeumannBoundary()
{
	WellFlow solver("task/config2d.xml");
	double p_bhp;
	double p_avg;

	auto t = measure_time(
		[&]() {
		VerticalNeumannBoundaries inclSum(solver.getProps(), solver.getWell());
		solver.setSummator(&inclSum);
		const Parameters* props = solver.getProps();
		p_bhp = solver.getP_bhp() * props->p_dim / BAR;
		p_avg = inclSum.getPresAvg() * props->p_dim / BAR;
	}, 1);

	print_test_results("PERF_TEST", t);

	cout << "P_bhp = " << p_bhp << endl;
	cout << "P_avg = " << p_avg << endl;
}

void test3D()
{
	WellFlow solver("task/config3d.xml");
	double p_bhp;

	auto t = measure_time(
		[&]() {
		Inclined3dSum inclSum(solver.getProps(), solver.getWell());
		solver.setSummator(&inclSum);
		const Parameters* props = solver.getProps();
		p_bhp = solver.getP_bhp() * props->p_dim / BAR;
	}, 1);

	print_test_results("PERF_TEST", t);

	cout << "P_bhp = " << p_bhp << endl;
}

void testDirichlet()
{
	WellFlow solver("task/config3d.xml");
	double p_bhp;

	auto t = measure_time(
		[&]() {
		InclinedSum inclSum(solver.getProps(), solver.getWell());
		solver.setSummator(&inclSum);
		const Parameters* props = solver.getProps();
		p_bhp = solver.getP_bhp() * props->p_dim / BAR;
	}, 1);

	print_test_results("PERF_TEST", t);

	cout << "P_bhp = " << p_bhp << endl;
}

int main(int argc, char* argv[])
{
	init_paralution();

	testNeumannBoundary();
	//testDirichlet();
	//test3D();

	stop_paralution();

	while (true)
 	{
	}

	return 0;
}