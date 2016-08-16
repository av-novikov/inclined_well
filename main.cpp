#include <iostream>
#include <iomanip>

#include "paralution.hpp"

#include "src/WellFlow.h"
#include "src/utils/perf-utils.hpp"

using namespace std;
using namespace paralution;

int main()
{
	typedef double T;

	init_paralution();

	WellFlow<T> solver("task/config.xml");
	double p_bhp;

	auto t = measure_time(
		[&]() {
			Inclined3dSum<T> inclSum(solver.getProps(), solver.getWell());
			solver.setSummator(&inclSum);
			const Parameters<T>* props = solver.getProps();
			p_bhp = solver.getP_bhp() * props->p_dim / BAR;
		}, 1);
	stop_paralution();

	print_test_results("CUDA_PERF_TEST", t);

	cout << setprecision(6);
	cout << "P_bhp = " << p_bhp << endl;

	while(true)
	{ }

	return 0;
}