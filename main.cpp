#include <iostream>
#include <string>
#include <iomanip>

#include "paralution.hpp"

#include "src/WellFlow.h"
#include "src/utils/perf-utils.hpp"

#include "src/inclined_sum/InclinedSum.hpp"
#include "src/inclined_sum/Inclined3dSum.h"
#include "src/inclined_sum/VerticalDirichlet.h"
#include "src/inclined_sum/VerticalNeumann.h"

using namespace std;
using namespace paralution;

int main(int argc, char* argv[])
{
	init_paralution();

	WellFlow solver ("task/config2d.xml");
	double p_bhp;

	auto t = measure_time(
		[&]() {
			VerticalNeumann inclSum(solver.getProps(), solver.getWell());
			solver.setSummator(&inclSum);
			const Parameters* props = solver.getProps();
			p_bhp = solver.getP_bhp() * props->p_dim / BAR;
		}, 1);
	stop_paralution();

	print_test_results("PERF_TEST", t);

	cout << "P_bhp = " << p_bhp << endl;		

	while (true)
 	{
	}

	return 0;
}
