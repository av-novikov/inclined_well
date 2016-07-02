#include <iostream>
#include <string>

#include "src/utils/perf-utils.hpp"
#include "src/WellFlow.hpp"
#include "src/inclined_sum/InclinedSum.hpp"

#define BAR 1.E+5

using std::cout;
using std::endl;
using std::to_string;

int main(int argc, char* argv[])
{
	WellFlow solver ("task/vertical.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator( &inclSum );

	const Parameters* props = solver.getProps();

	print_test_title("vertical test");

	double p_bhp;
	const double p_bhp_true = 38.9687;
	auto t = measure_time(
        [&](){ p_bhp = solver.getP_bhp() * props->p_dim / BAR; },
        1
    );

	cout << "P_bhp = " << p_bhp << endl;
	cout << "P_bhp_true = " << p_bhp_true << endl;
	if( fabs(p_bhp - p_bhp_true) > 0.5 )
		cout << "TEST FAILED !!!" << endl;
			
	print_test_results("vertical test", t);
	
	return 0;
}
