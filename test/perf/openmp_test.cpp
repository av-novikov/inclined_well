#include <iostream>
#include <string>

#include "src/WellFlow.hpp"
#include "src/inclined_sum/InclinedSum.hpp"
#include "src/utils/perf-utils.hpp"

#define BAR 1.E+5

using std::cout;
using std::endl;
using std::to_string;

int main(int argc, char* argv[])
{
	WellFlow solver ("task/config.xml");

	double p1, p2;
	auto t = measure_time(
        [&](){ 
				InclinedSum inclSum( solver.getProps(), solver.getWell() );
				solver.setSummator( &inclSum );
				const Parameters* props = solver.getProps();
				p1 = solver.getP_bhp() * props->p_dim / BAR; 
			},
        1
    );

	cout << "p = " << p1 << endl;	
		
	print_test_results("perf test", t);

	return 0;
}
