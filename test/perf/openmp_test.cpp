#include <iostream>
#include <unistd.h>
#include <string>
#include <cstdlib>
#include <omp.h>

#include "src/WellFlow.hpp"
#include "src/inclined_sum/InclinedSum.hpp"
#include "src/utils/perf-utils.hpp"

#define BAR 1.E+5

using std::cout;
using std::endl;
using std::to_string;

int main(int argc, char* argv[])
{
	int opt = 0;
	int n = 1;
	
	while ((opt = getopt(argc, argv, "n:")) != -1) {
		switch (opt) {
			case 'n' : n = atoi(optarg); break;
			case '?' : fprintf(stderr, "Usage: inclined -n [cores num]\n"); return -1;
		}
	}
	
	WellFlow solver ("task/config.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator( &inclSum );
	const Parameters* props = solver.getProps();

	double p1, p2;
	auto t = measure_time2(
		[&](){ DUMMY_FUNC; },
        [&](){ p1 = solver.getP_bhp() * props->p_dim / BAR; },
        [&](){ omp_set_num_threads( n ); },
        [&](){ p2 = solver.getP_bhp() * props->p_dim / BAR; },
        1
    );

	cout << "p[1 thread] = " << p1 << endl;		
	cout << "p[" << n << " threads] = " << p2 << endl;		

	print_test_results("1 thread", t.first, to_string(n) + " threads", t.second);

	return 0;
}
