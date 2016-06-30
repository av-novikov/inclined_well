#include <iostream>
#include <string>
#include <mpi.h>

#include "src/utils/perf-utils.hpp"
#include "src/WellFlow.hpp"
#include "src/inclined_sum/InclinedSum.hpp"

#define BAR 1.E+5

using std::cout;
using std::endl;
using std::to_string;

int main(int argc, char* argv[])
{
	MPI::Init();
	const int rank = MPI::COMM_WORLD.Get_rank();
	
	WellFlow solver ("task/vertical.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator( &inclSum );

	const Parameters* props = solver.getProps();
	Point p = solver.getObsPoint();

	if(rank == 0)
	{
		print_test_title("vertical test");
		cout << "Observation point: " << props->x_dim * p;
	}

	double p_bhp;
	const double p_bhp_true = 38.9687;
	auto t = measure_time(
        [&](){ p_bhp = solver.getP_bhp() * props->p_dim / BAR; },
        1
    );

	if(rank == 0)
	{
		cout << "P_bhp = " << p_bhp << endl;
		cout << "P_bhp_true = " << p_bhp_true << endl;
		if( fabs(p_bhp - p_bhp_true) > 0.5 )
			cout << "TEST FAILED !!!" << endl;
			
		print_test_results("vertical test", t);
	}
	
	MPI::Finalize();
	
	return 0;
}
