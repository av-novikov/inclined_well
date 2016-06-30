#include <iostream>
#include <mpi.h>

#include "src/utils/perf-utils.hpp"
#include "src/WellFlow.hpp"
#include "src/inclined_sum/InclinedSum.hpp"

#define BAR 1.E+5

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
	MPI::Init();
	const int rank = MPI::COMM_WORLD.Get_rank();
	
	WellFlow solver ("task/2d.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator( &inclSum );

	const Parameters* props = solver.getProps();
	Point p = solver.getObsPoint();

	if(rank == 0)
	{
		print_test_title("2d VS 3d");
		cout << "Observation point: " << props->x_dim * p;
	}

	double p_2d, p_3d;
	auto t = measure_time2(
        [&](){ p_2d = inclSum.get2D(p) * props->p_dim / BAR; },
        [&](){ p_3d = inclSum.get3D(p) * props->p_dim / BAR; },
        1
    );

	if(rank == 0)
	{
		cout << "p_2d = " << p_2d << endl;
		cout << "p_3d = " << p_3d << endl;
		cout << "sum  = " << p_2d + p_3d << endl;
		print_test_results("2D", t.first, "3D", t.second);
	}
	
	MPI::Finalize();
	
	return 0;
}
