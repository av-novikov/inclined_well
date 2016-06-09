#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <string>
#include <iomanip>

#include "src/WellFlow.hpp"
#include "src/VerticalWellFlow.hpp"
#include "src/Snapshotter.hpp"

#include "src/inclined_sum/InclinedSum.hpp"

using namespace std;
using namespace std::placeholders;

int main(int argc, char* argv[])
{
	int opt = 0;
	bool writeSnapshots = false;
	
	while ((opt = getopt(argc, argv, "s")) != -1) {
		switch (opt) {
			case 's' : writeSnapshots = true; break;
			case '?' : fprintf(stderr, "Usage: inclined -s[to write snaps]\n"); return -1;
		}
	}
	
	MPI::Init();
	const int rank = MPI::COMM_WORLD.Get_rank();
	//const int size = MPI::COMM_WORLD.Get_size();	
	
	WellFlow solver ("task/config.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator( &inclSum );

	const Parameters* props = solver.getProps();
	Point p = solver.getObsPoint();

	if(rank == 0)
	{
		cout << "P_bhp = " << solver.getP_bhp() * props->p_dim / BAR << endl;
	}
	
	MPI::Finalize();
	
	return 0;
}
