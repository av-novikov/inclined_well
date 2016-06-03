#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <string>

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
	if(rank == 0)
	{
		cout << "Well coords:" << endl 
			<< "\t" << props->x_dim * props->r1
			<< "\t" << props->x_dim * props->r2;
		cout << "P_old= " << solver.getP_bhp() << endl;
	}
	
	MPI::Finalize();
	
	return 0;
}
