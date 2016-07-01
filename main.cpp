#include <iostream>
#include <unistd.h>
#include <string>
#include <iomanip>
#include <omp.h>

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
	int n = 1;
	
	while ((opt = getopt(argc, argv, "n:s")) != -1) {
		switch (opt) {
			case 'n' : n = atoi(optarg); break;
			case 's' : writeSnapshots = true; break;
			case '?' : fprintf(stderr, "Usage: inclined -s[to write snaps]\n"); return -1;
		}
	}
	
	cout << "Run programm with OpenMP, " << n << " thereads" << endl;
	omp_set_num_threads( n );
	
	WellFlow solver ("task/config.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator( &inclSum );

	const Parameters* props = solver.getProps();
	Point p = solver.getObsPoint();

	cout << "P_bhp = " << solver.getP_bhp() * props->p_dim / BAR << endl;		
	
	if( writeSnapshots )
	{
		Snapshotter snapshotter (props->sizes, props->nx, props->ny, props->nz, props->x_dim);
		snapshotter.setPresFoo( bind(&InclinedSum::getPres, &inclSum, _1) );
		snapshotter.snapshot("snap_0.vts");
	}

	return 0;
}
