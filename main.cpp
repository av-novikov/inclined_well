#include <iostream>
#include <string>
#include <iomanip>
#include <omp.h>

#include "src/WellFlow.hpp"
#include "src/VerticalWellFlow.hpp"
#include "src/inclined_sum/InclinedSum.hpp"

using namespace std;

int main(int argc, char* argv[])
{
	WellFlow solver ("task/config.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator( &inclSum );

	const Parameters* props = solver.getProps();

	cout << "P_bhp = " << solver.getP_bhp() * props->p_dim / BAR << endl;		

	return 0;
}
