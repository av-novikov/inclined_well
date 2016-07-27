#include <iostream>

#include "src/WellFlow.h"

using namespace std;

int main()
{
	WellFlow solver ("task/config.xml");
	InclinedSum inclSum( solver.getProps(), solver.getWell() );
	solver.setSummator(&inclSum);

	const Parameters* props = solver.getProps();
	double p_bhp = solver.getP_bhp() * props->p_dim / BAR;

	while (true)
	{
	}

	return 0;
}