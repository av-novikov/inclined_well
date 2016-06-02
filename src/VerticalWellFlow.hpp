#ifndef VERTICALWELLFLOW_HPP_
#define VERTICALWELLFLOW_HPP_

#include "src/WellFlow.hpp"

class VerticalWellFlow : public WellFlow
{
public:
	VerticalWellFlow(const std::string fileName);
	~VerticalWellFlow();
	
	//double calcPressure(const Point& r);
	
};

#endif /* VERTICALWELLFLOW_HPP_ */
