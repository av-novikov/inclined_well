#ifndef GRID_HPP_
#define GRID_HPP_

#include <string>
#include <functional>

#include "Well.hpp"

typedef std::function<double(Point&, bool)> FooType;

class Grid
{
protected:
		Point sizes;
		Point steps;
		int nx, ny, nz;
		double x_dim;
		 
		int rank;
		int size;
		
		int startIdx;
		int finishIdx;
		
		FooType presFoo;
		
public:
		Grid(const Point& _sizes, const int _nx, const int _ny, const int _nz, double _x_dim);
		~Grid();
		
		void setPresFoo(const FooType& foo);
		void snapshot(std::string name);
};

#endif /* GRID_HPP_ */
