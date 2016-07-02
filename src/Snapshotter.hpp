#ifndef SNAPSHOTTER_HPP_
#define SNAPSHOTTER_HPP_

#include <string>
#include <functional>

#include "src/Well.hpp"

typedef std::function<double(int)> FooType;

class Snapshotter
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
		Snapshotter(const Point& _sizes, const int _nx, const int _ny, const int _nz, double _x_dim);
		~Snapshotter();
		
		void setPresFoo(const FooType& foo);
		void snapshot(std::string name);
};

#endif /* SNAPSHOTTER_HPP_ */
