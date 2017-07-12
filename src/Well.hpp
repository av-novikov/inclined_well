#ifndef WELL_HPP_
#define WELL_HPP_

#include <iostream>
#include <vector>
#include <algorithm>
#include <new>

#include "src/Properties.hpp"
#include "src/utils/Point.hpp"

#define BAR 1.E+5

enum WellType {VERTICAL, VERTICAL_NEUMANN, SLANTED, HORIZONTAL, FRAC};
struct WellSegment
{
	const int well_idx;
	const int seg_idx;

	const Point r1;
	const Point r2;
	const Point r_bhp;
	
	double tau1, tau2;
	double length;
	double pres;
	double pres2D;
	double pres3D;
	double rate;
	
	WellSegment(const Point& _r1, const Point& _r2, const Point _r_bhp, 
				double _tau1, double _tau2, const int _well_idx, const int _seg_idx) : r1(_r1), r2(_r2), r_bhp(_r_bhp), 
																						tau1(_tau1), tau2(_tau2), well_idx(_well_idx), seg_idx(_seg_idx)
	{
		length = sqrt((r2 - r1) * (r2 -r1));
		pres = pres2D = pres3D = rate = 0.0;
	};
};
class Well
{
	protected:
		const int well_idx;
		const WellGeomProperties props;
		const std::string name;
		const int num;
		double rate;
	public:
		Well(const WellGeomProperties& _props, const WellType& _type, const std::string _name, const int _well_idx);
		~Well();
		
		const WellType type;

		void setRate(double _rate);
		void setUniformRate();
		void setGaussRate(const double sigma);
		
		std::vector<WellSegment> segs;
		
		double pres_av;
		double pres_dev;
		
		const WellGeomProperties* getGeomProps() const;
		void printRates(const MainProperties* props) const;
		void writeRates(const MainProperties* props);
};

#endif /* WELL_HPP_ */

