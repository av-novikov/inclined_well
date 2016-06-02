#include "src/Well.hpp"

#include <functional>

using namespace std::placeholders;
using std::for_each;
using std::bind;

Well::Well(const Point& _r1, const Point& _r2, const int _num) : r1(_r1), r2(_r2), num(_num)
{
	length = sqrt((r2 - r1) * (r2 - r1));
	pres_av = rate = 0.0;
	
	Point tmp1 = r1;
	Point tmp2;
	for(int i = 0; i < num; i++)
	{
		tmp2 = r1 + (double)( i + 1 ) * (r2 - r1) / (double)( num );
		segs.push_back( WellSegment(tmp1, tmp2) );
		tmp1 = tmp2;
	};
}

Well::~Well()
{
}

void Well::setRate(double _rate)
{
	rate = _rate;
}

void Well::setUniformRate()
{
	for_each(segs.begin(), segs.end(), bind(&Well::uniformRate, this, _1) );
}

void Well::uniformRate(WellSegment& seg)
{
	seg.rate = seg.length / length * rate;
}
