#include "src/Well.hpp"

#include <functional>
#include <fstream>

using namespace std::placeholders;
using std::for_each;
using std::bind;
using std::cout;
using std::endl;
using std::ofstream;

Well::Well(const Point& _r1, const Point& _r2, const int _num, const double _r_w) : r1(_r1), r2(_r2), num(_num), r_w(_r_w)
{
	length = sqrt((r2 - r1) * (r2 - r1));
	alpha = atan( (r2.x - r1.x) / (r1.z - r2.z) );
	pres_av = pres_dev = rate = 0.0;
	
	Point tmp1 = r1;
	Point tmp2, tmp3;
	for(int i = 0; i < num; i++)
	{
		tmp2 = r1 + (double)( i + 1 ) * (r2 - r1) / (double)( num );
		tmp3 = (tmp1 + tmp2) / 2.0;		tmp3.y += r_w; // cos(alpha);
		segs.push_back( WellSegment(tmp1, tmp2, tmp3) );
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
	for_each(segs.begin(), segs.end(), [this](WellSegment& seg)
	{
		seg.rate = seg.length / length * rate;
	});
}

void Well::printRates(const Parameters* props)
{
	for(int i = 0; i < num; i++)
	{
		cout << "--- " << i << " ---\tRate = " << segs[i].rate * 86400.0 * props->x_dim * props->x_dim * props->x_dim / props->t_dim << 
									"\tPressure = " << segs[i].pres * props->p_dim / (double)(BAR) << endl;
	}
	cout << "Av. pressure = " << pres_av * props->p_dim / BAR << "\tDeviation = " << pres_dev * props->p_dim * props->p_dim / BAR / BAR << endl;
	
	
	ofstream file;
	file.open("rate.dat", ofstream::out);
	for(int i = 0; i < num; i++)
	{
		file << i << "\t" << 
			segs[i].rate * 86400.0 * props->x_dim * props->x_dim * props->x_dim / props->t_dim << "\t" <<
			segs[i].rate * 86400.0 * props->x_dim * props->x_dim / props->t_dim / segs[i].length << "\t" <<
			segs[i].pres * props->p_dim / (double)(BAR) << endl;
	}
	file << "Av. pressure = " << pres_av * props->p_dim / BAR << "\tDeviation = " << pres_dev * props->p_dim * props->p_dim / BAR / BAR << endl;
	
	file.close();
}
