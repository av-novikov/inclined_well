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
	double av2d = 0.0, av3d = 0.0;
	for(int i = 0; i < num; i++)
	{
		av2d += segs[i].pres2D;
		av3d += segs[i].pres3D;
		cout << "--- " << i << " ---\tRate = " << segs[i].rate * 86400.0 * props->x_dim * props->x_dim * props->x_dim / props->t_dim << 
									"\tPressure = " << segs[i].pres * props->p_dim / (double)(BAR) << "\t" << segs[i].pres2D * props->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres3D * props->p_dim / (double)(BAR) << endl;
	}
	av2d /= (double)(num);	av3d /= (double)(num);
	
	cout << "Av. pressure = " << pres_av * props->p_dim / BAR << endl;
	cout << "Av. 2D = " << av2d * props->p_dim / BAR << endl;
	cout << "Av. 3D = " << av3d * props->p_dim / BAR << endl;
	cout << "Deviation = " << pres_dev * props->p_dim * props->p_dim / BAR / BAR << endl;
}

void Well::writeRates(const Parameters* props)
{
	double av2d = 0.0, av3d = 0.0;	
	ofstream file;
	file.open("rate.dat", ofstream::out);
	for(int i = 0; i < num; i++)
	{
		av2d += segs[i].pres2D;
		av3d += segs[i].pres3D;

		file << i << "\t" << 
			segs[i].rate * 86400.0 * props->x_dim * props->x_dim * props->x_dim / props->t_dim << "\t" <<
			segs[i].pres * props->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres2D * props->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres3D * props->p_dim / (double)(BAR) << endl;
	}
	
	av2d /= (double)(num);	av3d /= (double)(num);
	
	file << "Av. pressure = " << pres_av * props->p_dim / BAR << endl;
	file << "Av. 2D = " << av2d * props->p_dim / BAR << endl;
	file << "Av. 3D = " << av3d * props->p_dim / BAR << endl;
	file << "Deviation = " << pres_dev * props->p_dim * props->p_dim / BAR / BAR << endl;
	
	file.close();
}
