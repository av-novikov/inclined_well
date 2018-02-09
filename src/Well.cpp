#define _USE_MATH_DEFINES
#include "src/Well.hpp"

#include <functional>
#include <fstream>
#include <cmath>

using namespace std::placeholders;
using std::for_each;
using std::bind;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::erf;

const double MainProperties::porosity = 0.1;
const double MainProperties::compressibility = 5.E-5;

Well::Well(const WellGeomProperties& _props, const WellType& _type, const string _name, const int _well_idx) : props(_props), num(_props.seg_num), type(_type), name(_name), well_idx(_well_idx)
{
	pres_av = pres_dev = rate = 0.0;
	
	Point tmp1 = props.r1;
	Point tmp2, tmp3;
	double tau1, tau2;
	for(int i = 0; i < num; i++)
	{
		tmp2 = props.r1 + (double)( i + 1 ) * (props.r2 - props.r1) / (double)( num );
		tmp3 = (tmp1 + tmp2) / 2.0;		
		tmp3.y += props.rw; // cos(alpha);
		tau1 = (double)(i) / (double)(num);
		tau2 = (double)(i+1) / (double)(num);
		segs.push_back( WellSegment(tmp1, tmp2, tmp3, tau1, tau2, well_idx, i) );
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
		seg.rate = seg.length / props.length * rate;
	});
}
void Well::setGaussRate(const double sigma)
{
	const double A = rate / erf(props.length / sqrt(8.0) / sigma);
	const double x_c = (segs[0].r1.x + segs[segs.size() - 1].r2.x ) / 2.0;

	double sum = 0.0;
	for (size_t i = 0; i < num; i++)
	{
		auto& seg = segs[i];
		if (i != num / 2 - 1)
			seg.rate = A / sigma / sqrt(2.0 * M_PI) * exp(-(seg.r_bhp.x - x_c) * (seg.r_bhp.x - x_c) / 2.0 / sigma / sigma) * seg.length;
		else
			seg.rate = rate / 2.0 - sum;
		sum += seg.rate;
	}
}
void Well::setParabolicRate(const double ratio)
{
	const double xf = props.length / 2.0;
	const double a = 6.0 * rate / xf / xf / xf / (1 + ratio);
	const double b = ratio / (1.0 + ratio) * rate / xf / 2.0;
	const double xc = (segs[0].r1.x + segs[segs.size() - 1].r2.x) / 2.0;

	double sum = 0.0;
	for (int i = 0; i < num / 2; i++)
	{
		auto& seg1 = segs[num / 2 + i];
		auto& seg2 = segs[num / 2 - i - 1];
		if (i != num / 2 - 1)
			seg1.rate = seg2.rate = (a * (seg1.r_bhp.x - xc - xf / 2.0) * (seg1.r_bhp.x - xc - xf / 2.0) + b) * seg1.length;
		else
			seg1.rate = seg2.rate = rate / 2.0 - sum;
		sum += seg1.rate;			
	}
}
void Well::setParabolicAllRate(const double ratio)
{
	const double xf = props.length / 2.0;
	const double a = 3.0 * rate / xf / xf / xf / (1 + ratio) / 2.0;
	const double b = ratio / (1.0 + ratio) * rate / xf / 2.0;
	const double xc = (segs[0].r1.x + segs[segs.size() - 1].r2.x) / 2.0;
	
	double sum = 0.0;
	for (int i = 0; i < num / 2; i++)
	{
		auto& seg1 = segs[num / 2 + i];
		auto& seg2 = segs[num / 2 - i - 1];
		if (i != num / 2 - 1)
			seg1.rate = seg2.rate = (a * (seg1.r_bhp.x - xc) * (seg1.r_bhp.x - xc) + b) * seg1.length;
		else
			seg1.rate = seg2.rate = rate / 2.0 - sum;
		sum += seg1.rate;
	}
}
void Well::printRates(const MainProperties* mprops) const
{
	double av2d = 0.0, av3d = 0.0;
	cout << name << endl;
	double sum_rate = 0.0;
	for(int i = 0; i < num; i++)
	{
		av2d += segs[i].pres2D;
		av3d += segs[i].pres3D;
		cout << "--- " << i << " ---\tRate = " << segs[i].rate * 86400.0 * mprops->x_dim * mprops->x_dim * mprops->x_dim / mprops->t_dim << 
									"\tPressure = " << segs[i].pres * mprops->p_dim / (double)(BAR) << "\t" << segs[i].pres2D * mprops->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres3D * mprops->p_dim / (double)(BAR) << endl;
		sum_rate += segs[i].rate * 86400.0 * mprops->x_dim * mprops->x_dim * mprops->x_dim / mprops->t_dim;
	}
	av2d /= (double)(num);	av3d /= (double)(num);
	
	cout << "Sum_rate = " << sum_rate << endl;
	cout << "Av. pressure = " << pres_av * mprops->p_dim / BAR << endl;
	cout << "Av. 2D = " << av2d * mprops->p_dim / BAR << endl;
	cout << "Av. 3D = " << av3d * mprops->p_dim / BAR << endl;
	cout << "Deviation = " << pres_dev * mprops->p_dim * mprops->p_dim / BAR / BAR << endl;
}
void Well::writeRates(const MainProperties* mprops)
{
	double av2d = 0.0, av3d = 0.0;	
	ofstream file;
	file.open("rate.dat", ofstream::out);
	double sum_rate = 0.0;
	for(int i = 0; i < num; i++)
	{
		av2d += segs[i].pres2D;
		av3d += segs[i].pres3D;

		file << i << "\t" << 
			segs[i].rate * 86400.0 * mprops->x_dim * mprops->x_dim * mprops->x_dim / mprops->t_dim << "\t" <<
			segs[i].pres * mprops->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres2D * mprops->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres3D * mprops->p_dim / (double)(BAR) << endl;
		sum_rate += segs[i].rate * 86400.0 * mprops->x_dim * mprops->x_dim * mprops->x_dim / mprops->t_dim;
	}
	
	av2d /= (double)(num);	av3d /= (double)(num);
	
	file << "Sum_rate = " << sum_rate << endl;
	file << "Av. pressure = " << pres_av * mprops->p_dim / BAR << endl;
	file << "Av. 2D = " << av2d * mprops->p_dim / BAR << endl;
	file << "Av. 3D = " << av3d * mprops->p_dim / BAR << endl;
	file << "Deviation = " << pres_dev * mprops->p_dim * mprops->p_dim / BAR / BAR << endl;
	
	file.close();
}
const WellGeomProperties* Well::getGeomProps() const
{
	return &props;
}
