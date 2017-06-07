#include "src/Well.hpp"

#include <functional>
#include <fstream>

using namespace std::placeholders;
using std::for_each;
using std::bind;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;

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
		if(type == WellType::FRAC)
			tmp3.y += props.rw / 10.0; // cos(alpha);
		else
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
	for(int i = 0; i < num; i++)
	{
		av2d += segs[i].pres2D;
		av3d += segs[i].pres3D;

		file << i << "\t" << 
			segs[i].rate * 86400.0 * mprops->x_dim * mprops->x_dim * mprops->x_dim / mprops->t_dim << "\t" <<
			segs[i].pres * mprops->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres2D * mprops->p_dim / (double)(BAR) <<
			"\t" << segs[i].pres3D * mprops->p_dim / (double)(BAR) << endl;
	}
	
	av2d /= (double)(num);	av3d /= (double)(num);
	
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
