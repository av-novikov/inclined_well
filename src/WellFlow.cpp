#include <cassert>
#include <new>
#include <iomanip>
#include <tinyxml.h>

#include "src/WellFlow.hpp"

using std::string;
using std::for_each;
using std::stoi;
using std::stod;

WellFlow::WellFlow(const string fileName)
{
	loadTask(fileName);
	assert( fabs(props.alpha) > EQUALITY_TOLERANCE);
	
	well = new Well(props.r1, props.r2, props.K, props.rw);
	well->setRate(props.rate);
	well->setUniformRate();
	
	inclSum = new InclinedSum(&props, well);
}

WellFlow::~WellFlow()
{
	//delete well;
	//delete inclSum;
	
	if( props.K > 1 )
	{
		gsl_vector_free( q_gsl );
		gsl_vector_free( dq_gsl );
		gsl_vector_free( b_gsl );
		gsl_vector_free( x_gsl );
		gsl_matrix_free( a_gsl );
		gsl_matrix_free( dpdq_gsl );
		gsl_permutation_free( perm_gsl );	
	}
}

void WellFlow::findRateDistribution()
{
	// Fills the vector of rates
	auto fill_q = [this]() {
		for(int i = 0; i < props.K; i++)
			gsl_vector_set( q_gsl, i, well->segs[i].rate );
 	};
 	
 	// Fills the vector of rate's deviations with zeros
	auto fill_dq = [this]() {
		for(int i = 0; i < props.K; i++)
			gsl_vector_set( q_gsl, i, 0.0 );
 	};
 	
 	// Set rate deviation
 	auto setRateDev = [this](int seg_idx, double ratio) {
		well->segs[ seg_idx ].rate += props.rate * ratio;
 	};
 	
 	// Fills dpdq matrix
 	auto fill_dpdq = [&,this](double mult) {
		double p1, p2, ratio;
		ratio = mult * 0.001 / (double)(props.K);
		
		for(int i = 0; i < props.K; i++)
		{
			for(int j = 1; j < props.K; j++)
			{
				setRateDev(j, -ratio);	setRateDev(0, ratio);
				p1 = inclSum->getPres( well->segs[i].r_bhp );
				
				setRateDev(j, 2.0 * ratio);	setRateDev(0, -2.0 * ratio);
				p2 = inclSum->getPres( well->segs[i].r_bhp );
				
				gsl_matrix_set(dpdq_gsl, i, j-1, (p2 - p1) / (2.0 * ratio * props.rate) );
			}
		}
	};
	
	auto solve_sys = [this]() {
		double s, p1, p2;
		
		for(int i = 0; i < props.K-1; i++)
		{
			for(int j = 0; j < props.K-1; j++)	
			{
				s = 0.0;
				for(int k = 0; k < props.K-1; k++)
					s += ( gsl_matrix_get(dpdq_gsl, k+1, j) - gsl_matrix_get(dpdq_gsl, k, j) ) * 
						( gsl_matrix_get(dpdq_gsl, k+1, i) - gsl_matrix_get(dpdq_gsl, k, i) );
				
				gsl_matrix_set(a_gsl, i, j, s);
			}
			
			s = 0.0;
			for(int k = 0; k < props.K-1; k++)
			{
				p1 = inclSum->getPres( well->segs[k].r_bhp );
				p2 = inclSum->getPres( well->segs[k+1].r_bhp );
				s += ( p2 - p1 ) * ( gsl_matrix_get(dpdq_gsl, k+1, i) - gsl_matrix_get(dpdq_gsl, k, i) );
			}
			gsl_vector_set(b_gsl, i, -s);
		}
		
		int ss;
		gsl_linalg_LU_decomp(a_gsl, perm_gsl, &ss);
		gsl_linalg_LU_solve (a_gsl, perm_gsl, b_gsl, x_gsl);
		
		s = 0.0;
		double tmp;
		for(int i = 0; i < props.K-1; i++)
		{
			tmp = gsl_vector_get(x_gsl, i);
			gsl_vector_set(dq_gsl, i+1, tmp);
			s += tmp;
		}
		gsl_vector_set(dq_gsl, 0, -s);
	};	
	
	// Finds dq
 	auto solve_dq = [&,this](double mult) {
		fill_dq();
		fill_dpdq(mult);
		calcPressure();
		solve_sys();		
 	};
	
	// Body of function
	well->setUniformRate();
	calcPressure();
	
	double H0 = well->pres_dev;
	if(H0 > 0.1)
	{
		well->printRates();
		fill_q();
		
		double mult = 0.9;
		double H = H0;	
		
		while(H > H0 / 50.0 || H > 0.05)
		{
				solve_dq(mult);
			
				double tmp;
				for(int i = 0; i < props.K; i++)
				{
					tmp = gsl_vector_get(q_gsl, i) + mult * gsl_vector_get(dq_gsl, i);
					gsl_vector_set(q_gsl, i, tmp);
					
					well->segs[i].rate = tmp;
				}
				
				calcPressure();
				well->printRates();

				H = well->pres_dev;
		}
	}
}

void WellFlow::calcPressure()
{
	well->pres_av = well->pres_dev = 0.0;
	
	for_each(well->segs.begin(), well->segs.end(), [this](WellSegment& seg)
	{
		seg.pres = inclSum->getPres( seg.r_bhp );
		
		well->pres_av += seg.pres;
	});
	well->pres_av /= well->segs.size();
	
	// Evaluation of pressure deviation
	for(int i = 1; i < props.K; i++)
		well->pres_dev += (well->segs[i].pres - well->segs[i-1].pres) * 
							(well->segs[i].pres - well->segs[i-1].pres);

	well->pres_dev /= 2.0;
}

/*double WellFlow::highIntegral2D(const Point& r, double xi_c)
{
	double sum = 0.0;
	double buf;
	
	for(int k = 0; k < props.K; k++)
		for(int m = 1; m <= props.M; m++)
			for(int n = 1; n <= props.N; n++)
			{
				buf = ((double)(m) * (double)(m) / props.sizes.x / props.sizes.x + (double)(n) * (double)(n) / props.sizes.y / props.sizes.y);
				
				sum += well->segs[k].rate / well->segs[k].length / (double)(m) / buf * 
					sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x)) * 
					sin(M_PI * (double)(n) * r.y / props.sizes.y) * sin(M_PI * (double)(n) * props.r1.y / props.sizes.y) * 
					exp(-buf * M_PI * M_PI * xi_c);
			}
	sum *= (4.0 * props.visc / M_PI / M_PI / M_PI / props.sizes.y / props.sizes.z / props.perm / sin(props.alpha));
	
	return sum;
}

double WellFlow::lowIntegral2D(const Point& r, double xi_c)
{
	double sum = 0.0;
	double Iminus, Iplus;
	double tmp1, tmp2;
	
	for(int k = 0; k < props.K; k++)
		for(int m = 1; m <= props.M; m++)
			for(int i = -props.I; i <= props.I; i++)
			{
				Iminus = exp(-M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
							(1.0 - erf((fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) - 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) + 
							exp(2.0 * M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
							(erf((fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) + 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) - 1.0));
				
				Iplus = exp(-M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
									(1.0 - erf((fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) - 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) + 
							exp(2.0 * M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) * 
							(erf((fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) + 2.0 * M_PI * (double)(m) / props.sizes.x * xi_c) / 2.0 / sqrt(xi_c)) - 1.0));
				
				
				tmp1 = (fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / 2.0 - M_PI * (double)(m) / props.sizes.x * xi_c) / sqrt(xi_c);
				tmp2 = (fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / 2.0 + M_PI * (double)(m) / props.sizes.x * xi_c) / sqrt(xi_c);
				
				Iminus = (2.0 * exp( -M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x)- 
						 exp(-tmp1*tmp1 - M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x)) / tmp1 - 
						 exp(-tmp2*tmp2 + M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) / tmp2;
				
				tmp1 = (fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / 2.0 - M_PI * (double)(m) / props.sizes.x * xi_c) / sqrt(xi_c);
				tmp2 = (fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / 2.0 + M_PI * (double)(m) / props.sizes.x * xi_c) / sqrt(xi_c);
				
				Iplus = (2.0 * exp( -M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x)-
						exp(-tmp1*tmp1 - M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x)) / tmp1 - 
						exp(-tmp2*tmp2 + M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) / tmp2;
				
				sum += well->segs[k].rate / well->segs[k].length / (double)(m) / (double)(m) * 
					sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x)) * 
					(Iminus - Iplus);
					
				//std::cout << "Iminus = " << Iminus << "\t" << "Iplus = " << Iplus << std::endl;
				//std::cout << "sum = " << sum << std::endl;
			}
	
	sum *= (props.visc * props.sizes.x / 2.0 / M_PI / M_PI / props.sizes.z / props.perm / sin(props.alpha));
	
	return sum;	
}

double WellFlow::calc3D(const Point& r)
{
	double sum = 0.0;
	double F, buf;
	
	for(int k = 0; k < props.K; k++)
		for(int m = 1; m <= props.M; m++)
			for(int l = 1; l <= props.L; l++)
			{
				buf = sqrt((double)(m) * (double)(m) / props.sizes.x / props.sizes.x + (double)(l) * (double)(l) / props.sizes.z / props.sizes.z);
				
				for(int i = -props.I; i <= props.I; i++)
				{
					F = ((cos(M_PI * (double)(m) * (well->segs[k].r1.x + 
							tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x + 
							M_PI * (double)(l) * well->segs[k].r1.z / props.sizes.z) -
						cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x +
							M_PI * (double)(l) * well->segs[k].r2.z / props.sizes.z)) /
							(M_PI * (double)(l) / props.sizes.z - M_PI * (double)(m) * tan(props.alpha) / props.sizes.x) - 
						(cos(M_PI * (double)(m) * (well->segs[k].r1.x + 
							tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x - 
							M_PI * (double)(l) * well->segs[k].r1.z / props.sizes.z) -
						cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x -
							M_PI * (double)(l) * well->segs[k].r2.z / props.sizes.z)) /
							(M_PI * (double)(l) / props.sizes.z + M_PI * (double)(m) * tan(props.alpha) / props.sizes.x)
							) / 2.0;
						
						sum += well->segs[k].rate / well->segs[k].length / buf * 
							( exp(-M_PI * buf * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y)) - 
							exp(-M_PI * buf * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y)) ) *
							sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
							cos(M_PI * (double)(l) * r.z / props.sizes.z);
				}
			}
				
	sum *= (2.0 * props.visc / M_PI / props.sizes.x / props.sizes.z / props.perm / cos(props.alpha));
				
	return sum;
}

double WellFlow::calcPressure(const Point& r, bool debug)
{	
	double sum = 0.0;

	if(debug)
	{
		for(int k = 0; k < props.K; k++)
			for(int m = 1; m <= props.M; m++)
				for(int i = -props.I; i <= props.I; i++)
				{
					sum += well->segs[k].rate / well->segs[k].length / (double)(m) / (double)(m) * 
						( exp(-M_PI * (double)(m) * fabs(r.y - props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) - 
						exp(-M_PI * (double)(m) * fabs(r.y + props.r1.y + 2.0 * (double)(i) * props.sizes.y) / props.sizes.x) ) *
						sin(M_PI * (double)(m) * r.x / props.sizes.x) * 
						(cos(M_PI * (double)(m) * well->segs[k].r1.x / props.sizes.x) -
						cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props.alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props.sizes.x));
				}
		sum *= (props.visc * props.sizes.x / M_PI / M_PI / props.sizes.z / props.perm / sin(props.alpha));
		
		sum += calc3D(r);
		
	} else
	{
		sum = lowIntegral2D(r, props.xi_c) + highIntegral2D(r, props.xi_c);
	}
		
	return sum * props.p_dim;
}*/

void WellFlow::loadTask(const string fileName)
{
	TiXmlDocument* xml_file = new TiXmlDocument( fileName.c_str() );
	if (!xml_file->LoadFile())
		throw ("Specified taskfile is invalid");
	TiXmlElement* xml_task = xml_file->FirstChildElement("task");

	TiXmlElement* xml_dims = xml_task->FirstChildElement("dimensions");
	props.x_dim = stod( xml_dims->Attribute("x_dim") );
	props.t_dim = stod( xml_dims->Attribute("t_dim") );
	props.p_dim = stod( xml_dims->Attribute("p_dim") );
	
	TiXmlElement* xml_geometry = xml_task->FirstChildElement("geometry");
	TiXmlElement* xml_sizes = xml_geometry->FirstChildElement("sizes");
	props.sizes.x = stod( xml_sizes->Attribute("sx") ) / props.x_dim;
	props.sizes.y = stod( xml_sizes->Attribute("sy") ) / props.x_dim;
	props.sizes.z = stod( xml_sizes->Attribute("sz") ) / props.x_dim;
	TiXmlElement* xml_r1 = xml_geometry->FirstChildElement("r1");
	props.r1.x = stod( xml_r1->Attribute("x") ) / props.x_dim;
	props.r1.y = stod( xml_r1->Attribute("y") ) / props.x_dim;
	props.r1.z = stod( xml_r1->Attribute("z") ) / props.x_dim;
	TiXmlElement* xml_alpha = xml_geometry->FirstChildElement("alpha");
	props.alpha = stod( xml_alpha->Attribute("value") ) * M_PI / 180.0;
	TiXmlElement* xml_rw = xml_geometry->FirstChildElement("rw");
	props.rw = stod( xml_rw->Attribute("value") ) / props.x_dim;
	
	TiXmlElement* xml_visc = xml_task->FirstChildElement("viscosity");
	props.visc = 1.E-3 * stod( xml_visc->Attribute("value") ) / (props.p_dim * props.t_dim);
	
	TiXmlElement* xml_perm = xml_task->FirstChildElement("permeability");
	props.perm = 0.986923 * 1.E-15 * stod( xml_perm->Attribute("value") ) / (props.x_dim * props.x_dim);
	
	TiXmlElement* xml_rate = xml_task->FirstChildElement("rate");
	props.rate = stod( xml_rate->Attribute("value") ) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
	
	TiXmlElement* xml_idxs = xml_task->FirstChildElement("indexes");
	props.K = stoi( xml_idxs->Attribute("K") );
	props.M = stoi( xml_idxs->Attribute("M") );
	props.N = stoi( xml_idxs->Attribute("N") );
	props.L = stoi( xml_idxs->Attribute("L") );
	props.I = stoi( xml_idxs->Attribute("I") );

	TiXmlElement* xml_grid = xml_task->FirstChildElement("grid");
	props.nx = stoi( xml_grid->Attribute("nx") );
	props.ny = stoi( xml_grid->Attribute("ny") );
	props.nz = stoi( xml_grid->Attribute("nz") );
	
	TiXmlElement* xml_xi_c = xml_task->FirstChildElement("xi_c");
	props.xi_c = stod( xml_xi_c->Attribute("value") ) / props.x_dim / props.x_dim;
	
	props.r2 = props.r1;
	props.r2.x += tan(props.alpha) * props.sizes.z;
	props.r2.z = -props.sizes.z;
	
	if( props.K > 1 )
	{
		q_gsl = gsl_vector_alloc( props.K );
		dq_gsl = gsl_vector_alloc( props.K );
		b_gsl = gsl_vector_alloc( props.K - 1 );
		x_gsl = gsl_vector_alloc( props.K - 1 );
		a_gsl = gsl_matrix_alloc( props.K - 1, props.K - 1 );
		dpdq_gsl = gsl_matrix_alloc( props.K, props.K - 1 );
		perm_gsl = gsl_permutation_alloc( props.K - 1 );	
	}
}

void WellFlow::setSummator(InclinedSum* _inclSum)
{
	inclSum = _inclSum;
}

const Parameters* WellFlow::getProps() const
{
	return &props;
}

const Well* WellFlow::getWell() const
{
	return well;
}

Point WellFlow::getObsPoint() const
{
	return well->segs[ int(props.K / 2) ].r_bhp;
}

double WellFlow::getP_bhp()
{
	findRateDistribution();
	
	Point p = getObsPoint();
	std::cout << std::setprecision(10);
	std::cout << "Observation point: "<< props.x_dim * p;
	return inclSum->getPres(p);
}
