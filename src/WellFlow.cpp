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
	assert( fabs(props.alpha) > EQUALITY_TOLERANCE); // non-vertical
	assert( fabs(props.alpha) < M_PI_2 - EQUALITY_TOLERANCE); // non-horizontal
	
	well = new Well(props.r1, props.r2, props.K, props.rw);
	well->setRate(props.rate);
	well->setUniformRate();
	
	props.r_obs = well->segs[ int(props.K / 2) ].r_bhp;
}

WellFlow::~WellFlow()
{
	delete well;
	
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
			gsl_vector_set( dq_gsl, i, 0.0 );
 	};
 	
 	// Set rate deviation
 	auto setRateDev = [this](int seg_idx, double ratio) {
		well->segs[ seg_idx ].rate += props.rate * ratio;
 	};
 	
 	// Fills dpdq matrix
 	auto fill_dpdq = [&,this](double mult) {
		double p1, p2, ratio;
		ratio = mult * 0.0001 / (double)(props.K);
		
		for(int i = 0; i < props.K; i++)
		{
			for(int j = 1; j < props.K; j++)
			{
				setRateDev(j, -ratio);	setRateDev(0, ratio);
				p1 = inclSum->getPres( i );
				
				setRateDev(j, 2.0 * ratio);	setRateDev(0, -2.0 * ratio);
				p2 = inclSum->getPres( i );
				
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
				p1 = inclSum->getPres( k );
				p2 = inclSum->getPres( k + 1 );
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
		well->printRates(&props);
		fill_q();
		
		double mult = 0.98;
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
				well->printRates(&props);

				H = well->pres_dev;
		}
	}
	else 
		well->printRates(&props);
}

void WellFlow::calcPressure()
{
	well->pres_av = well->pres_dev = 0.0;
	
	for(int k = 0; k < props.K; k++)
	{
		WellSegment& seg = well->segs[k];
		seg.pres2D = inclSum->get2D( k );
		seg.pres3D = inclSum->get3D( k );
		seg.pres = seg.pres2D + seg.pres3D;
	
		//std::cout << std::setprecision(10) << props.x_dim * seg.r_bhp << "2d = " << seg.pres2D << "\t3d = " << seg.pres3D << std::endl;
		
		well->pres_av += seg.pres;
	}
	
	well->pres_av /= well->segs.size();
	
	// Evaluation of pressure deviation
	for(int i = 1; i < props.K; i++)
		well->pres_dev += (well->segs[i].pres - well->segs[i-1].pres) * 
							(well->segs[i].pres - well->segs[i-1].pres);

	well->pres_dev /= 2.0;
}

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
	TiXmlElement* xml_rc = xml_geometry->FirstChildElement("rc");
	props.rc.x = stod( xml_rc->Attribute("x") ) / props.x_dim;
	props.rc.y = stod( xml_rc->Attribute("y") ) / props.x_dim;
	props.rc.z = stod( xml_rc->Attribute("z") ) / props.x_dim;
	TiXmlElement* xml_length = xml_geometry->FirstChildElement("length");
	props.length = stod( xml_length->Attribute("value") ) / props.x_dim;
	TiXmlElement* xml_alpha = xml_geometry->FirstChildElement("alpha");
	props.alpha = stod( xml_alpha->Attribute("value") ) * M_PI / 180.0;
	TiXmlElement* xml_rw = xml_geometry->FirstChildElement("rw");
	props.rw = stod( xml_rw->Attribute("value") ) / props.x_dim;
	
	TiXmlElement* xml_visc = xml_task->FirstChildElement("viscosity");
	props.visc = 1.E-3 * stod( xml_visc->Attribute("value") ) / (props.p_dim * props.t_dim);
	
	TiXmlElement* xml_perm = xml_task->FirstChildElement("permeability");
	//props.perm = 0.986923 * 1.E-15 * stod( xml_perm->Attribute("value") ) / (props.x_dim * props.x_dim);
	props.kx = 0.986923 * 1.E-15 * stod( xml_perm->Attribute("kx") ) / (props.x_dim * props.x_dim);
	props.kz = 0.986923 * 1.E-15 * stod( xml_perm->Attribute("kz") ) / (props.x_dim * props.x_dim);
	
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
	
	delete xml_file;
	 
	double alpha = props.alpha;
	props.alpha = atan( tan(alpha) * sqrt(props.kz / props.kx) );
	props.length *= sin(alpha) / sin(props.alpha);
	props.sizes.z *= sqrt(props.kx / props.kz);
	props.rc.z *= sqrt(props.kx / props.kz);
	
	props.r1 = props.r2 = props.rc;
	props.r1.x -= props.length * sin(props.alpha) / 2.0;
	props.r2.x += props.length * sin(props.alpha) / 2.0;
	props.r1.z += props.length * cos(props.alpha) / 2.0;
	props.r2.z -= props.length * cos(props.alpha) / 2.0;
	
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

double WellFlow::getP_bhp()
{
	findRateDistribution();
	
	return well->pres_av;
}
