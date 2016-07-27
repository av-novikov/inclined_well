#include <cassert>
#include <new>
#include <iomanip>
#include <tinyxml.h>

#include "src/WellFlow.h"

using std::string;
using std::stoi;
using std::stod;

WellFlow::WellFlow(const string fileName)
{
	loadTask(fileName);

	assert(fabs(props.alpha) > EQUALITY_TOLERANCE); // non-vertical
	assert(fabs(props.alpha) < M_PI_2 - EQUALITY_TOLERANCE); // non-horizontal

	well = new Well(props.r1, props.r2, props.K, props.rw);
	well->setRate(props.rate);
	well->setUniformRate();

	props.r_obs = well->segs[int(props.K / 2)].r_bhp;
}

WellFlow::~WellFlow()
{
	delete well;
}

void WellFlow::loadTask(const string fileName)
{
	TiXmlDocument* xml_file = new TiXmlDocument(fileName.c_str());
	if (!xml_file->LoadFile())
		throw ("Specified taskfile is invalid");
	TiXmlElement* xml_task = xml_file->FirstChildElement("task");

	TiXmlElement* xml_dims = xml_task->FirstChildElement("dimensions");
	props.x_dim = stod(xml_dims->Attribute("x_dim"));
	props.t_dim = stod(xml_dims->Attribute("t_dim"));
	props.p_dim = stod(xml_dims->Attribute("p_dim"));

	TiXmlElement* xml_geometry = xml_task->FirstChildElement("geometry");
	TiXmlElement* xml_sizes = xml_geometry->FirstChildElement("sizes");
	props.sizes.x = stod(xml_sizes->Attribute("sx")) / props.x_dim;
	props.sizes.y = stod(xml_sizes->Attribute("sy")) / props.x_dim;
	props.sizes.z = stod(xml_sizes->Attribute("sz")) / props.x_dim;
	TiXmlElement* xml_rc = xml_geometry->FirstChildElement("rc");
	props.rc.x = stod(xml_rc->Attribute("x")) / props.x_dim;
	props.rc.y = stod(xml_rc->Attribute("y")) / props.x_dim;
	props.rc.z = stod(xml_rc->Attribute("z")) / props.x_dim;
	TiXmlElement* xml_length = xml_geometry->FirstChildElement("length");
	props.length = stod(xml_length->Attribute("value")) / props.x_dim;
	TiXmlElement* xml_alpha = xml_geometry->FirstChildElement("alpha");
	props.alpha = stod(xml_alpha->Attribute("value")) * M_PI / 180.0;
	TiXmlElement* xml_rw = xml_geometry->FirstChildElement("rw");
	props.rw = stod(xml_rw->Attribute("value")) / props.x_dim;

	TiXmlElement* xml_visc = xml_task->FirstChildElement("viscosity");
	props.visc = 1.E-3 * stod(xml_visc->Attribute("value")) / (props.p_dim * props.t_dim);

	TiXmlElement* xml_perm = xml_task->FirstChildElement("permeability");
	//props.perm = 0.986923 * 1.E-15 * stod( xml_perm->Attribute("value") ) / (props.x_dim * props.x_dim);
	props.kx = 0.986923 * 1.E-15 * stod(xml_perm->Attribute("kx")) / (props.x_dim * props.x_dim);
	props.kz = 0.986923 * 1.E-15 * stod(xml_perm->Attribute("kz")) / (props.x_dim * props.x_dim);

	TiXmlElement* xml_rate = xml_task->FirstChildElement("rate");
	props.rate = stod(xml_rate->Attribute("value")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);

	TiXmlElement* xml_idxs = xml_task->FirstChildElement("indexes");
	props.K = stoi(xml_idxs->Attribute("K"));
	props.M = stoi(xml_idxs->Attribute("M"));
	props.N = stoi(xml_idxs->Attribute("N"));
	props.L = stoi(xml_idxs->Attribute("L"));
	props.I = stoi(xml_idxs->Attribute("I"));

	TiXmlElement* xml_grid = xml_task->FirstChildElement("grid");
	props.nx = stoi(xml_grid->Attribute("nx"));
	props.ny = stoi(xml_grid->Attribute("ny"));
	props.nz = stoi(xml_grid->Attribute("nz"));

	TiXmlElement* xml_xi_c = xml_task->FirstChildElement("xi_c");
	props.xi_c = stod(xml_xi_c->Attribute("value")) / props.x_dim / props.x_dim;

	double alpha = props.alpha;
	props.alpha = atan(tan(alpha) * sqrt(props.kz / props.kx));
	props.length *= sin(alpha) / sin(props.alpha);
	props.sizes.z *= sqrt(props.kx / props.kz);
	props.rc.z *= sqrt(props.kx / props.kz);

	props.r1 = props.r2 = props.rc;
	props.r1.x -= props.length * sin(props.alpha) / 2.0;
	props.r2.x += props.length * sin(props.alpha) / 2.0;
	props.r1.z += props.length * cos(props.alpha) / 2.0;
	props.r2.z -= props.length * cos(props.alpha) / 2.0;

	/*if (props.K > 1)
	{
		q_gsl = gsl_vector_alloc(props.K);
		dq_gsl = gsl_vector_alloc(props.K);
		b_gsl = gsl_vector_alloc(props.K - 1);
		x_gsl = gsl_vector_alloc(props.K - 1);
		a_gsl = gsl_matrix_alloc(props.K - 1, props.K - 1);
		dpdq_gsl = gsl_matrix_alloc(props.K, props.K - 1);
		perm_gsl = gsl_permutation_alloc(props.K - 1);
	}*/
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

void WellFlow::findRateDistribution()
{
}

double WellFlow::getP_bhp()
{
	findRateDistribution();

	return well->pres_av;
}