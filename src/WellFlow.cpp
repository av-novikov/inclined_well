#include <cassert>
#include <new>
#include <iomanip>
#include <tinyxml.h>

#include "src/WellFlow.h"

using namespace paralution;
using std::string;
using std::stoi;
using std::stod;

WellFlow::WellFlow(const string fileName, const WellType type)
{
	loadTask(fileName, type);

	if(type == SLANTED)
		assert(fabs(props.alpha) > EQUALITY_TOLERANCE); // non-vertical
	//assert(fabs(props.alpha) < M_PI_2 - EQUALITY_TOLERANCE); // non-horizontal

	well = new Well(props.r1, props.r2, props.K, props.rw);
	well->setRate(props.rate);
	well->setUniformRate();

	props.r_obs = well->segs[int(props.K / 2)].r_bhp;

	firstTime = true;
}
WellFlow::~WellFlow()
{
	delete well;

	if (props.K > 1)
	{
		delete q;
		delete dq;
		delete b;
		delete x;
		delete a;
		delete ind_i;
		for (int i = 0; i < props.K; i++)
			delete dpdq[i];
		delete dpdq;
	}
}
void WellFlow::loadTask(const string fileName, const WellType type)
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

	TiXmlElement* xml_fluxes = xml_task->FirstChildElement("boundary_fluxes");
	props.fx1 = stod(xml_fluxes->Attribute("fx1")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
	props.fx2 = stod(xml_fluxes->Attribute("fx2")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
	props.fy1 = stod(xml_fluxes->Attribute("fy1")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
	props.fy2 = stod(xml_fluxes->Attribute("fy2")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);

	props.fx1 *= (props.visc / props.kx / props.sizes.y / props.sizes.z);
	props.fx2 *= (props.visc / props.kx / props.sizes.y / props.sizes.z);
	props.fy1 *= (props.visc / props.kx / props.sizes.x / props.sizes.z);
	props.fy2 *= (props.visc / props.kx / props.sizes.x / props.sizes.z);
	
	double alpha = props.alpha;
	switch (type)
	{
	case SLANTED:
		props.alpha = atan(tan(alpha) * sqrt(props.kz / props.kx));
		props.length *= sin(alpha) / sin(props.alpha);
		props.sizes.z *= sqrt(props.kx / props.kz);
		props.rc.z *= sqrt(props.kx / props.kz);

		props.r1 = props.r2 = props.rc;
		props.r1.x -= props.length * sin(props.alpha) / 2.0;
		props.r2.x += props.length * sin(props.alpha) / 2.0;
		props.r1.z += props.length * cos(props.alpha) / 2.0;
		props.r2.z -= props.length * cos(props.alpha) / 2.0;
		break;
	case FRAC:
		props.sizes.z *= sqrt(props.kx / props.kz);
		props.rc.z *= sqrt(props.kx / props.kz);

		props.r1 = props.r2 = props.rc;
		props.r1.x -= props.length * cos(props.alpha) / 2.0;
		props.r2.x += props.length * cos(props.alpha) / 2.0;
		props.r1.y -= props.length * sin(props.alpha) / 2.0;
		props.r2.y += props.length * sin(props.alpha) / 2.0;
		break;
	}

	if (props.K > 1)
	{
		matSize = (props.K - 1) * (props.K - 1);
		q = new double[props.K];
		dq = new double[props.K];
		b = new double[props.K - 1];
		x = new double[props.K - 1];

		a = new double[ matSize ];
		ind_i = new int[ matSize ];
		ind_j = new int[ matSize ];
		for (int i = 0; i < matSize; i++)
		{
			ind_i[i] = i / (props.K - 1);
			ind_j[i] = i % (props.K - 1);
		}

		dpdq = new double* [props.K];
		for (int i = 0; i < props.K; i++)
			dpdq[i] = new double[props.K - 1];

		//A.AllocateDENSE("A", props.K - 1, props.K - 1);
		X.Allocate("X", props.K - 1);
		B.Allocate("B", props.K - 1);
	}
}
void WellFlow::setSummator(BaseSum* _inclSum)
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
	// Fills the vector of rates
	auto fill_q = [this]() {
		for (int i = 0; i < props.K; i++)
			q[i] = well->segs[i].rate;
	};

	// Fills the vector of rate's deviations with zeros
	auto fill_dq = [this]() {
		for (int i = 0; i < props.K; i++)
			dq[i] = 0.0;
	};

	// Set rate deviation
	auto setRateDev = [this](int seg_idx, double ratio) {
		well->segs[seg_idx].rate += props.rate * ratio;
	};

	// Fills dpdq matrix
	auto fill_dpdq = [&, this](double mult) {
		double p1, p2, ratio;
		ratio = mult * 0.00001 / (double)(props.K);

		for (int i = 0; i < props.K; i++)
		{
			for (int j = 1; j < props.K; j++)
			{
				setRateDev(j, -ratio);	setRateDev(0, ratio);
				p1 = inclSum->getPres(i);

				setRateDev(j, 2.0 * ratio);	setRateDev(0, -2.0 * ratio);
				p2 = inclSum->getPres(i);

				dpdq[i][j - 1] = (p2 - p1) / (2.0 * ratio * props.rate);
			}
		}
	};

	auto solve_sys = [this]() {
		double s, p1, p2;

		for (int i = 0; i < props.K - 1; i++)
		{
			for (int j = 0; j < props.K - 1; j++)
			{
				s = 0.0;
				for (int k = 0; k < props.K - 1; k++)
					s += (dpdq[k + 1][j] - dpdq[k][j]) * (dpdq[k + 1][i] - dpdq[k][i]);

				a[i * (props.K - 1) + j] = s;

			}

			s = 0.0;
			for (int k = 0; k < props.K - 1; k++)
			{
				p1 = inclSum->getPres(k);
				p2 = inclSum->getPres(k + 1);
				s += (p2 - p1) * (dpdq[k + 1][i] - dpdq[k][i]);
			}
			B[i] = -s;
		}

		X.Zeros();
		A.Assemble(ind_i, ind_j, a, (props.K - 1) * (props.K - 1), "A", props.K - 1, props.K - 1);
		A.MoveToAccelerator();	X.MoveToAccelerator();	B.MoveToAccelerator();

		ls.SetOperator(A);
		ls.Build();
		ls.Solve(B, &X);
		//A.WriteFileMTX("mat.mtx");
		//B.WriteFileASCII("rhs.dat");
		ls.Clear();

		X.MoveToHost();

		s = 0.0;
		double tmp;
		for (int i = 0; i < props.K - 1; i++)
		{
			tmp = X[i];
			dq[i + 1] = tmp;
			s += tmp;
		}
		dq[0] = -s;
	};

	// Finds dq
	auto solve_dq = [&, this](double mult) {
		fill_dq();
		fill_dpdq(mult);
		calcPressure();
		solve_sys();
	};

	calcPressure();

	double H0 = well->pres_dev;
	if (sqrt(H0) > 0.001 * fabs(well->pres_av))
	{
		well->printRates(&props);
		fill_q();

		double mult = 0.9;
		double H = H0;

		while (H > H0 / 200.0 || (sqrt(H) > 0.002 * fabs(well->pres_av)))
		{
			solve_dq(mult);

			double tmp;
			for (int i = 0; i < props.K; i++)
			{
				tmp = q[i] + mult * dq[i];
				well->segs[i].rate = q[i] = tmp;
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

	for (int k = 0; k < props.K; k++)
	{
		WellSegment& seg = well->segs[k];
		seg.pres2D = inclSum->get2D(k);
		seg.pres3D = inclSum->get3D(k);
		seg.pres = seg.pres2D + seg.pres3D;

		well->pres_av += seg.pres;
	}

	well->pres_av /= props.K;

	// Evaluation of pressure deviation
	if (props.K > 1)
	{
		for (int i = 1; i < props.K; i++)
			well->pres_dev += (well->segs[i].pres - well->segs[i - 1].pres) *
			(well->segs[i].pres - well->segs[i - 1].pres);
	}

	well->pres_dev /= 2.0;
}
double WellFlow::getP_bhp()
{
	inclSum->prepare();

	if (firstTime)
	{
		// Body of function
		well->setUniformRate();
		firstTime = false;
		findRateDistribution();
	}

	calcPressure();

	return well->pres_av;
}