#include <new>
#include <iomanip>
#include <new>
#include <cassert>
#include "tinyxml2.h"

#include "src/WellFlow.h"
#include "src/inclined_sum/VerticalDirichlet.h"
#include "src/inclined_sum/neumann/VerticalNeumannUniformBoundaries.h"
#include "src/inclined_sum/neumann/VerticalNeumannGaussBoundaries.h"
#include "src/inclined_sum/Frac2dSum.hpp"
#include "src/inclined_sum/InclinedSum.hpp"
#include "src/inclined_sum/Inclined3dSum.h"
#include "src/inclined_sum/welltests/HorizontalLogDerivation.hpp"
#include <array>

#define MAX_WELL_SIZE 10

using namespace paralution;
using namespace tinyxml2;
using std::string;
using std::stoi;
using std::stod;
using std::array;
using std::vector;

#define FRAC_COEF 0.999

WellFlow::WellFlow(const string fileName)
{
	perf_length = perf_surf = 0.0;
	load(fileName);

	bool isFracAmongWells = false;
	for (const auto& well : wells)
		if (well.type == WellType::FRAC)
		{
			isFracAmongWells = true;
			break;
		}
	isFracAmongWells *= (wells.size() > 1);
	for (auto& well : wells)
	{
		if (isFracAmongWells)
		{
			if (well.type == WellType::FRAC)
				well.setRate(FRAC_COEF * props.rate);
			else
				well.setRate((1.0 - FRAC_COEF) * props.rate);
		}
		else
			well.setRate(props.rate * well.getGeomProps()->length / perf_length);
		well.setUniformRate();
		for (int i = 0; i < well.segs.size(); i++)
			segs.push_back(&well.segs[i]);
	}

	for (auto& sum : summators)
		sum->setSegments(&segs);
	if(insideSum != nullptr)
		insideSum->setSegments(&segs);

	//props.r_obs = well->segs[int(props.K / 2)].r_bhp;

	firstTime = true;
}
WellFlow::~WellFlow()
{
	wells.clear();
	fractures.clear();
	for (auto& summator : summators)
		delete summator;
	summators.clear();
	segs.clear();

	if (seg_num > 1)
	{
		delete q;
		delete dq;
		delete b;
		delete x;
		delete a;
		delete ind_i;
		for (int i = 0; i < seg_num; i++)
			delete dpdq[i];
		delete dpdq;
	}
}
void WellFlow::load(const string fileName)
{
	XMLDocument xml_file;
	XMLError eResult = xml_file.LoadFile(fileName.c_str());
	if (eResult != tinyxml2::XML_SUCCESS) exit(-1);

	XMLNode* xml_task = xml_file.FirstChildElement("task");
	XMLElement* xml_main = xml_task->FirstChildElement("main");
	XMLElement* xml_dims = xml_main->FirstChildElement("dimensions");
	props.x_dim = stod(xml_dims->Attribute("x_dim"));
	props.t_dim = stod(xml_dims->Attribute("t_dim"));
	props.p_dim = stod(xml_dims->Attribute("p_dim"));
	XMLElement* xml_visc = xml_main->FirstChildElement("viscosity");
	props.visc = 1.E-3 * stod(xml_visc->Attribute("value")) / (props.p_dim * props.t_dim);
	XMLElement* xml_perm = xml_main->FirstChildElement("permeability");
	//props.perm = 0.986923 * 1.E-15 * stod( xml_perm->Attribute("value") ) / (props.x_dim * props.x_dim);
	props.kx = 0.986923 * 1.E-15 * stod(xml_perm->Attribute("kx")) / (props.x_dim * props.x_dim);
	props.kz = 0.986923 * 1.E-15 * stod(xml_perm->Attribute("kz")) / (props.x_dim * props.x_dim);
	XMLElement* xml_rate = xml_main->FirstChildElement("rate");
	props.rate = stod(xml_rate->Attribute("value")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
	XMLElement* xml_sizes = xml_main->FirstChildElement("sizes");
	props.sizes.x = stod(xml_sizes->Attribute("sx")) / props.x_dim;
	props.sizes.y = stod(xml_sizes->Attribute("sy")) / props.x_dim;
	props.sizes.z = stod(xml_sizes->Attribute("sz")) / props.x_dim;

	seg_num = 0;
	int seg_idx = 0;
	wells.reserve(MAX_WELL_SIZE);
	for (XMLElement* xml_well = xml_task->FirstChildElement("well"); 
			xml_well != NULL; 
			xml_well = xml_well->NextSiblingElement("well"))
	{
		WellGeomProperties geom_props;
		string w_name = xml_well->Attribute("name");
		string w_type = xml_well->Attribute("type");
		geom_props.seg_num = stoi(xml_well->Attribute("K"));
		seg_num += geom_props.seg_num;
		if(w_type == "slanted")
			loadSlantedGeomProps(xml_well, &geom_props);
		else
			loadGeomProps(xml_well, &geom_props);

		SummatorProperties sum_props;
		sum_props.K = geom_props.seg_num;
		loadSumProps(xml_well->FirstChildElement("sum_props"), &sum_props);

		double alpha = geom_props.alpha;
		if(w_type == "vertical")
		{
			assert(geom_props.alpha == 0.0);

			geom_props.r1 = geom_props.r2 = geom_props.rc;
			geom_props.r1.x -= geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r2.x += geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r1.z += geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r2.z -= geom_props.length * cos(geom_props.alpha) / 2.0;

			perf_length += geom_props.length;
			wells.push_back(Well(geom_props, WellType::VERTICAL, w_name, seg_idx++));
			summators.push_back(static_cast<BaseSum*>(new VerticalDirichlet(sum_props, &props, &wells[wells.size()-1])));
		}
		else if (w_type == "vertical_neumann")
		{
			geom_props.alpha = atan(tan(alpha) * sqrt(props.kz / props.kx));
			geom_props.length *= sin(alpha) / sin(geom_props.alpha);
			props.sizes.z *= sqrt(props.kx / props.kz);
			geom_props.rc.z *= sqrt(props.kx / props.kz);

			geom_props.r1 = geom_props.r2 = geom_props.rc;
			geom_props.r1.x -= geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r2.x += geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r1.z += geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r2.z -= geom_props.length * cos(geom_props.alpha) / 2.0;

			XMLElement* xml_fluxes = xml_main->FirstChildElement("boundary_fluxes");
			props.fx1 = stod(xml_fluxes->Attribute("fx1")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
			props.fx2 = stod(xml_fluxes->Attribute("fx2")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
			props.fy1 = stod(xml_fluxes->Attribute("fy1")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);
			props.fy2 = stod(xml_fluxes->Attribute("fy2")) / 86400.0 / (props.x_dim * props.x_dim * props.x_dim / props.t_dim);

			props.fx1 *= (props.visc / props.kx / props.sizes.y / props.sizes.z);
			props.fx2 *= (props.visc / props.kx / props.sizes.y / props.sizes.z);
			props.fy1 *= (props.visc / props.kx / props.sizes.x / props.sizes.z);
			props.fy2 *= (props.visc / props.kx / props.sizes.x / props.sizes.z);

			perf_length += geom_props.length;
			wells.push_back(Well(geom_props, WellType::VERTICAL_NEUMANN, w_name, seg_idx++));
			summators.push_back(static_cast<BaseSum*>(new VerticalNeumannUniformBoundaries(sum_props, &props, &wells[wells.size() - 1])));
		}
		else if (w_type == "slanted")
		{
			perf_length += geom_props.length;
			wells.push_back(Well(geom_props, WellType::SLANTED, w_name, seg_idx++));
			summators.push_back(static_cast<BaseSum*>(new Inclined3dSum(sum_props, &props, &wells[wells.size() - 1])));
		}
		else if (w_type == "slanted_old")
		{
			geom_props.alpha = atan(tan(alpha) * sqrt(props.kz / props.kx));
			geom_props.length *= sin(alpha) / sin(geom_props.alpha);
			props.sizes.z *= sqrt(props.kx / props.kz);
			geom_props.rc.z *= sqrt(props.kx / props.kz);

			geom_props.r1 = geom_props.r2 = geom_props.rc;
			geom_props.r1.x -= geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r2.x += geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r1.z += geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r2.z -= geom_props.length * cos(geom_props.alpha) / 2.0;

			perf_length += geom_props.length;
			wells.push_back(Well(geom_props, WellType::SLANTED, w_name, seg_idx++));
			summators.push_back(static_cast<BaseSum*>(new InclinedSum(sum_props, &props, &wells[wells.size() - 1])));
		}
		else if (w_type == "horizontal")
		{
			geom_props.alpha = atan(tan(alpha) * sqrt(props.kz / props.kx));
			geom_props.length *= sin(alpha) / sin(geom_props.alpha);
			props.sizes.z *= sqrt(props.kx / props.kz);
			geom_props.rc.z *= sqrt(props.kx / props.kz);

			geom_props.r1 = geom_props.r2 = geom_props.rc;
			geom_props.r1.x -= geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r2.x += geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r1.z += geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r2.z -= geom_props.length * cos(geom_props.alpha) / 2.0;

			perf_length += geom_props.length;
			wells.push_back(Well(geom_props, WellType::HORIZONTAL, w_name, seg_idx++));
			summators.push_back(static_cast<BaseSum*>(new HorizontalLogDerivation(sum_props, &props, &wells[wells.size() - 1])));
		}
		else if (w_type == "frac2d")
		{
			props.sizes.z *= sqrt(props.kx / props.kz);
			geom_props.rc.z *= sqrt(props.kx / props.kz);

			geom_props.r1 = geom_props.r2 = geom_props.rc;
			geom_props.r1.x -= geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r2.x += geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r1.y -= geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r2.y += geom_props.length * sin(geom_props.alpha) / 2.0;

			perf_length += geom_props.length;
			wells.push_back(Well(geom_props, WellType::FRAC, w_name, seg_idx++));
			summators.push_back(static_cast<BaseSum*>(new Frac2dSum(sum_props, &props, &wells[wells.size() - 1])));
		}
		else if (w_type == "frac2d_inside")
		{
			props.sizes.z *= sqrt(props.kx / props.kz);
			geom_props.rc.z *= sqrt(props.kx / props.kz);

			geom_props.r1 = geom_props.r2 = geom_props.rc;
			geom_props.r1.x -= geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r2.x += geom_props.length * cos(geom_props.alpha) / 2.0;
			geom_props.r1.y -= geom_props.length * sin(geom_props.alpha) / 2.0;
			geom_props.r2.y += geom_props.length * sin(geom_props.alpha) / 2.0;

			perf_length += geom_props.length;
			wells.push_back(Well(geom_props, WellType::FRAC, w_name, seg_idx++));
			summators.push_back(static_cast<BaseSum*>(new Frac2dSum(sum_props, &props, &wells[wells.size() - 1])));
			insideSum = new InsideFrac1d(sum_props, &props, &wells[wells.size() - 1]);
		}

	}

	if (seg_num > 1)
	{
			matSize = (seg_num - 1) * (seg_num - 1);
			q = new double[seg_num];
			dq = new double[seg_num];
			b = new double[seg_num - 1];
			x = new double[seg_num - 1];

			a = new double[matSize];
			ind_i = new int[matSize];
			ind_j = new int[matSize];
			for (int i = 0; i < matSize; i++)
			{
				ind_i[i] = i / (seg_num - 1);
				ind_j[i] = i % (seg_num - 1);
			}

			dpdq = new double*[seg_num];
			if (insideSum != nullptr)
				dpdq2 = new double*[seg_num];
			for (int i = 0; i < seg_num; i++)
			{
				dpdq[i] = new double[seg_num - 1];
				if (insideSum != nullptr)
					dpdq2[i] = new double[seg_num - 1];
			}

			//A.AllocateDENSE("A", props.K - 1, props.K - 1);
			X.Allocate("X", seg_num - 1);
			B.Allocate("B", seg_num - 1);
	}
}
void WellFlow::loadGeomProps(tinyxml2::XMLElement* xml_well_props, WellGeomProperties* geom_props)
{
	XMLElement* xml_geom = xml_well_props->FirstChildElement("geometry");
	XMLElement* xml_rc = xml_geom->FirstChildElement("rc");
	geom_props->rc.x = stod(xml_rc->Attribute("x")) / props.x_dim;
	geom_props->rc.y = stod(xml_rc->Attribute("y")) / props.x_dim;
	geom_props->rc.z = stod(xml_rc->Attribute("z")) / props.x_dim;
	XMLElement* xml_length = xml_geom->FirstChildElement("length");
	geom_props->length = stod(xml_length->Attribute("value")) / props.x_dim;
	XMLElement* xml_alpha = xml_geom->FirstChildElement("alpha");
	geom_props->alpha = stod(xml_alpha->Attribute("value")) * M_PI / 180.0;
	XMLElement* xml_rw = xml_geom->FirstChildElement("rw");
	geom_props->rw = stod(xml_rw->Attribute("value")) / props.x_dim;
}
void WellFlow::loadSlantedGeomProps(tinyxml2::XMLElement* xml_well_props, WellGeomProperties* geom_props)
{
	XMLElement* xml_geom = xml_well_props->FirstChildElement("geometry");
	XMLElement* xml_r1 = xml_geom->FirstChildElement("r1");
	geom_props->r1.x = stod(xml_r1->Attribute("x")) / props.x_dim;
	geom_props->r1.y = stod(xml_r1->Attribute("y")) / props.x_dim;
	geom_props->r1.z = stod(xml_r1->Attribute("z")) / props.x_dim;
	XMLElement* xml_r2 = xml_geom->FirstChildElement("r2");
	geom_props->r2.x = stod(xml_r2->Attribute("x")) / props.x_dim;
	geom_props->r2.y = stod(xml_r2->Attribute("y")) / props.x_dim;
	geom_props->r2.z = stod(xml_r2->Attribute("z")) / props.x_dim;
	geom_props->length = sqrt((geom_props->r2 - geom_props->r1) * (geom_props->r2 - geom_props->r1));
	geom_props->rc = (geom_props->r1 + geom_props->r2) / 2.0;
	geom_props->alpha = 0.0;
	XMLElement* xml_rw = xml_geom->FirstChildElement("rw");
	geom_props->rw = stod(xml_rw->Attribute("value")) / props.x_dim;
}
void WellFlow::loadSumProps(tinyxml2::XMLElement* xml_sum_props, SummatorProperties* sum_props)
{
	sum_props->M = stoi(xml_sum_props->Attribute("M"));
	sum_props->N = stoi(xml_sum_props->Attribute("N"));
	sum_props->L = stoi(xml_sum_props->Attribute("L"));
	sum_props->I = stoi(xml_sum_props->Attribute("I"));
	sum_props->xi_c = stod(xml_sum_props->Attribute("xi_c")) / props.x_dim / props.x_dim;
}

const MainProperties* WellFlow::getProps() const
{
	return &props;
}
const Well* WellFlow::getWell(const int i) const
{
	return &wells[i];
}
const BaseSum* WellFlow::getSummator(const int i) const
{
	return summators[i];
}

void WellFlow::findRateDistribution()
{
	// Fills the vector of rates
	auto fill_q = [this]() 
	{
		for (int i = 0; i < seg_num; i++)
			q[i] = segs[i]->rate;
	};
	// Fills the vector of rate's deviations with zeros
	auto fill_dq = [this]() {
		for (int i = 0; i < seg_num; i++)
			dq[i] = 0.0;
	};
	// Set rate deviation
	auto setRateDev = [this](const int seg_idx, double ratio) {
		segs[seg_idx]->rate += props.rate * ratio;
	};
	// Fills dpdq matrix
	auto fill_dpdq = [&, this](double mult) {
		double p1, p2, ratio;
		ratio = mult * 0.001 / (double)(seg_num);

		for (int i = 0; i < seg_num; i++)
		{
			for (int j = 1; j < seg_num; j++)
			{
				setRateDev(j, -ratio);	setRateDev(0, ratio);
				p1 = getPres(i);

				setRateDev(j, 2.0 * ratio);	setRateDev(0, -2.0 * ratio);
				p2 = getPres(i);

				setRateDev(j, -ratio);	setRateDev(0, ratio);

				dpdq[i][j - 1] = (p2 - p1) / (2.0 * ratio * props.rate);
			}
		}
	};
	auto solve_sys = [this]() {
		double s, p1, p2;

		for (int i = 0; i < seg_num - 1; i++)
		{
			for (int j = 0; j < seg_num - 1; j++)
			{
				s = 0.0;
				for (int k = 0; k < seg_num - 1; k++)
					s += (dpdq[k + 1][j] - dpdq[k][j]) * (dpdq[k + 1][i] - dpdq[k][i]);

				a[i * (seg_num - 1) + j] = s;

			}

			s = 0.0;
			for (int k = 0; k < seg_num - 1; k++)
			{
				p1 = getPres(k);
				p2 = getPres(k+1);
				s += (p2 - p1) * (dpdq[k + 1][i] - dpdq[k][i]);
			}
			B[i] = -s;
		}

		X.Zeros();
		A.Assemble(ind_i, ind_j, a, (seg_num - 1) * (seg_num - 1), "A", seg_num - 1, seg_num - 1);
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
		for (int i = 0; i < seg_num - 1; i++)
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

	double H0 = pres_dev;
	if (sqrt(H0) > 0.001 * fabs(pres_av))
	{
		for(const auto& well : wells)
			well.printRates(&props);

		fill_q();

		double mult = 1;
		double H = H0;

		while (H > H0 / 200.0 || (sqrt(H) > 0.002 * fabs(pres_av)))
		{
			solve_dq(mult);

			double tmp;
			for (int i = 0; i < seg_num; i++)
			{
				tmp = q[i] + mult * dq[i];
				segs[i]->rate = q[i] = tmp;
			}

			calcPressure();
			for (const auto& well : wells)
				well.printRates(&props);

			H = pres_dev;
		}
	}
	else
		for (const auto& well : wells)
			well.printRates(&props);
}
void WellFlow::findRateDistributionInside()
{
	// Fills the vector of rates
	auto fill_q = [this]()
	{
		for (int i = 0; i < seg_num; i++)
			q[i] = segs[i]->rate;
	};
	// Fills the vector of rate's deviations with zeros
	auto fill_dq = [this]() {
		for (int i = 0; i < seg_num; i++)
			dq[i] = 0.0;
	};
	// Set rate deviation
	auto setRateDev = [this](const int seg_idx, double ratio) {
		segs[seg_idx]->rate += props.rate * ratio;
	};
	// Fills dpdq matrix
	auto fill_dpdq = [&, this](double mult) {
		double p1, p2, ratio;
		double p11, p22;
		ratio = mult * 0.001 / (double)(seg_num);

		for (int i = 0; i < seg_num; i++)
		{
			for (int j = 1; j < seg_num; j++)
			{
				setRateDev(j, -ratio);	setRateDev(0, ratio);
				p1 = getPres(i);
				p11 = getPresInside(i);

				setRateDev(j, 2.0 * ratio);	setRateDev(0, -2.0 * ratio);
				p2 = getPres(i);
				p22 = getPresInside(i);

				setRateDev(j, -ratio);	setRateDev(0, ratio);

				dpdq[i][j - 1] = (p2 - p1) / (2.0 * ratio * props.rate);
				dpdq2[i][j - 1] = (p22 - p11) / (2.0 * ratio * props.rate);
			}
		}
	};
	auto solve_sys = [this]() {
		double s, p1, p2;

		for (int i = 0; i < seg_num - 1; i++)
		{
			for (int j = 0; j < seg_num - 1; j++)
			{
				s = 0.0;
				for (int k = 0; k < seg_num; k++)
					s += (dpdq[k][j] - dpdq2[k][j]) * (dpdq[k][i] - dpdq2[k][i]);

				a[i * (seg_num - 1) + j] = s;

			}

			s = 0.0;
			for (int k = 0; k < seg_num; k++)
			{
				p1 = getPres(k);
				p2 = getPresInside(k);
				s += (p1 - p2) * (dpdq[k][i] - dpdq2[k][i]);
			}
			B[i] = -s;
		}

		X.Zeros();
		A.Assemble(ind_i, ind_j, a, (seg_num - 1) * (seg_num - 1), "A", seg_num - 1, seg_num - 1);
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
		for (int i = 0; i < seg_num - 1; i++)
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

	double H0 = pres_dev;
	if (sqrt(H0) > 0.001 * fabs(pres_av))
	{
		for (const auto& well : wells)
			well.printRates(&props);

		fill_q();

		double mult = 1;
		double H = H0;

		while (H > H0 / 200.0 || (sqrt(H) > 0.002 * fabs(pres_av)))
		{
			solve_dq(mult);

			double tmp;
			for (int i = 0; i < seg_num; i++)
			{
				tmp = q[i] + mult * dq[i];
				segs[i]->rate = q[i] = tmp;
			}

			calcPressure();
			for (const auto& well : wells)
				well.printRates(&props);

			H = pres_dev;
		}
	}
	else
		for (const auto& well : wells)
			well.printRates(&props);
}
void WellFlow::calcPressure()
{
	pres_av = pres_dev = 0.0;
	double p2d, p3d;
	double* inside2d = NULL;
	for (int seg_idx = 0; seg_idx < seg_num; seg_idx++)
	{
		WellSegment* seg = segs[seg_idx];
		seg->pres2D = seg->pres3D = seg->pres = 0.0;
		for (auto& summator : summators)
		{
			p2d = summator->get2D(seg_idx);
			p3d = summator->get3D(seg_idx);
			seg->pres2D += summator->get2D(seg_idx);
			seg->pres3D += summator->get3D(seg_idx);
			seg->pres += p2d + p3d;
		}

		pres_av += seg->pres;
	}
	pres_av /= seg_num;

	if (insideSum != nullptr)
	{
		Point pt1({ segs[0]->r1.x - wells[0].getGeomProps()->rw, segs[0]->r1.y, segs[0]->r1.z });
		Point pt2({ segs[segs.size()-1]->r2.x + wells[0].getGeomProps()->rw, segs[segs.size()-1]->r2.y, segs[segs.size() - 1]->r2.z });
		insideSum->setBounds(summators[0]->getPressure(pt1) , summators[0]->getPressure(pt2));
		inside2d = new double[seg_num];
		for (size_t k = 0; k < seg_num; k++)
			inside2d[k] = insideSum->get2D(k);
		for (int i = 0; i < seg_num; i++)
			std::cout << inside2d[i] << std::endl;
	}

	for (auto& summator : summators)
	{
		auto well = summator->getWell();
		well->pres_av = well->pres_dev = 0.0;
		const int K = summator->getSumProps()->K;
		if (K > 1)
		{
			if (insideSum == nullptr)
			{
				for (int i = 1; i < K; i++)
					well->pres_dev += (well->segs[i].pres - well->segs[i - 1].pres) *
						(well->segs[i].pres - well->segs[i - 1].pres);
			} 
			else
			{
				for (int i = 0; i < K; i++)
					well->pres_dev += (well->segs[i].pres - inside2d[i]) *
						(well->segs[i].pres - inside2d[i]);
			}
		}
		well->pres_dev /= 2.0;
	}
	if (seg_num > 1)
	{
		if (insideSum == nullptr)
		{
			for (int seg_idx = 1; seg_idx < seg_num; seg_idx++)
			{
				pres_dev += (segs[seg_idx]->pres - segs[seg_idx - 1]->pres) *
					(segs[seg_idx]->pres - segs[seg_idx - 1]->pres);
			}
		}
		else
		{
			for (int seg_idx = 0; seg_idx < seg_num; seg_idx++)
			{
				pres_dev += (segs[seg_idx]->pres - inside2d[seg_idx]) *
					(segs[seg_idx]->pres - inside2d[seg_idx]);
			}
		}
		pres_dev /= 2.0;
	}
}
double WellFlow::getP_bhp()
{
	for(auto& sum : summators)
		sum->prepare();
	if(insideSum != nullptr)
		insideSum->prepare();

	if (firstTime)
	{
		// Body of function
		for (auto& well : wells)
			well.setUniformRate();
		firstTime = false;
		if (insideSum == nullptr)
			findRateDistribution();
		else
			findRateDistributionInside();
	}

	calcPressure();
	wells[0].writeRates(&props);

	return pres_av;
}
double WellFlow::getPres(const int seg_idx)
{
	double sum = 0.0;
	for (const auto summator : summators)
		sum += summator->getPres(seg_idx);
	return sum;
}
double WellFlow::getPresInside(const int seg_idx)
{
	assert(insideSum != nullptr);
	Point pt1({ segs[0]->r1.x - wells[0].getGeomProps()->rw, segs[0]->r1.y, segs[0]->r1.z });
	Point pt2({ segs[segs.size() - 1]->r2.x + wells[0].getGeomProps()->rw, segs[segs.size() - 1]->r2.y, segs[segs.size() - 1]->r2.z });
	insideSum->setBounds(summators[0]->getPressure(pt1), summators[0]->getPressure(pt2));
	return insideSum->getPres(seg_idx);
}