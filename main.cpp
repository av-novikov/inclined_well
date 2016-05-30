#include <iostream>
#include <mpi.h>
#include <tinyxml.h>
#include <string>

#include "WellFlow.hpp"
#include "VerticalWellFlow.hpp"
#include "Grid.hpp"

using namespace std;
using namespace std::placeholders;

Parameters loadTask(string fileName)
{
	Parameters props;
	
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
	props.perm = 1.E-15 * stod( xml_perm->Attribute("value") ) / (props.x_dim * props.x_dim);
	
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
	
	return props;
}

double getPbhp(WellFlow* solver, Parameters* props, bool debug)
{
	const double h_phi = 2.0 * M_PI / 36.0;
	const double h_z = props->sizes.z / props->K;
	
	double z = 0.0;
	double p_bhp = 0.0;
	for(int kz = 0; kz <= props->K; kz++)
		for(int phi = 0; phi < 36; phi++)
		{
			z = (double)(kz) * h_z;
			Point point (props->r1.x + tan(props->alpha) * z + props->rw * cos((double)(phi) * h_phi), props->r1.y + props->rw * sin((double)(phi) * h_phi), -z);
			p_bhp += solver->calcPressure(point, debug);
		}
		
	return p_bhp / (props->K + 1) / 36.0 / BAR;
}

int main()
{
	MPI::Init();
	const int rank = MPI::COMM_WORLD.Get_rank();
	//const int size = MPI::COMM_WORLD.Get_size();
	
	Parameters props = loadTask("config.xml");
	
	WellFlow solver (props);
	
	Grid grid (props.sizes, props.nx, props.ny, props.nz, props.x_dim);
	grid.setPresFoo( bind(&WellFlow::calcPressure, &solver, _1, _2) );
	//grid.snapshot("snap_" + to_string(rank) + ".vts");
	
	if(rank == 0)
	{
		cout << "Well coords:" << endl 
			<< "\t" << props.x_dim * props.r1
			<< "\t" << props.x_dim * props.r2;
		cout << "P_old= " << getPbhp(&solver, &props, true) << endl;
	}
	
	MPI::Finalize();
	
	return 0;
}
