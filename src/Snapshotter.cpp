#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>

#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>

#include <mpi.h>

#include "src/Snapshotter.hpp"
#include "src/WellFlow.hpp"

using std::string;

Snapshotter::Snapshotter(const Point& _sizes, const int _nx, const int _ny, const int _nz, double _x_dim) :
			sizes(_sizes), nx(_nx), ny(_ny), nz(_nz), x_dim(_x_dim)
{
	steps.x = sizes.x / (double)(nx);
	steps.y = sizes.y / (double)(ny);
	steps.z = sizes.z / (double)(nz);
	
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();
	
	startIdx = int((double)(rank) / (double)(size) * (double)(nx));
	finishIdx = int((double)(rank + 1) / (double)(size) * (double)(nx));
}

Snapshotter::~Snapshotter()
{
}

void Snapshotter::snapshot(string name)
{
	vtkSmartPointer<vtkStructuredGrid> grid =
		vtkSmartPointer<vtkStructuredGrid>::New();
		
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	
	
	grid->SetDimensions(finishIdx - startIdx + 1, ny + 1, nz + 1);
	for(int k = 0; k < nz+1; k++)		
		for(int j = 0; j < ny+1; j++)		
			for(int i = startIdx; i < finishIdx+1; i++)		
			{
				points->InsertNextPoint( x_dim * (double)(i) * steps.x, x_dim * (double)(j) * steps.y, -x_dim * (double)(k) * steps.z );
			}		
	grid->SetPoints(points);
	
	vtkSmartPointer<vtkDoubleArray> pres = 
		vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	
	for(int k = 0; k < nz+1; k++)		
		for(int j = 0; j < ny+1; j++)		
			for(int i = startIdx; i < finishIdx+1; i++)		
			{
				Point point( (double)(i) * steps.x, (double)(j) * steps.y, -(double)(k) * steps.z );
				//pres->InsertNextValue( presFoo(point) / BAR );
			}
			
	vtkPointData* fd = grid->GetPointData();
	fd->AddArray( pres );
	
	vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
	
	writer->SetFileName( name.c_str() );
	writer->SetInputData(grid);
	writer->Write();
}

void Snapshotter::setPresFoo(const FooType& foo)
{
	presFoo = foo;
}
