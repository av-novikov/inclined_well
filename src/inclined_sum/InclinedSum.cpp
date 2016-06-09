#include "src/inclined_sum/InclinedSum.hpp"

#define TAG_EMPTY_MES 100000000

using namespace MPI;

InclinedSum::InclinedSum(const Parameters* _props, const Well* _well) : props(_props), well(_well)
{
	rank = COMM_WORLD.Get_rank();
	size = COMM_WORLD.Get_size();
	
	startIdx = 1 + int((double)(rank) / (double)(size) * (double)(props->M));
	finishIdx = int((double)(rank + 1) / (double)(size) * (double)(props->M));
}

InclinedSum::~InclinedSum()
{
}

double InclinedSum::get2D(const Point& r)
{
	double sum = 0.0;
	double sum_prev = 0.0;
	int break_idx = 0;

	for(int k = 0; k < props->K; k++)
	{
		break_idx = 0;
		
		for(int m = startIdx; m <= finishIdx; m++)
		{
			for(int i = -props->I; i <= props->I; i++)
			{
				sum += well->segs[k].rate / well->segs[k].length / (double)(m) / (double)(m) * 
					( exp(-M_PI * (double)(m) * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) - 
					exp(-M_PI * (double)(m) * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y) / props->sizes.x) ) *
					sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
					(cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x) -
					cos(M_PI * (double)(m) * (well->segs[k].r1.x + tan(props->alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props->sizes.x));
			}
			
			if( fabs(sum - sum_prev) > sum * EQUALITY_TOLERANCE )
			{
				sum_prev = sum;
				break_idx = 0;
			} else
				break_idx++;
				
			if(break_idx > 1)
			{
				std::cout << m << std::endl;
				break;
			}
			
		}
	}
	
	sum *= (props->visc * props->sizes.x / M_PI / M_PI / props->sizes.z / props->perm / sin(props->alpha));
	
	if(size > 1)
	{
		if(rank == 0)
		{
			MPI_Status status;
			double sum_buf;
		
			for(int i = 1; i < size; i++)
			{
				MPI_Recv(&sum_buf, 1, MPI_DOUBLE, i, TAG_EMPTY_MES + i, MPI_COMM_WORLD, &status);
				sum += sum_buf;
			}	
		} else
			MPI_Send(&sum, 1, MPI_DOUBLE, 0, TAG_EMPTY_MES + rank, MPI_COMM_WORLD);
	}
		
	COMM_WORLD.Barrier();

	return sum;	
}

double InclinedSum::get3D(const Point& r)
{
	double sum = 0.0;
	double sum_prev1 = 0.0, sum_prev2 = 0.0;
	int break_idx1, break_idx2;
	double F, buf;
	
	for(int k = 0; k < props->K; k++)
	{
		break_idx2 = 0;
		
		for(int m = startIdx; m <= finishIdx; m++)
		{
			break_idx1 = 0;
			
			for(int l = 1; l <= props->L; l++)
			{
				buf = sqrt((double)(m) * (double)(m) / props->sizes.x / props->sizes.x + (double)(l) * (double)(l) / props->sizes.z / props->sizes.z);
				
				F = ((cos(M_PI * (double)(m) * well->segs[k].r2.x / props->sizes.x + M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
					cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x + M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							(M_PI * (double)(l) / props->sizes.z - M_PI * (double)(m) * tan(props->alpha) / props->sizes.x) - 
					(cos(M_PI * (double)(m) * well->segs[k].r2.x / props->sizes.x - M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
					cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x - M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							(M_PI * (double)(l) / props->sizes.z + M_PI * (double)(m) * tan(props->alpha) / props->sizes.x)
					) / 2.0;
				
				for(int i = -props->I; i <= props->I; i++)
				{
					/*F = ((cos(M_PI * (double)(m) * (well->segs[k].r1.x + 
							tan(props->alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props->sizes.x + 
							M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
						cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x +
							M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							(M_PI * (double)(l) / props->sizes.z - M_PI * (double)(m) * tan(props->alpha) / props->sizes.x) - 
						(cos(M_PI * (double)(m) * (well->segs[k].r1.x + 
							tan(props->alpha) * (well->segs[k].r1.z - well->segs[k].r2.z)) / props->sizes.x - 
							M_PI * (double)(l) * well->segs[k].r1.z / props->sizes.z) -
						cos(M_PI * (double)(m) * well->segs[k].r1.x / props->sizes.x -
							M_PI * (double)(l) * well->segs[k].r2.z / props->sizes.z)) /
							(M_PI * (double)(l) / props->sizes.z + M_PI * (double)(m) * tan(props->alpha) / props->sizes.x)
							) / 2.0;*/
						
					sum += F * well->segs[k].rate / well->segs[k].length / buf * 
						( exp(-M_PI * buf * fabs(r.y - props->r1.y + 2.0 * (double)(i) * props->sizes.y)) - 
						exp(-M_PI * buf * fabs(r.y + props->r1.y + 2.0 * (double)(i) * props->sizes.y)) ) *
						sin(M_PI * (double)(m) * r.x / props->sizes.x) * 
						cos(M_PI * (double)(l) * r.z / props->sizes.z);
				}
				
				if( fabs(sum - sum_prev1) > sum * EQUALITY_TOLERANCE )
				{
					sum_prev1 = sum;
					break_idx1 = 0;
				} else
					break_idx1++;
				
				if(break_idx1 > 1)
				{
					//std::cout << "l=" << l << std::endl;
					break;
				}
			}
			
			if( fabs(sum - sum_prev2) > sum * EQUALITY_TOLERANCE )
			{
				sum_prev2 = sum;
				break_idx2 = 0;
			} else
				break_idx2++;
				
			if(break_idx2 > 1)
			{
				std::cout << "m=" << m << std::endl;
				break;
			}
		}
	}
	
	if(size > 1)
	{
		if(rank == 0)
		{
			MPI_Status status;
			double sum_buf;
		
			for(int i = 1; i < size; i++)
			{
				MPI_Recv(&sum_buf, 1, MPI_DOUBLE, i, TAG_EMPTY_MES + i, MPI_COMM_WORLD, &status);
				sum += sum_buf;
			}	
		} else
			MPI_Send(&sum, 1, MPI_DOUBLE, 0, TAG_EMPTY_MES + rank, MPI_COMM_WORLD);
	}
		
	COMM_WORLD.Barrier();
				
	sum *= (2.0 * props->visc / M_PI / props->sizes.x / props->sizes.z / props->perm / cos(props->alpha));
				
	return sum;	
}

double InclinedSum::getPres(const Point& r)
{
	//return get2D(r);
	double s1, s2;
	s1 = get2D(r);
	s2 = get3D(r);
	//std::cout << "2d = " << s1 << "\t3d = " << s2 << std::endl;
	return s1 + s2;
}
