#include "src/inclined_sum/welltests/HorizontalLogDerivation.hpp"

HorizontalLogDerivation::HorizontalLogDerivation(const Parameters* _props, const Well* _well) : BaseSum(_props, _well)
{
}

HorizontalLogDerivation::~HorizontalLogDerivation()
{
}

double HorizontalLogDerivation::getPres(const Point& point)
{
	return 0.0;
}

void HorizontalLogDerivation::prepare()
{
	double bufx, bufy, lambda;

	for (int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[arr_idx % props->K];
		const Point& r = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;

		F2d[arr_idx] = 0.0;

		for (int m = 1; m < props->M; m++)
		{
			bufx = sin(M_PI * (double)(m)* r.x / props->sizes.x) / (double)(m) *
				(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x) -
					cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x));
			for (int n = 1; n < props->N; n++)
			{
				bufy = sin(M_PI * (double)(n)* r.y / props->sizes.y) *
					sin(M_PI * (double)(n)* props->rc.y / props->sizes.y);

				lambda = M_PI * M_PI * props->kx / props->visc / Parameters::compressibility *
					((double)(m * m) / props->sizes.x / props->sizes.x +
					(double)(n * n) / props->sizes.y / props->sizes.y);
				
				F2d[arr_idx] += bufx * bufy * exp(-lambda * time) / lambda;
			}
		}
	}

	for (int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[arr_idx % props->K];
		const Point& r = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;

		F3d[arr_idx] = 0.0;

		for (int m = 1; m < props->M; m++)
		{
			bufx = sin(M_PI * (double)(m)* r.x / props->sizes.x) / (double)(m) * 
				(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x) -
				cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x));
			for (int n = 1; n < props->N; n++)
			{
				bufy = sin(M_PI * (double)(n)* r.y / props->sizes.y) *
					sin(M_PI * (double)(n)* props->rc.y / props->sizes.y);
				for (int l = 1; l < props->L; l++)
				{
					lambda = M_PI * M_PI * props->kx / props->visc / Parameters::compressibility * 
						((double)(m * m) / props->sizes.x / props->sizes.x +
						(double)(n * n) / props->sizes.y / props->sizes.y +
						(double)(l * l) / props->sizes.z / props->sizes.z);
					F3d[arr_idx] += bufx * bufy * exp(-lambda * time) / lambda * 
							cos(M_PI * (double)(l)* r.z / props->sizes.z) * 
							cos(M_PI * (double)(l)* props->rc.z / props->sizes.z);
				}
			}
		}
	}
}

double HorizontalLogDerivation::get2D(int seg_idx)
{
	double sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F2d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (-4.0 / M_PI / props->sizes.y / props->sizes.z / Parameters::compressibility);

	return sum;
}

double HorizontalLogDerivation::get3D(int seg_idx)
{
	double sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F3d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (-8.0 / M_PI / props->sizes.y / props->sizes.z / Parameters::compressibility);

	return sum;
}

double HorizontalLogDerivation::getLogDerivative()
{
	double sum = 0.0;
	double bufx, bufy, lambda;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment seg = well->segs[k];
		const Point& r = well->segs[k].r_bhp;
		
		for (int m = 1; m < props->M; m++)
		{
			bufx = sin(M_PI * (double)(m)* r.x / props->sizes.x) / (double)(m) *
				(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x) -
					cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x));
			for (int n = 1; n < props->N; n++)
			{
				bufy = sin(M_PI * (double)(n)* r.y / props->sizes.y) *
					sin(M_PI * (double)(n)* props->rc.y / props->sizes.y);
				for (int l = 0; l < props->L; l++)
				{
					lambda = M_PI * M_PI * props->kx / props->visc / Parameters::compressibility *
						((double)(m * m) / props->sizes.x / props->sizes.x +
						(double)(n * n) / props->sizes.y / props->sizes.y +
							(double)(l * l) / props->sizes.z / props->sizes.z);
					sum	+= seg.rate / seg.length * 
						bufx * bufy * exp(-lambda * time) * time *
						cos(M_PI * (double)(l)* r.z / props->sizes.z) *
						cos(M_PI * (double)(l)* props->rc.z / props->sizes.z);
				}
			}
		}
	}

	sum *= (8.0 / M_PI / props->sizes.y / props->sizes.z / Parameters::compressibility);

	return sum;
}

void HorizontalLogDerivation::setTime(const double _time)
{
	time = _time / props->t_dim;
}