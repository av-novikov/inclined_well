#include "src/inclined_sum/welltests/HorizontalLogDerivation.hpp"

HorizontalLogDerivation::HorizontalLogDerivation(const SummatorProperties& _sprops, const MainProperties* _props, const Well* _well) : BaseSum(_sprops, _props, _well)
{
}
HorizontalLogDerivation::~HorizontalLogDerivation()
{
}
void HorizontalLogDerivation::prepare()
{
	double bufx, bufy, lambda;

	for (int arr_idx = 0; arr_idx < size; arr_idx++)
	{
		const WellSegment& seg = well->segs[arr_idx % sprops.K];
		const Point& r = (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->r_bhp;

		F2d[arr_idx] = 0.0;

		for (int m = 1; m < sprops.M; m++)
		{
			bufx = sin(M_PI * (double)(m)* r.x / props->sizes.x) / (double)(m) *
				(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x) -
					cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x));
			for (int n = 1; n < sprops.N; n++)
			{
				bufy = sin(M_PI * (double)(n)* r.y / props->sizes.y) *
					sin(M_PI * (double)(n)* gprops->rc.y / props->sizes.y);

				lambda = M_PI * M_PI * props->kx / props->visc / MainProperties::compressibility *
					((double)(m * m) / props->sizes.x / props->sizes.x +
					(double)(n * n) / props->sizes.y / props->sizes.y);
				
				F2d[arr_idx] += bufx * bufy * exp(-lambda * time) / lambda;
			}
		}
	}

	for (int arr_idx = 0; arr_idx < size; arr_idx++)
	{
		const WellSegment& seg = well->segs[arr_idx % sprops.K];
		const Point& r = (*segs)[int((double)(arr_idx) / (double)(sprops.K))]->r_bhp;

		F3d[arr_idx] = 0.0;

		for (int m = 1; m < sprops.M; m++)
		{
			bufx = sin(M_PI * (double)(m)* r.x / props->sizes.x) / (double)(m) * 
				(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x) -
				cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x));
			for (int n = 1; n < sprops.N; n++)
			{
				bufy = sin(M_PI * (double)(n)* r.y / props->sizes.y) *
					sin(M_PI * (double)(n)* gprops->rc.y / props->sizes.y);
				for (int l = 1; l < sprops.L; l++)
				{
					lambda = M_PI * M_PI * props->kx / props->visc / MainProperties::compressibility * 
						((double)(m * m) / props->sizes.x / props->sizes.x +
						(double)(n * n) / props->sizes.y / props->sizes.y +
						(double)(l * l) / props->sizes.z / props->sizes.z);
					F3d[arr_idx] += bufx * bufy * exp(-lambda * time) / lambda * 
							cos(M_PI * (double)(l)* r.z / props->sizes.z) * 
							cos(M_PI * (double)(l)* gprops->rc.z / props->sizes.z);
				}
			}
		}
	}
}
double HorizontalLogDerivation::get2D(int seg_idx) const
{
	double sum = 0.0;

	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F2d[seg_idx * sprops.K + k] * seg.rate / seg.length;
	}

	sum *= (4.0 / M_PI / props->sizes.y / props->sizes.z / MainProperties::compressibility);

	return sum;
}
double HorizontalLogDerivation::get3D(int seg_idx) const
{
	double sum = 0.0;

	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F3d[seg_idx * sprops.K + k] * seg.rate / seg.length;
	}

	sum *= (8.0 / M_PI / props->sizes.y / props->sizes.z / MainProperties::compressibility);

	return sum;
}
double HorizontalLogDerivation::getLogDerivative()
{
	double sum = 0.0;
	double bufx, bufy, lambda;

	for (int k = 0; k < sprops.K; k++)
	{
		const WellSegment seg = well->segs[k];
		const Point& r = well->segs[k].r_bhp;
		
		for (int m = 1; m < sprops.M; m++)
		{
			bufx = sin(M_PI * (double)(m)* r.x / props->sizes.x) / (double)(m) *
				(cos(M_PI * (double)(m)* seg.r1.x / props->sizes.x) -
					cos(M_PI * (double)(m)* seg.r2.x / props->sizes.x));
			for (int n = 1; n < sprops.N; n++)
			{
				bufy = sin(M_PI * (double)(n)* r.y / props->sizes.y) *
					sin(M_PI * (double)(n)* gprops->rc.y / props->sizes.y);
				for (int l = 0; l < sprops.L; l++)
				{
					lambda = M_PI * M_PI * props->kx / props->visc / MainProperties::compressibility *
						((double)(m * m) / props->sizes.x / props->sizes.x +
						(double)(n * n) / props->sizes.y / props->sizes.y +
							(double)(l * l) / props->sizes.z / props->sizes.z);
					sum	+= seg.rate / seg.length * 
						bufx * bufy * exp(-lambda * time) * time *
						cos(M_PI * (double)(l)* r.z / props->sizes.z) *
						cos(M_PI * (double)(l)* gprops->rc.z / props->sizes.z);
				}
			}
		}
	}

	sum *= (8.0 / M_PI / props->sizes.y / props->sizes.z / MainProperties::compressibility);

	return sum;
}
void HorizontalLogDerivation::setTime(const double _time)
{
	time = _time / props->t_dim;
}