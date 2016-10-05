#include <new>

#include "src/inclined_sum/Inclined3dSum.h"

Inclined3dSum::Inclined3dSum(const Parameters* _props, const Well* _well) : BaseSum(_props, _well)
{
}

Inclined3dSum::~Inclined3dSum()
{
}

double Inclined3dSum::getPres(const Point& point)
{
	return 0.0;
}

double Inclined3dSum::get2D(int seg_idx)
{
	double sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F2d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	return sum;
}

double Inclined3dSum::get3D(int seg_idx)
{
	double sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment& seg = well->segs[k];
		sum += F3d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	return sum;
}

void Inclined3dSum::prepare()
{
	prepare2D();	prepare3D();
}

void Inclined3dSum::prepare2D()
{
	for (int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		F2d[arr_idx] = 0.0;
	}
}

void Inclined3dSum::prepare3D()
{
	double a[2], b[2], c[2], d[2], e[2], f[2];
	double A, B1, B2;
	double sum1, sum2;
	double buf1, buf2, buf3;
	int k;

	auto intExpansion = [this](double A, double B) {
		return 2.0 / sqrt(M_PI) * (1.0 / 2.0 * sqrt(M_PI * B * B / A) * erf(sqrt(A / props->xi_c)) * (2.0 + B * B / A / 3.0 + 3.0 * B * B * B * B / A / A / 20.0) - 
				B * B * B / A / sqrt(props->xi_c) * exp(-A / props->xi_c) * (1.0 / 3.0 + B * B / props->xi_c / 10.0 + 3.0 * B * B / A / 20.0));
		/*return 2.0 / sqrt(M_PI) * B / sqrt(props->xi_c) *
			(2.0 - 2.0 * B * B / 9.0 / props->xi_c * (3.0 * A / B / B + 1.0) +
			B * B * B * B / 75.0 / props->xi_c / props->xi_c * (15.0 * A * A / B / B / B / B + 10.0 * A / B / B + 3.0));*/
	};

	for (int arr_idx = 0; arr_idx < props->K * props->K; arr_idx++)
	{
		const WellSegment seg = well->segs[arr_idx % props->K];
		const Point& point = well->segs[int((double)(arr_idx) / (double)(props->K))].r_bhp;

		F3d[arr_idx] = sum1 = sum2 = 0.0;
	
		for (int p = -props->I; p <= props->I; p++)
		{
			for (int q = -props->I; q <= props->I; q++)
			{
				for (int r = -props->I; r <= props->I; r++)
				{
					a[0] = props->r1.x - props->r2.x;	b[0] = point.x + 2.0 * (double)(p)* props->sizes.x - props->r1.x;
					a[1] = props->r2.x - props->r1.x;	b[1] = point.x + 2.0 * (double)(p)* props->sizes.x + props->r1.x;
					c[0] = props->r1.y - props->r2.y;	d[0] = point.y + 2.0 * (double)(q)* props->sizes.y - props->r1.y;
					c[1] = props->r2.y - props->r1.y;	d[1] = point.y + 2.0 * (double)(q)* props->sizes.y + props->r1.y;
					e[0] = props->r1.z - props->r2.z;	f[0] = point.z + 2.0 * (double)(r)* props->sizes.z - props->r1.z;
					e[1] = props->r2.z - props->r1.z;	f[1] = point.z + 2.0 * (double)(r)* props->sizes.z + props->r1.z;

					// 0-0-0
					A = (a[0] * a[0] * (d[0] * d[0] + f[0] * f[0]) -
						2.0 * a[0] * b[0] * (c[0] * d[0] + e[0] * f[0]) +
						b[0] * b[0] * (c[0] * c[0] + e[0] * e[0]) +
						(e[0] * d[0] - c[0] * f[0]) * (e[0] * d[0] - c[0] * f[0])) / 
						(4.0 * (a[0] * a[0] + c[0] * c[0] + e[0] * e[0]));
					B1 = a[0] * (a[0] * seg.tau1 + b[0]) + c[0] * (c[0] * seg.tau1 + d[0]) + e[0] * (e[0] * seg.tau1 + f[0]) / 
						(2.0 * sqrt(a[0] * a[0] + c[0] * c[0] + e[0] * e[0]));
					B2 = a[0] * (a[0] * seg.tau2 + b[0]) + c[0] * (c[0] * seg.tau2 + d[0]) + e[0] * (e[0] * seg.tau2 + f[0]) /
						(2.0 * sqrt(a[0] * a[0] + c[0] * c[0] + e[0] * e[0]));

					sum1 += sqrt(M_PI / (a[0] * a[0] + c[0] * c[0] + e[0] * e[0])) *
							(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

					// 0-0-1
					A = (a[0] * a[0] * (d[0] * d[0] + f[1] * f[1]) -
						2.0 * a[0] * b[0] * (c[0] * d[0] + e[1] * f[1]) +
						b[0] * b[0] * (c[0] * c[0] + e[1] * e[1]) +
						(e[1] * d[0] - c[0] * f[1]) * (e[1] * d[0] - c[0] * f[1])) /
						(4.0 * (a[0] * a[0] + c[0] * c[0] + e[1] * e[1]));
					B1 = a[0] * (a[0] * seg.tau1 + b[0]) + c[0] * (c[0] * seg.tau1 + d[0]) + e[1] * (e[1] * seg.tau1 + f[1]) /
						(2.0 * sqrt(a[0] * a[0] + c[0] * c[0] + e[1] * e[1]));
					B2 = a[0] * (a[0] * seg.tau2 + b[0]) + c[0] * (c[0] * seg.tau2 + d[0]) + e[1] * (e[1] * seg.tau2 + f[1]) /
						(2.0 * sqrt(a[0] * a[0] + c[0] * c[0] + e[1] * e[1]));

					sum1 += sqrt(M_PI / (a[0] * a[0] + c[0] * c[0] + e[1] * e[1])) *
						(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

					// 0-1-0
					A = (a[0] * a[0] * (d[1] * d[1] + f[0] * f[0]) -
						2.0 * a[0] * b[0] * (c[1] * d[1] + e[0] * f[0]) +
						b[0] * b[0] * (c[1] * c[1] + e[0] * e[0]) +
						(e[0] * d[1] - c[1] * f[0]) * (e[0] * d[1] - c[1] * f[0])) /
						(4.0 * (a[0] * a[0] + c[1] * c[1] + e[0] * e[0]));
					B1 = a[0] * (a[0] * seg.tau1 + b[0]) + c[1] * (c[1] * seg.tau1 + d[1]) + e[0] * (e[0] * seg.tau1 + f[0]) /
						(2.0 * sqrt(a[0] * a[0] + c[1] * c[1] + e[0] * e[0]));
					B2 = a[0] * (a[0] * seg.tau2 + b[0]) + c[1] * (c[1] * seg.tau2 + d[1]) + e[0] * (e[0] * seg.tau2 + f[0]) /
						(2.0 * sqrt(a[0] * a[0] + c[1] * c[1] + e[0] * e[0]));

					sum1 += sqrt(M_PI / (a[0] * a[0] + c[1] * c[1] + e[0] * e[0])) *
						(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

					// 0-1-1
					A = (a[0] * a[0] * (d[1] * d[1] + f[1] * f[1]) -
						2.0 * a[0] * b[0] * (c[1] * d[1] + e[1] * f[1]) +
						b[0] * b[0] * (c[1] * c[1] + e[1] * e[1]) +
						(e[1] * d[1] - c[1] * f[1]) * (e[1] * d[1] - c[1] * f[1])) /
						(4.0 * (a[0] * a[0] + c[1] * c[1] + e[1] * e[1]));
					B1 = a[0] * (a[0] * seg.tau1 + b[0]) + c[1] * (c[1] * seg.tau1 + d[1]) + e[1] * (e[1] * seg.tau1 + f[1]) /
						(2.0 * sqrt(a[0] * a[0] + c[1] * c[1] + e[1] * e[1]));
					B2 = a[0] * (a[0] * seg.tau2 + b[0]) + c[1] * (c[1] * seg.tau2 + d[1]) + e[1] * (e[1] * seg.tau2 + f[1]) /
						(2.0 * sqrt(a[0] * a[0] + c[1] * c[1] + e[1] * e[1]));

					sum1 += sqrt(M_PI / (a[0] * a[0] + c[1] * c[1] + e[1] * e[1])) *
						(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

					// 1-0-0
					A = (a[1] * a[1] * (d[0] * d[0] + f[0] * f[0]) -
						2.0 * a[1] * b[1] * (c[0] * d[0] + e[0] * f[0]) +
						b[1] * b[1] * (c[0] * c[0] + e[0] * e[0]) +
						(e[0] * d[0] - c[0] * f[0]) * (e[0] * d[0] - c[0] * f[0])) /
						(4.0 * (a[1] * a[1] + c[0] * c[0] + e[0] * e[0]));
					B1 = a[1] * (a[1] * seg.tau1 + b[1]) + c[0] * (c[0] * seg.tau1 + d[0]) + e[0] * (e[0] * seg.tau1 + f[0]) /
						(2.0 * sqrt(a[1] * a[1] + c[0] * c[0] + e[0] * e[0]));
					B2 = a[1] * (a[1] * seg.tau2 + b[1]) + c[0] * (c[0] * seg.tau2 + d[0]) + e[0] * (e[0] * seg.tau2 + f[0]) /
						(2.0 * sqrt(a[1] * a[1] + c[0] * c[0] + e[0] * e[0]));

					sum1 += sqrt(M_PI / (a[1] * a[1] + c[0] * c[0] + e[0] * e[0])) *
						(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

					// 1-0-1
					A = (a[1] * a[1] * (d[0] * d[0] + f[1] * f[1]) -
						2.0 * a[1] * b[1] * (c[0] * d[0] + e[1] * f[1]) +
						b[1] * b[1] * (c[0] * c[0] + e[1] * e[1]) +
						(e[1] * d[0] - c[0] * f[1]) * (e[1] * d[0] - c[0] * f[1])) /
						(4.0 * (a[1] * a[1] + c[0] * c[0] + e[1] * e[1]));
					B1 = a[1] * (a[1] * seg.tau1 + b[1]) + c[0] * (c[0] * seg.tau1 + d[0]) + e[1] * (e[1] * seg.tau1 + f[1]) /
						(2.0 * sqrt(a[1] * a[1] + c[0] * c[0] + e[1] * e[1]));
					B2 = a[1] * (a[1] * seg.tau2 + b[1]) + c[0] * (c[0] * seg.tau2 + d[0]) + e[1] * (e[1] * seg.tau2 + f[1]) /
						(2.0 * sqrt(a[1] * a[1] + c[0] * c[0] + e[1] * e[1]));

					sum1 += sqrt(M_PI / (a[1] * a[1] + c[0] * c[0] + e[1] * e[1])) *
						(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

					// 1-1-0
					A = (a[1] * a[1] * (d[1] * d[1] + f[0] * f[0]) -
						2.0 * a[1] * b[1] * (c[1] * d[1] + e[0] * f[0]) +
						b[1] * b[1] * (c[1] * c[1] + e[0] * e[0]) +
						(e[0] * d[1] - c[1] * f[0]) * (e[0] * d[1] - c[1] * f[0])) /
						(4.0 * (a[1] * a[1] + c[1] * c[1] + e[0] * e[0]));
					B1 = a[1] * (a[1] * seg.tau1 + b[1]) + c[1] * (c[1] * seg.tau1 + d[1]) + e[0] * (e[0] * seg.tau1 + f[0]) /
						(2.0 * sqrt(a[1] * a[1] + c[1] * c[1] + e[0] * e[0]));
					B2 = a[1] * (a[1] * seg.tau2 + b[1]) + c[1] * (c[1] * seg.tau2 + d[1]) + e[0] * (e[0] * seg.tau2 + f[0]) /
						(2.0 * sqrt(a[1] * a[1] + c[1] * c[1] + e[0] * e[0]));

					sum1 += sqrt(M_PI / (a[1] * a[1] + c[1] * c[1] + e[0] * e[0])) *
						(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

					// 1-1-1
					A = (a[1] * a[1] * (d[1] * d[1] + f[1] * f[1]) -
						2.0 * a[1] * b[1] * (c[1] * d[1] + e[1] * f[1]) +
						b[1] * b[1] * (c[1] * c[1] + e[1] * e[1]) +
						(e[1] * d[1] - c[1] * f[1]) * (e[1] * d[1] - c[1] * f[1])) /
						(4.0 * (a[1] * a[1] + c[1] * c[1] + e[1] * e[1]));
					B1 = a[1] * (a[1] * seg.tau1 + b[1]) + c[1] * (c[1] * seg.tau1 + d[1]) + e[1] * (e[1] * seg.tau1 + f[1]) /
						(2.0 * sqrt(a[1] * a[1] + c[1] * c[1] + e[1] * e[1]));
					B2 = a[1] * (a[1] * seg.tau2 + b[1]) + c[1] * (c[1] * seg.tau2 + d[1]) + e[1] * (e[1] * seg.tau2 + f[1]) /
						(2.0 * sqrt(a[1] * a[1] + c[1] * c[1] + e[1] * e[1]));

					sum1 += sqrt(M_PI / (a[1] * a[1] + c[1] * c[1] + e[1] * e[1])) *
						(asinh(B2 / sqrt(A)) - asinh(B1 / sqrt(A)) + intExpansion(A, B2) - intExpansion(A, B1));

				}
			}
		}

		sum1 *= (props->visc * props->length / 8.0 / sqrt(M_PI) / M_PI / props->kx);		

		for (int m = 1; m <= props->M; m++)
		{
			a[0] = M_PI * (double)(m)* seg.r1.x / props->sizes.x;
			a[1] = M_PI * (double)(m)* seg.r2.x / props->sizes.x;
			b[0] = M_PI * (double)(m)* (props->r2.x - props->r1.x) / props->sizes.x;

			for (int n = 1; n <= props->N; n++)
			{
				c[0] = M_PI * (double)(n)* seg.r1.y / props->sizes.y;
				c[1] = M_PI * (double)(n)* seg.r2.y / props->sizes.y;
				d[0] = M_PI * (double)(n)* (props->r2.y - props->r1.y) / props->sizes.y;

				for (int l = 0; l <= props->L; l++)
				{
					e[0] = M_PI * (double)(l)* seg.r1.z / props->sizes.z;
					e[1] = M_PI * (double)(l)* seg.r2.z / props->sizes.z;
					f[0] = M_PI * (double)(l)* (props->r2.z - props->r1.z) / props->sizes.z;

					if (fabs(b[0]) + fabs(d[0]) + fabs(f[0]) > 0.0)
					{
						buf1 = (sin(a[0] - c[0] + e[0]) / (b[0] - d[0] + f[0]) +
							sin(a[0] - c[0] - e[0]) / (b[0] - d[0] - f[0]) -
							sin(a[0] + c[0] - e[0]) / (b[0] + d[0] - f[0]) -
							sin(a[0] + c[0] + e[0]) / (b[0] + d[0] + f[0]));

						buf2 = (sin(a[1] - c[1] + e[1]) / (b[0] - d[0] + f[0]) +
							sin(a[1] - c[1] - e[1]) / (b[0] - d[0] - f[0]) -
							sin(a[1] + c[1] - e[1]) / (b[0] + d[0] - f[0]) -
							sin(a[1] + c[1] + e[1]) / (b[0] + d[0] + f[0]));

						buf3 = ((double)(m) * (double)(m) / props->sizes.x / props->sizes.x +
							(double)(n) * (double)(n) / props->sizes.y / props->sizes.y +
							(double)(l) * (double)(l) / props->sizes.z / props->sizes.z);

						sum2 += exp(-M_PI * M_PI * buf3 * props->xi_c) / M_PI / M_PI / buf3 *
							sin(M_PI * (double)(m)* point.x / props->sizes.x) *
							sin(M_PI * (double)(n)* point.y / props->sizes.y) *
							cos(M_PI * (double)(l)* point.z / props->sizes.z) *
							(buf2 - buf1) / 4.0;
					}
					else
						std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAA!" << std::endl;
				}
			}
		}

		sum2 *= (8.0 * props->visc * props->length / props->sizes.x / props->sizes.y / props->sizes.z / props->kx);

		F3d[arr_idx] = sum1 + sum2;
	}
}