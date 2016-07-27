#include <cstdio>

#include "src/InclinedSum.cuh"
#include "math_constants.h"

__global__ void prepare2D(double* F, const Parameters props, const WellSegment* segs)
{
	//int idx = blockDim.x * blockIdx.x + threadIdx.x;
	
	double buf;
	double sum_prev = 0.0;
	int break_idx = 0;
	int k;

	for (int arr_idx = 0; arr_idx < props.K * props.K; arr_idx++)
	{
		k = arr_idx % props.K;
		const Point& r = segs[int((double)(arr_idx) / (double)(props.K))].r_bhp;

		F[arr_idx] = sum_prev = 0.0;

		break_idx = 0;

		for (int m = 1; m <= props.M; m++)
		{
			buf = sin(CUDART_PI * (double)(m)* r.x / props.sizes.x) *
				(cos(CUDART_PI * (double)(m) * segs[k].r1.x / props.sizes.x) -
					cos(CUDART_PI * (double)(m) * segs[k].r2.x / props.sizes.x));

			for (int i = -props.I; i <= props.I; i++)
			{
				F[arr_idx] += buf / segs[k].length / (double)(m) / (double)(m)*
					(exp(-CUDART_PI * (double)(m)* fabs(r.y - props.r1.y + 2.0 * (double)(i)* props.sizes.y) / props.sizes.x) -
						exp(-CUDART_PI * (double)(m)* fabs(r.y + props.r1.y + 2.0 * (double)(i)* props.sizes.y) / props.sizes.x));
			}

			if (fabs(F[arr_idx] - sum_prev) > F[arr_idx] * EQUALITY_TOLERANCE)
			{
				sum_prev = F[arr_idx];
				break_idx = 0;
			}
			else
				break_idx++;

			if (break_idx > 1)
				break;
		}
	}
}

__global__ void prepare3D()
{
}

InclinedSum::InclinedSum(const Parameters* _props, const Well* _well) : props(_props), well(_well)
{
	size = props->K * props->K;
	F2d = new double [ size ];
	F3d = new double [ size ];

	prepare();
	//prepare2D();
	//prepare3D();
}

InclinedSum::~InclinedSum()
{
	delete [] F2d;
	delete [] F3d;
}

void InclinedSum::prepare()
{
	// Allocate memory on device
	cudaErrorsChecker(cudaMalloc((void**)&F2d_dev, sizeof(double) * size));
	cudaErrorsChecker(cudaMalloc((void**)&F3d_dev, sizeof(double) * size));
	cudaErrorsChecker(cudaMalloc((void**)&segs, sizeof(WellSegment) * props->K));
	
	cudaErrorsChecker(cudaMemcpy(segs, well->segs, sizeof(WellSegment) * props->K, cudaMemcpyHostToDevice));

	// Get device properties
	cudaDeviceProp deviceProp;
	cudaErrorsChecker(cudaGetDeviceProperties(&deviceProp, 0));

	dim3 blockSize(size, 1, 1);
	dim3 gridSize(1, 1, 1);

	// Perform calculations
	prepare2D<<<gridSize, blockSize>>>(F2d_dev, *props, segs);
	cudaErrorsChecker(cudaGetLastError());
	cudaDeviceSynchronize();

	// Transfer results on host memody
	cudaErrorsChecker(cudaMemcpy(F2d, F2d_dev, sizeof(double) * size, cudaMemcpyDeviceToHost));
	cudaErrorsChecker(cudaMemcpy(F3d, F3d_dev, sizeof(double) * size, cudaMemcpyDeviceToHost));

	// Free device memory
	cudaErrorsChecker(cudaFree(F2d_dev));
	cudaErrorsChecker(cudaFree(F3d_dev));
	cudaErrorsChecker(cudaFree(segs));
}

