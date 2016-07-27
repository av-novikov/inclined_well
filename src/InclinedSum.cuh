#ifndef INCLINEDSUM_CUH_
#define INCLINEDSUM_CUH_

#include "src/Well.cuh"

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

template <typename cudaError>
inline void cudaErrorsChecker(cudaError error_t)
{
	unsigned int err = static_cast<unsigned int>(error_t);
	switch (err)
	{
	case 0:
		break;

	default:
		fprintf(stderr, "device error: %s\n", cudaGetErrorString(error_t));
		break;
	}
}

__global__ void prepare2D(double* F, const Parameters props, const WellSegment segs[]);
__global__ void prepare3D();

class InclinedSum
{
protected:

	const Parameters* props;
	const Well* well;

	int size;
	double* F2d;
	double* F3d;
	double* F2d_dev;
	double* F3d_dev;
	WellSegment* segs;

	void prepare();
	//__global__ void prepare2D(const Parameters props_reg);
	//__global__ void prepare3D();

public:
	InclinedSum(const Parameters* _props, const Well* well);
	~InclinedSum();
};

#endif /* INCLINEDSUM_CUH_ */