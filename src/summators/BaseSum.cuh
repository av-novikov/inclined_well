#ifndef BASESUM_CUH_
#define BASESUM_CUH_

#include <cstdio>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include "math_constants.h"

#include "src/Well.cuh"

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

template <class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

template<>
struct SharedMemory<double>
{
    __device__ inline operator double *()
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }

    __device__ inline operator const double *() const
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }
};

template <class T>
class BaseSum
{
protected:

	const Parameters<T>* props;
	const Well<T>* well;

	int size;
	T* F2d;
	T* F3d;
	T* F2d_dev;
	T* F3d_dev;
	T* F2d_buf;
	T* F3d_buf;
	WellSegment<T>* segs;

public:
	BaseSum(const Parameters<T>* _props, const Well<T>* _well);
	virtual ~BaseSum();

	virtual void prepare() = 0;
	virtual T get2D(int seg_idx) = 0;
	virtual T get3D(int seg_idx) = 0;
	T getPres(int seg_idx);
};

#endif /* BASESUM_CUH_ */
