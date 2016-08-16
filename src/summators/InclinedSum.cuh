#ifndef INCLINEDSUM_CUH_
#define INCLINEDSUM_CUH_

#include "device_functions.h"

#include "src/summators/BaseSum.cuh"

template <class T, int blockSize>
__global__ void prep2D(T* Fbuf, const Parameters<T> props, const WellSegment<T>* _segs)
{
	if (blockIdx.x < gridDim.x && blockIdx.y < gridDim.y)
	{
		const int tid = threadIdx.x;

		// Allocate an array to fill
		T* sdata = SharedMemory<T>();
		T* sdata_cur = &sdata[threadIdx.x];

		// Allocate an array for well segments
		WellSegment<T>* segs = (WellSegment<T>*)(&sdata[blockSize]);
		if (tid < props.K)
			segs[tid] = _segs[tid];

		__syncthreads();

		const int itersPerBlock = int(props.M / gridDim.y);
		const int startIdx = 1 + itersPerBlock * blockIdx.y + int(itersPerBlock / blockSize * tid);
		const int finishIdx = itersPerBlock * blockIdx.y + int(itersPerBlock / blockSize * (tid + 1));
		const int k = blockIdx.x % props.K;
		const Point<T>& r = segs[int((T)(blockIdx.x) / (T)(props.K))].r_bhp;
		T buf;

		// Perform calculations
		*sdata_cur = 0.0;
		for (int m = startIdx; m <= finishIdx; m++)
		{
			buf = sin(CUDART_PI * (T)(m)* r.x / props.sizes.x) *
				(cos(CUDART_PI * (T)(m)* segs[k].r1.x / props.sizes.x) -
					cos(CUDART_PI * (T)(m)* segs[k].r2.x / props.sizes.x));

			for (int i = -props.I; i <= props.I; i++)
			{
				*sdata_cur += buf / (T)(m) / (T)(m)*
					(exp(-CUDART_PI * (T)(m)* fabs(r.y - props.r1.y + 2.0 * (T)(i)* props.sizes.y) / props.sizes.x) -
						exp(-CUDART_PI * (T)(m)* fabs(r.y + props.r1.y + 2.0 * (T)(i)* props.sizes.y) / props.sizes.x));
			}
		}

		__syncthreads();

		// Reduction in block
		if (blockSize >= 512) { if (tid < 256) { *sdata_cur += sdata[tid + 256]; } __syncthreads(); }
		if (blockSize >= 256) { if (tid < 128) { *sdata_cur += sdata[tid + 128]; } __syncthreads(); }
		if (blockSize >= 128) { if (tid < 64) { *sdata_cur += sdata[tid + 64]; } __syncthreads(); }
		if (blockSize >= 64) { if (tid < 32) { *sdata_cur += sdata[tid + 32]; } __syncthreads(); }
		if (blockSize >= 32) { if (tid < 16) { *sdata_cur += sdata[tid + 16]; } __syncthreads(); }
		if (blockSize >= 16) { if (tid < 8) { *sdata_cur += sdata[tid + 8]; } __syncthreads(); }
		if (blockSize >= 8) { if (tid < 4) { *sdata_cur += sdata[tid + 4]; } __syncthreads(); }
		if (blockSize >= 4) { if (tid < 2) { *sdata_cur += sdata[tid + 2]; } __syncthreads(); }
		if (blockSize >= 2) { if (tid < 1) { *sdata_cur += sdata[tid + 1]; } __syncthreads(); }

		if (tid == 0)
			Fbuf[blockIdx.x * gridDim.y + blockIdx.y] = sdata[0];
	}
};

template <class T, int blockSize>
__global__ void prep3D(T* Fbuf, const Parameters<T> props, const WellSegment<T>* _segs)
{
	if (blockIdx.x < gridDim.x && blockIdx.y < gridDim.y)
	{
		const int tid = threadIdx.x;

		// Allocate an array to fill
		T* sdata = SharedMemory<T>();
		T* sdata_cur = &sdata[threadIdx.x];

		// Allocate an array for well segments
		WellSegment<T>* segs = (WellSegment<T>*)(&sdata[blockDim.x]);
		if (threadIdx.x < props.K)
			segs[threadIdx.x] = _segs[threadIdx.x];

		__syncthreads();

		const int itersPerBlock = int(props.M / gridDim.y);
		const int startIdx = 1 + itersPerBlock * blockIdx.y + int(itersPerBlock / blockSize * tid);
		const int finishIdx = itersPerBlock * blockIdx.y + int(itersPerBlock / blockSize * (tid + 1));
		const int k = blockIdx.x % props.K;
		const Point<T>& r = segs[int((T)(blockIdx.x) / (T)(props.K))].r_bhp;
		T buf1, buf2;

		// Perform calculations
		*sdata_cur = 0.0;
		for (int m = startIdx; m <= finishIdx; m++)
		{
			for (int l = 1; l <= props.L; l++)
			{
				buf1 = sin(CUDART_PI * (T)(m)* r.x / props.sizes.x) * cos(CUDART_PI * (T)(l)* r.z / props.sizes.z) *
					((cos(CUDART_PI * (T)(m)* segs[k].r1.x / props.sizes.x - CUDART_PI * (T)(l)* segs[k].r1.z / props.sizes.z) -
						cos(CUDART_PI * (T)(m)* segs[k].r2.x / props.sizes.x - CUDART_PI * (T)(l)* segs[k].r2.z / props.sizes.z)) /
						(CUDART_PI * (T)(m)* tan(props.alpha) / props.sizes.x + CUDART_PI * (T)(l) / props.sizes.z) +
						(cos(CUDART_PI * (T)(m)* segs[k].r1.x / props.sizes.x + CUDART_PI * (T)(l)* segs[k].r1.z / props.sizes.z) -
							cos(CUDART_PI * (T)(m)* segs[k].r2.x / props.sizes.x + CUDART_PI * (T)(l)* segs[k].r2.z / props.sizes.z)) /
						(CUDART_PI * (T)(m)* tan(props.alpha) / props.sizes.x - CUDART_PI * (T)(l) / props.sizes.z)
						) / 2.0;

				buf2 = sqrt((T)(m)* (T)(m) / props.sizes.x / props.sizes.x + (T)(l)* (T)(l) / props.sizes.z / props.sizes.z);

				for (int i = -props.I; i <= props.I; i++)
				{
					*sdata_cur += buf1 / buf2 *
						(exp(-CUDART_PI * buf2 * fabs(r.y - props.r1.y + 2.0 * (T)(i)* props.sizes.y)) -
							exp(-CUDART_PI * buf2 * fabs(r.y + props.r1.y + 2.0 * (T)(i)* props.sizes.y)));
				}
			}
		}

		__syncthreads();

		// Reduction in block
		if (blockSize >= 512) { if (tid < 256) { *sdata_cur += sdata[tid + 256]; } __syncthreads(); }
		if (blockSize >= 256) { if (tid < 128) { *sdata_cur += sdata[tid + 128]; } __syncthreads(); }
		if (blockSize >= 128) { if (tid < 64) { *sdata_cur += sdata[tid + 64]; } __syncthreads(); }
		if (blockSize >= 64) { if (tid < 32) { *sdata_cur += sdata[tid + 32]; } __syncthreads(); }
		if (blockSize >= 32) { if (tid < 16) { *sdata_cur += sdata[tid + 16]; } __syncthreads(); }
		if (blockSize >= 16) { if (tid < 8) { *sdata_cur += sdata[tid + 8]; } __syncthreads(); }
		if (blockSize >= 8) { if (tid < 4) { *sdata_cur += sdata[tid + 4]; } __syncthreads(); }
		if (blockSize >= 4) { if (tid < 2) { *sdata_cur += sdata[tid + 2]; } __syncthreads(); }
		if (blockSize >= 2) { if (tid < 1) { *sdata_cur += sdata[tid + 1]; } __syncthreads(); }

		if (tid == 0)
			Fbuf[blockIdx.x * gridDim.y + blockIdx.y] = sdata[0];
	}
};

template <class T>
__global__ void reduce(T* Fbuf, T* F, int blocksPerNum)
{
	if (blocksPerNum > 1)
	{
		const int startIdx = blockIdx.x * blocksPerNum;
		const int finishIdx = (blockIdx.x + 1) * blocksPerNum - 1;
		//printf("%d\t%d\t%d\n", threadIdx.x, startIdx, finishIdx);
		T sum = 0.0;

		for (int i = finishIdx; i >= startIdx; --i)
		{
			sum += Fbuf[i];
		}

		F[blockIdx.x] = sum;
	}
	else {
		F[blockIdx.x] = Fbuf[blockIdx.x * blocksPerNum];
	}

};

template <class T>
class InclinedSum : public BaseSum<T>
{
protected:
	using BaseSum<T>::props;
	using BaseSum<T>::well;
	
	using BaseSum<T>::size;
	using BaseSum<T>::F2d;
	using BaseSum<T>::F3d;
	using BaseSum<T>::F2d_dev;
	using BaseSum<T>::F3d_dev;
	using BaseSum<T>::F2d_buf;
	using BaseSum<T>::F3d_buf;
	using BaseSum<T>::segs;

public:
	InclinedSum(const Parameters<T>* _props, const Well<T>* well);
	~InclinedSum();
	
	void prepare();
	T get2D(int seg_idx);
	T get3D(int seg_idx);
};

#endif /* INCLINEDSUM_CUH_ */
