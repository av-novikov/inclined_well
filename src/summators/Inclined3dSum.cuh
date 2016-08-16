#ifndef INCLINED3DSUM_CUH_
#define iNCLINED3DSUM_CUH_

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_functions.h>

#include "src/summators/BaseSum.cuh"

template <class T>
struct Boof
{
	T x;
	T x1;
	T x2;
	T y;
	T y1;
	T y2;
	T z;
	T z1;
	T z2;
};

template <class T, int blockSize>
__global__ void prep(T* Fbuf, const Parameters<T> props, const WellSegment<T>* _segs)
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


		const int lastIdx = props.M * props.N * (props.L + 1);
		const int itersPerBlock = (T)(lastIdx) / (T)(gridDim.y);
		const int startIdx = itersPerBlock * blockIdx.y + (T)(itersPerBlock) / (T)(blockSize)* tid;
		const int finishIdx = itersPerBlock * blockIdx.y + (T)(itersPerBlock) / (T)(blockSize)* (tid + 1);
		const int k = blockIdx.x % props.K;
		const Point<T>& r = segs[int((T)(blockIdx.x) / (T)(props.K))].r_bhp;
		WellSegment<T>* seg = &segs[k];
		Boof<T> buf;
		int l = startIdx % (props.L + 1);
		int n = 1 + (((startIdx - l) / (props.L + 1)) % props.N);
		int m = 1 + startIdx / (props.N * (props.L + 1));

		buf.x = (T)(m)* (props.r2.x - props.r1.x) / props.sizes.x;
		buf.x1 = (T)(m)* seg->r1.x / props.sizes.x * CUDART_PI;
		buf.x2 = (T)(m)* seg->r2.x / props.sizes.x * CUDART_PI;

		buf.y = (T)(n)* (props.r2.y - props.r1.y) / props.sizes.y;
		buf.y1 = (T)(n)* seg->r1.y / props.sizes.y * CUDART_PI;
		buf.y2 = (T)(n)* seg->r2.y / props.sizes.y * CUDART_PI;

		// Perform calculations
		*sdata_cur = 0.0;
		for (int i = startIdx; i < finishIdx; i++)
		{
			if (l == props.L + 1)
			{
				l = 0;
				n++;

				if (n == props.N + 1)
				{
					n = 1;
					m++;

					buf.x = (T)(m)* (props.r2.x - props.r1.x) / props.sizes.x;
					buf.x1 = (T)(m)* seg->r1.x / props.sizes.x * CUDART_PI;
					buf.x2 = (T)(m)* seg->r2.x / props.sizes.x * CUDART_PI;
				}

				buf.y = (T)(n)* (props.r2.y - props.r1.y) / props.sizes.y;
				buf.y1 = (T)(n)* seg->r1.y / props.sizes.y * CUDART_PI;
				buf.y2 = (T)(n)* seg->r2.y / props.sizes.y * CUDART_PI;
			}

			buf.z = (T)(l)* (props.r2.z - props.r1.z) / props.sizes.z;
			buf.z1 = (T)(l)* seg->r1.z / props.sizes.z * CUDART_PI;
			buf.z2 = (T)(l)* seg->r2.z / props.sizes.z * CUDART_PI;

			*sdata_cur += ((sin(buf.x2 - buf.y2 + buf.z2) - sin(buf.x1 - buf.y1 + buf.z1)) /
				(buf.x - buf.y + buf.z) +
				(sin(buf.x2 - buf.y2 - buf.z2) - sin(buf.x1 - buf.y1 - buf.z1)) /
				(buf.x - buf.y - buf.z) -
				(sin(buf.x2 + buf.y2 - buf.z2) - sin(buf.x1 + buf.y1 - buf.z1)) /
				(buf.x + buf.y - buf.z) -
				(sin(buf.x2 + buf.y2 + buf.z2) - sin(buf.x1 + buf.y1 + buf.z1)) /
				(buf.x + buf.y + buf.z)) / 4.0 *
				sin(CUDART_PI * (T)(m)* r.x / props.sizes.x) *
				sin(CUDART_PI * (T)(n)* r.y / props.sizes.y) *
				cos(CUDART_PI * (T)(l)* r.z / props.sizes.z) /
				((T)(m)* (T)(m) / props.sizes.x / props.sizes.x +
					(T)(n)* (T)(n) / props.sizes.y / props.sizes.y +
					(T)(l)* (T)(l) / props.sizes.z / props.sizes.z);

			l++;
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
__global__ void reduce3D(T* Fbuf, T* F, int blocksPerNum)
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
class Inclined3dSum : public BaseSum<T>
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
	Inclined3dSum(const Parameters<T>* _props, const Well<T>* well);
	~Inclined3dSum();

	void prepare();
	T get2D(int seg_idx);
	T get3D(int seg_idx);
};

#endif /* INCLINED3DSUM_CUH_ */
