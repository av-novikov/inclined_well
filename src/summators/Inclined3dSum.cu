#include "src/summators/Inclined3dSum.cuh"
#include <cmath>

template <class T>
Inclined3dSum<T>::Inclined3dSum(const Parameters<T>* _props, const Well<T>* _well) : BaseSum<T>(_props, _well)
{
}

template <class T>
Inclined3dSum<T>::~Inclined3dSum()
{
}

template <class T>
void Inclined3dSum<T>::prepare()
{
	const int threadsPerBlock = 256;
	const int blocksPerGrid = 1 * size;

	// Allocate memory on device
	cudaErrorsChecker(cudaMalloc((void**)&F3d_dev, sizeof(T) * size));
	cudaErrorsChecker(cudaMalloc((void**)&F3d_buf, sizeof(T) * blocksPerGrid));
	cudaErrorsChecker(cudaMalloc((void**)&segs, sizeof(WellSegment<T>) * props->K));

	cudaErrorsChecker(cudaMemcpy(segs, well->segs, sizeof(WellSegment<T>) * props->K, cudaMemcpyHostToDevice));

	// Get device properties
	cudaDeviceProp deviceProp;
	cudaErrorsChecker(cudaGetDeviceProperties(&deviceProp, 0));

	dim3 blockSize(threadsPerBlock, 1, 1);			dim3 blockSizeR(1, 1, 1);
	dim3 gridSize(size, blocksPerGrid / size, 1);	dim3 gridSizeR(size, 1, 1);
	int sharedMem = threadsPerBlock * sizeof(T) + props->K * sizeof(WellSegment<T>);

	// Perform calculations
	prep<T, threadsPerBlock> << <gridSize, blockSize, sharedMem >> >(F3d_buf, *props, segs);

	cudaDeviceSynchronize();

	reduce3D<T> << <gridSizeR, blockSizeR >> >(F3d_buf, F3d_dev, blocksPerGrid / size);

	cudaErrorsChecker(cudaGetLastError());
	cudaDeviceSynchronize();

	// Transfer results on host memory
	cudaErrorsChecker(cudaMemcpy(F3d, F3d_dev, sizeof(T) * size, cudaMemcpyDeviceToHost));

	// Free device memory
	cudaErrorsChecker(cudaFree(F3d_dev));
	cudaErrorsChecker(cudaFree(F3d_buf));
	cudaErrorsChecker(cudaFree(segs));
}

template <class T>
T Inclined3dSum<T>::get2D(int seg_idx)
{
	/*T sum = 0.0;

	for(int k = 0; k < props->K; k++)
	{
	const WellSegment<T>& seg = well->segs[k];
	sum += F2d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (props->visc * props->sizes.x / CUDART_PI / CUDART_PI / props->sizes.z / props->kx / sin(props->alpha));

	return sum;*/
	return 0.0;
}

template <class T>
T Inclined3dSum<T>::get3D(int seg_idx)
{
	T sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment<T>& seg = well->segs[k];
		sum += F3d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (8.0 * props->visc * props->length / CUDART_PI / CUDART_PI / CUDART_PI /
		props->sizes.x / props->sizes.y / props->sizes.z / props->kx);

	return sum;
}

template class Inclined3dSum<float>;
template class Inclined3dSum<double>;
