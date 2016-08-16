#include "src/summators/InclinedSum.cuh"
#include <cmath>

template <class T>
InclinedSum<T>::InclinedSum(const Parameters<T>* _props, const Well<T>* _well) : BaseSum<T>(_props, _well)
{
}

template <class T>
InclinedSum<T>::~InclinedSum()
{
}

template <class T>
void InclinedSum<T>::prepare()
{
	const int threadsPerBlock = 128;
	const int blocksPerGrid = 16 * size;

	// Allocate memory on device
	cudaErrorsChecker(cudaMalloc((void**)&F2d_dev, sizeof(T) * size));
	cudaErrorsChecker(cudaMalloc((void**)&F2d_buf, sizeof(T) * blocksPerGrid));
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
	prep2D<T, threadsPerBlock><<<gridSize, blockSize, sharedMem>>>(F2d_buf, *props, segs);
	prep3D<T, threadsPerBlock><<<gridSize, blockSize, sharedMem>>>(F3d_buf, *props, segs);

	cudaDeviceSynchronize();

	reduce<T><<<gridSizeR, blockSizeR>>>(F2d_buf, F2d_dev, blocksPerGrid / size);
	reduce<T><<<gridSizeR, blockSizeR>>>(F3d_buf, F3d_dev, blocksPerGrid / size);

	cudaErrorsChecker(cudaGetLastError());
	cudaDeviceSynchronize();

	// Transfer results on host memory
	cudaErrorsChecker(cudaMemcpy(F2d, F2d_dev, sizeof(T) * size, cudaMemcpyDeviceToHost));
	cudaErrorsChecker(cudaMemcpy(F3d, F3d_dev, sizeof(T) * size, cudaMemcpyDeviceToHost));

	// Free device memory
	cudaErrorsChecker(cudaFree(F2d_dev));
	cudaErrorsChecker(cudaFree(F2d_buf));
	cudaErrorsChecker(cudaFree(F3d_dev));
	cudaErrorsChecker(cudaFree(F3d_buf));
	cudaErrorsChecker(cudaFree(segs));
}

template <class T>
T InclinedSum<T>::get2D(int seg_idx)
{
	T sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment<T>& seg = well->segs[k];
		sum += F2d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (props->visc * props->sizes.x / CUDART_PI / CUDART_PI / props->sizes.z / props->kx / sin(props->alpha));

	return sum;
}

template <class T>
T InclinedSum<T>::get3D(int seg_idx)
{
	T sum = 0.0;

	for (int k = 0; k < props->K; k++)
	{
		const WellSegment<T>& seg = well->segs[k];
		sum += F3d[seg_idx * props->K + k] * seg.rate / seg.length;
	}

	sum *= (2.0 * props->visc / CUDART_PI / props->sizes.x / props->sizes.z / props->kx / cos(props->alpha));

	return sum;
}

template class InclinedSum<float>;
template class InclinedSum<double>;
