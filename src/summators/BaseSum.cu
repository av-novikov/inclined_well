#include "src/summators/BaseSum.cuh"

template <class T>
BaseSum<T>::BaseSum(const Parameters<T>* _props, const Well<T>* _well) : props(_props), well(_well)
{
	size = props->K * props->K;
	F2d = new T [ size ];
	F3d = new T [ size ];
}

template <class T>
BaseSum<T>::~BaseSum()
{
	delete [] F2d;
	delete [] F3d;
}

template <class T>
T BaseSum<T>::getPres(int seg_idx)
{
	return get2D(seg_idx) + get3D(seg_idx);
}

template class BaseSum<float>;
template class BaseSum<double>;
