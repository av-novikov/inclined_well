#include "src/Well.cuh"

#include <iostream>
#include <functional>
#include <fstream>

using std::cout;
using std::endl;
using std::ofstream;

template <class T>
Well<T>::Well(const Point<T>& _r1, const Point<T>& _r2, const int _num, const T _r_w) : r1(_r1), r2(_r2), num(_num), r_w(_r_w)
{
	length = sqrt((r2 - r1) * (r2 - r1));
	alpha = atan((r2.x - r1.x) / (r1.z - r2.z));
	pres_av = pres_dev = rate = 0.0;

	segs = new WellSegment<T>[num];

	Point<T> tmp1 = r1;
	Point<T> tmp2, tmp3;
	for (int i = 0; i < num; i++)
	{
		tmp2 = r1 + (T)(i + 1) * (r2 - r1) / (T)(num);
		tmp3 = (tmp1 + tmp2) / (T)(2.0);		tmp3.y += r_w; // cos(alpha);
		segs[i].r1 = tmp1;
		segs[i].r2 = tmp2;
		segs[i].r_bhp = tmp3;
		segs[i].length = sqrt((tmp2 - tmp1) * (tmp2 - tmp1));
		segs[i].pres = segs[i].pres2D = segs[i].pres3D = segs[i].rate = 0.0;

		tmp1 = tmp2;
	};
}

template <class T>
Well<T>::~Well()
{
}

template <class T>
void Well<T>::setRate(T _rate)
{
	rate = _rate;
}

template <class T>
void Well<T>::setUniformRate()
{
	for (int i = 0; i < num; i++)
	{
		WellSegment<T>& seg = segs[i];
		seg.rate = seg.length / length * rate;
	}
}

template <class T>
void Well<T>::printRates(const Parameters<T>* props)
{
	T av2d = 0.0, av3d = 0.0;
	for (int i = 0; i < num; i++)
	{
		av2d += segs[i].pres2D;
		av3d += segs[i].pres3D;
		cout << "--- " << i << " ---\tRate = " << segs[i].rate * 86400.0 * props->x_dim * props->x_dim * props->x_dim / props->t_dim <<
			"\tPressure = " << segs[i].pres * props->p_dim / (T)(BAR) << "\t" << segs[i].pres2D * props->p_dim / (T)(BAR) <<
			"\t" << segs[i].pres3D * props->p_dim / (T)(BAR) << endl;
	}
	av2d /= (T)(num);	av3d /= (T)(num);

	cout << "Av. pressure = " << pres_av * props->p_dim / BAR << endl;
	cout << "Av. 2D = " << av2d * props->p_dim / BAR << endl;
	cout << "Av. 3D = " << av3d * props->p_dim / BAR << endl;
	cout << "Deviation = " << pres_dev * props->p_dim * props->p_dim / BAR / BAR << endl;
}

template <class T>
void Well<T>::writeRates(const Parameters<T>* props)
{
	T av2d = 0.0, av3d = 0.0;
	ofstream file;
	file.open("rate.dat", ofstream::out);
	for (int i = 0; i < num; i++)
	{
		av2d += segs[i].pres2D;
		av3d += segs[i].pres3D;

		file << i << "\t" <<
			segs[i].rate * 86400.0 * props->x_dim * props->x_dim * props->x_dim / props->t_dim << "\t" <<
			segs[i].pres * props->p_dim / (T)(BAR) <<
			"\t" << segs[i].pres2D * props->p_dim / (T)(BAR) <<
			"\t" << segs[i].pres3D * props->p_dim / (T)(BAR) << endl;
	}

	av2d /= (T)(num);	av3d /= (T)(num);

	file << "Av. pressure = " << pres_av * props->p_dim / BAR << endl;
	file << "Av. 2D = " << av2d * props->p_dim / BAR << endl;
	file << "Av. 3D = " << av3d * props->p_dim / BAR << endl;
	file << "Deviation = " << pres_dev * props->p_dim * props->p_dim / BAR / BAR << endl;

	file.close();
}

template <class T>
Well<T>& Well<T>::operator=(const Well<T>& well)
{
	for (int i = 0; i < num; i++)
		segs[i] = well.segs[i];

	return *this;
};

template struct Point<float>;
template struct Point<double>;
template struct WellSegment<float>;
template struct WellSegment<double>;
template struct Parameters<float>;
template struct Parameters<double>;
template class Well<float>;
template class Well<double>;