#ifndef WELLFLOW_HPP_
#define WELLFLOW_HPP_

#include <string>
#include <gsl/gsl_linalg.h>

#include "src/inclined_sum/InclinedSum.hpp"
#include "src/Well.hpp"

class WellFlow
{
	protected:
		Parameters props;
		Well* well;
		
		void loadTask(const std::string fileName);
		
		gsl_vector* q_gsl;
		gsl_vector* dq_gsl;
		gsl_vector* b_gsl;
		gsl_vector* x_gsl;
		gsl_matrix* a_gsl;
		gsl_matrix* dpdq_gsl;
		gsl_permutation* perm_gsl;
		
		/*double highIntegral2D(const Point& r, double xi_c);
		double lowIntegral2D(const Point& r, double xi_c);
		double calc3D(const Point& r);*/
		
	public:
		WellFlow(const std::string fileName);
		~WellFlow();
		
		InclinedSum* inclSum;
		
		void setSummator(InclinedSum* _inclSum);
		const Parameters* getProps() const;
		const Well* getWell() const;
		void findRateDistribution();
				
		void calcPressure();
		double getP_bhp();
		Point getObsPoint() const;
};

#endif /* WELLFLOW_HPP_ */
