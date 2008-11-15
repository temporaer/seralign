/*       Created   :  10/06/2008 12:26:10 AM
 *       Last Change: Tue Oct 07 12:00 AM 2008 CEST
 */

#ifndef __SDP_PROB_HPP__
#define __SDP_PROB_HPP__

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
struct SDPProb{
	typedef boost::numeric::ublas::matrix<double,
			boost::numeric::ublas::column_major>  MatT;
	typedef std::vector<MatT>                     MatVecT;
	typedef boost::numeric::ublas::vector<double> VecT;
	MatT    C;
	MatVecT F;
	VecT    b;
};

#endif /* __SDP_PROB_HPP__ */
