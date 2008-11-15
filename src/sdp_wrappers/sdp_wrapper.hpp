/*       Created   :  10/06/2008 12:54:34 AM
 *       Last Change: Fri Oct 31 10:00 AM 2008 CET
 */

#ifndef __SDP_WRAPPER_HPP__
#define __SDP_WRAPPER_HPP__

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
class SDPProb;

class SDPWrapper{
	public:
		typedef boost::numeric::ublas::matrix<double,
				boost::numeric::ublas::row_major> AnswerT;
		virtual AnswerT operator()(const SDPProb&);
		virtual void configure();
		virtual ~SDPWrapper();
};
#endif /* __SDP_WRAPPER_HPP__ */
