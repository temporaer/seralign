/*       Created   :  10/06/2008 12:54:34 AM
 *       Last Change: Sat Dec 06 01:00 PM 2008 CET
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

		static inline void setVerbose(bool b){mVerbose=b;}
		static inline bool  isVerbose()      {return mVerbose;}
	private:
		static bool mVerbose;
};
#endif /* __SDP_WRAPPER_HPP__ */
