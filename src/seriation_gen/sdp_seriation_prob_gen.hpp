/*       Created   :  10/06/2008 12:23:53 AM
 *       Last Change: Sat Nov 15 06:00 PM 2008 CET
 */

#ifndef __SDP_SERIATION_PROB_GEN_HPP__
#define __SDP_SERIATION_PROB_GEN_HPP__
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/shared_ptr.hpp>
#include <AdjMat.hpp>
class SDPProb;
class SDPSeriationProbGen{
	public:
		typedef AdjMat::AdjMatT AdjMatT;
		SDPSeriationProbGen(boost::shared_ptr<AdjMatT>);
		void operator()(SDPProb&);

    public:
		boost::shared_ptr<AdjMatT>  mAdjMat;
};
#endif /* __SDP_SERIATION_PROB_GEN_HPP__ */

