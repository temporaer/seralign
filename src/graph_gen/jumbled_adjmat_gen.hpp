#ifndef __JUMBLED_ADJMAT_GEN_HPP__
#define __JUMBLED_ADJMAT_GEN_HPP__
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include "adjmat_gen.hpp"


class JumbledAdjMatGen : public AdjMatGen{
	private:
		ProbAdjPerm& mPap;
	public:
		JumbledAdjMatGen(ProbAdjPerm& pap);

		virtual ProbAdjPerm operator()();
		virtual ~JumbledAdjMatGen();
};

#endif /* __JUMBLED_ADJMAT_GEN_HPP__ */
