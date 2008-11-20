#ifndef __ADJMAT_GEN_HPP__
#define __ADJMAT_GEN_HPP__

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <ProbData.hpp>

class AdjMatGen{
	public:
		virtual void configure();
		virtual ProbAdjPerm operator()();
		virtual ~AdjMatGen();
		virtual bool hasNext();
		virtual std::string getPrologDescription(int idx);
		virtual std::string getGraphID();
};

#endif /* __ADJMAT_GEN_HPP__ */

