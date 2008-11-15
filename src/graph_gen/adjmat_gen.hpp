#ifndef __ADJMAT_GEN_HPP__
#define __ADJMAT_GEN_HPP__

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/fwd.hpp>

class AdjMatGen{
	public:
		virtual void configure();
		typedef boost::numeric::ublas::matrix<double,
				boost::numeric::ublas::column_major
			> AdjMatT;
		virtual boost::shared_ptr<AdjMatT> operator()();
		virtual ~AdjMatGen();
};

#endif /* __ADJMAT_GEN_HPP__ */

