#ifndef __SIMMAT_HPP__
#define __SIMMAT_HPP__

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/shared_ptr.hpp>
#include <nana.h>
#undef C

/**
 * @brief SimMat Similarity Matrix -- describes similarity between nodes
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 * for two nodes i and j in a graph, represent their similarity in SimMat(i,j)
 */

class SimMat
{
  public:

    /**
     * Default constructor
     */
    SimMat();

    /**
     * Destructor
     */
    virtual ~SimMat();

	typedef boost::numeric::ublas::matrix<double> SimMatT;
	void setSimMat(boost::shared_ptr<SimMatT>);
	boost::shared_ptr<SimMatT> getSimMat();

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __SIMMAT_HPP__ */
