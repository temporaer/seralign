#ifndef __ADJMAT_HPP__
#define __ADJMAT_HPP__

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/shared_ptr.hpp>
#include <nana.h>
#undef C

/**
 * @brief AdjMat contains an Adjacency Matrix
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 * This class represents an Adjacency matrix.
 */

class AdjMat
{
  public:

    /**
     * Default constructor
     */
    AdjMat();

    /**
     * Destructor
     */
    virtual ~AdjMat();

	typedef boost::numeric::ublas::matrix<double> AdjMatT;
	void setAdjMat(boost::shared_ptr<AdjMatT>);
	boost::shared_ptr<AdjMatT> getAdjMat();

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __ADJMAT_HPP__ */
