#ifndef __LAPLACIAN_HPP__
#define __LAPLACIAN_HPP__

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/shared_ptr.hpp>

/**
 * @brief Laplacian Matrix -- combines Adjacency matrix and Degree of nodes in graph
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 */

class Laplacian
{
  public:

    /**
     * Default constructor
     */
    Laplacian();

    /**
     * Destructor
     */
    virtual ~Laplacian();

	typedef boost::numeric::ublas::matrix<double> LaplacianT;
	void setLaplacian(boost::shared_ptr<LaplacianT>);
	boost::shared_ptr<LaplacianT> getLaplacian();

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __LAPLACIAN_HPP__ */
