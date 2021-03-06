#ifndef __PERMMAT_HPP__
#define __PERMMAT_HPP__

#include <boost/numeric/ublas/fwd.hpp>
namespace boost{ namespace numeric{namespace ublas{
	template<class T, class A> class permutation_matrix;
}}}
#include <boost/shared_ptr.hpp>

/**
 * @brief PermMat Matrix -- describes an ordering of nodes
 * @author Hannes Schulz <mail at hannes-schulz dot de>
 * @date 2008-11-15
 */

class PermMat
{
  public:

    /**
     * Default constructor
     */
    PermMat();

    /**
     * Destructor
     */
    virtual ~PermMat();

	typedef boost::numeric::ublas::matrix<int> PermMatT;
	void setPermMat(boost::shared_ptr<PermMatT>);
	int getOriginalIndex(int i);
	boost::shared_ptr<PermMatT> getPermMat()const;

  private:

  	struct Impl;
	boost::shared_ptr<Impl> mImpl;

};

#endif /* #ifndef __PERMMAT_HPP__ */
