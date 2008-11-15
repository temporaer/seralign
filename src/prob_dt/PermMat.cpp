// 2008-11-15 Hannes Schulz <mail at hannes-schulz dot de>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "PermMat.hpp"

/**********************************************************
 *          PermMat Implementation                *
 **********************************************************/
struct PermMat::Impl{
	Impl();
	~Impl();
	boost::shared_ptr<PermMat::PermMatT> mPermMat;
};

// Impl constructor
PermMat::Impl::Impl(){
}

// Impl destructor
PermMat::Impl::~Impl(){
}


/**********************************************************
 *          PermMat Interface                     *
 **********************************************************/
PermMat::PermMat()
	:mImpl(new Impl())
{
}

PermMat::~PermMat()
{
  // cleanup
}

boost::shared_ptr<PermMat::PermMatT> PermMat::getPermMat()
{
	return mImpl->mPermMat;
}

void PermMat::setPermMat(boost::shared_ptr<PermMat::PermMatT> ptr)
{
	mImpl->mPermMat = ptr;
}
