// 2008-11-15 Hannes Schulz <mail at hannes-schulz dot de>

#include <boost/numeric/ublas/matrix.hpp>
#include "AdjMat.hpp"
#include <nana.h>
#undef C

/**********************************************************
 *          AdjMat Implementation                *
 **********************************************************/
struct AdjMat::Impl{
	Impl();
	~Impl();
	boost::shared_ptr<AdjMat::AdjMatT> mAdjMat;
};

// Impl constructor
AdjMat::Impl::Impl(){
}

// Impl destructor
AdjMat::Impl::~Impl(){
}


/**********************************************************
 *          AdjMat Interface                     *
 **********************************************************/
AdjMat::AdjMat()
	:mImpl(new Impl())
{
}

AdjMat::~AdjMat()
{
  // cleanup
}
boost::shared_ptr<AdjMat::AdjMatT> AdjMat::getAdjMat()const
{
	return mImpl->mAdjMat;
}

void AdjMat::setAdjMat(boost::shared_ptr<AdjMat::AdjMatT> ptr)
{
	mImpl->mAdjMat = ptr;
}
