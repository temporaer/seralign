// 2008-11-15 Hannes Schulz <mail at hannes-schulz dot de>

#include <boost/numeric/ublas/matrix.hpp>
#include "SimMat.hpp"
#include <nana.h>
#undef C

/**********************************************************
 *          SimMat Implementation                *
 **********************************************************/
struct SimMat::Impl{
	Impl();
	~Impl();
	boost::shared_ptr<SimMat::SimMatT> mSimMat;
};

// Impl constructor
SimMat::Impl::Impl(){
}

// Impl destructor
SimMat::Impl::~Impl(){
}


/**********************************************************
 *          SimMat Interface                     *
 **********************************************************/
SimMat::SimMat()
	:mImpl(new Impl())
{
}

SimMat::~SimMat()
{
  // cleanup
}

boost::shared_ptr<SimMat::SimMatT> SimMat::getSimMat()
{
	return mImpl->mSimMat;
}

void SimMat::setSimMat(boost::shared_ptr<SimMat::SimMatT> ptr)
{
	mImpl->mSimMat = ptr;
}
