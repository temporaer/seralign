// 2008-11-15 Hannes Schulz <mail at hannes-schulz dot de>

#include <boost/numeric/ublas/matrix.hpp>
#include "Laplacian.hpp"
#include <nana.h>
#undef C

/**********************************************************
 *          Laplacian Implementation                *
 **********************************************************/
struct Laplacian::Impl{
	Impl();
	~Impl();
	boost::shared_ptr<Laplacian::LaplacianT> mLaplacian;
};

// Impl constructor
Laplacian::Impl::Impl(){
}

// Impl destructor
Laplacian::Impl::~Impl(){
}


/**********************************************************
 *          Laplacian Interface                     *
 **********************************************************/
Laplacian::Laplacian()
	:mImpl(new Impl())
{
}

Laplacian::~Laplacian()
{
  // cleanup
}

boost::shared_ptr<Laplacian::LaplacianT> Laplacian::getLaplacian()
{
	return mImpl->mLaplacian;
}

void Laplacian::setLaplacian(boost::shared_ptr<Laplacian::LaplacianT> ptr)
{
	mImpl->mLaplacian = ptr;
}
