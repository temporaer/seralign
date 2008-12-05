// 2008-11-15 Hannes Schulz <mail at hannes-schulz dot de>

#include <stdexcept>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "PermMat.hpp"
#include <nana.h>
#undef C

namespace ublas = boost::numeric::ublas;
using namespace std;

/**********************************************************
 *          PermMat Implementation                *
 **********************************************************/
struct PermMat::Impl{
	Impl();
	~Impl();
	boost::shared_ptr<PermMat::PermMatT> mPermMat;
	int getOriginalIndex(int i);
};

// Impl constructor
PermMat::Impl::Impl(){
}

// Impl destructor
PermMat::Impl::~Impl(){
}

int PermMat::Impl::getOriginalIndex(int i)
{
	ublas::vector<int> v = ublas::column(*mPermMat,i);
	ublas::vector<int>::iterator it = find(v.begin(),v.end(),1);
	if(it==v.end())
		throw runtime_error("PermMat::getOriginalIndex(): Could not map index!");
	return distance(v.begin(),it);
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

boost::shared_ptr<PermMat::PermMatT> PermMat::getPermMat()const
{
	return mImpl->mPermMat;
}

void PermMat::setPermMat(boost::shared_ptr<PermMat::PermMatT> ptr)
{
	mImpl->mPermMat = ptr;
}
int PermMat::getOriginalIndex(int i)
{
	return mImpl->getOriginalIndex(i);
}

