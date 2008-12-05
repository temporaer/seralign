#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "jumbled_adjmat_gen.hpp"
#include <factory/factory.h>
using namespace std;
using namespace boost::numeric::ublas;

JumbledAdjMatGen::JumbledAdjMatGen(ProbAdjPerm& pap)
	: mPap(pap)
{
}

ProbAdjPerm JumbledAdjMatGen::operator()()
{
	int n = mPap.getAdjMat()->size1();
	PermMat::PermMatT perm(zero_matrix<int>(n,n));
	std::vector<int> order(n,0);
	for(int i=0;i<n;i++)
		order[i] = i;
	random_shuffle(order.begin(),order.end());
	for(int i=0;i<n;i++)
		perm(i, order[i]) = 1;
	matrix<double> tmp(n,n);
	noalias(tmp) = prod(perm, *mPap.getAdjMat());
	(*mPap.getAdjMat())  = prod(tmp, trans(perm));
	
	// save perm mat in problem
	if(mPap.getPermMat().get() == NULL){
		perm = trans(perm);
		mPap.setPermMat(boost::shared_ptr<PermMat::PermMatT>(new PermMat::PermMatT(perm)));
	}
	else{
		PermMat::PermMatT tmp(n,n);
		noalias(tmp) = prod(trans(perm), *mPap.getPermMat());
		*mPap.getPermMat() = prod(trans(perm),*mPap.getPermMat() );
	}

	noalias(tmp) = prod(perm, *mPap.getPermMat());
	(*mPap.getPermMat()) = prod(tmp, trans(perm));
	return mPap;
}
JumbledAdjMatGen::~JumbledAdjMatGen()
{
}

