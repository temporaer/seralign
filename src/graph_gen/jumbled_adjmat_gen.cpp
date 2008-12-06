#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "jumbled_adjmat_gen.hpp"
#include <factory/factory.h>
#include <nana.h>
using namespace std;
using namespace boost::numeric::ublas;

JumbledAdjMatGen::JumbledAdjMatGen(ProbAdjPerm& pap)
	: mPap(pap)
{
}

ProbAdjPerm JumbledAdjMatGen::operator()()
{
	ProbAdjPerm pap;
	int n = mPap.getAdjMat()->size1();
	boost::shared_ptr<PermMat::PermMatT> perm(new PermMat::PermMatT(zero_matrix<int>(n,n)));
	std::vector<int> order(n,0);
	for(int i=0;i<n;i++)
		order[i] = i;
	random_shuffle(order.begin(),order.end());
	for(int i=0;i<n;i++)
		(*perm)(i, order[i]) = 1;
	matrix<double> tmp(n,n);
	noalias(tmp) = prod(*perm, *mPap.getAdjMat());
	pap.setAdjMat( boost::shared_ptr<AdjMat::AdjMatT>(new AdjMat::AdjMatT(prod(tmp, trans(*perm)))));
	
	// save perm mat in problem
	if(mPap.getPermMat().get() == NULL){
		*perm = trans(*perm);
		pap.setPermMat(perm);
	}
	else{
		//*perm = prod(trans(*perm),*mPap.getPermMat());
		*perm = prod(*mPap.getPermMat(), trans(*perm));
		pap.setPermMat(perm);
	}

//#ifndef NDEBUG
	//tmp = prod(*pap.getPermMat(), *pap.getAdjMat());
	//tmp = prod(tmp, trans(*pap.getPermMat()));
	//matrix<double> tmp2(n,n);
	//tmp2 = prod(*mPap.getPermMat(), *mPap.getAdjMat());
	//tmp2 = prod(tmp2, trans(*mPap.getPermMat()));
	//for(int i=0;i<n;i++)
		//for(int j=0;j<n;j++){
			//I(tmp(i,j) - tmp2(i,j) < 0.001);
		//}
//#endif

	return pap;
}
JumbledAdjMatGen::~JumbledAdjMatGen()
{
}

