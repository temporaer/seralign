
#include	<algorithm>
#include	<numeric>
#include	<vector>
#include	<map>
#include	<boost/numeric/ublas/matrix.hpp>
#include	<boost/numeric/ublas/vector.hpp>
#include	<boost/numeric/ublas/matrix_proxy.hpp>

#include	<matlab_io.hpp>

#include "DegreeSort.hpp"
#include <nana.h>

using namespace std;
namespace ublas = boost::numeric::ublas;

// 2008-11-17 Hannes Schulz <mail at hannes-schulz dot de>

DegreeSort::DegreeSort()
{
}

DegreeSort::~DegreeSort()
{
  // cleanup
}

struct DegSort{
	map<int,double>& idx2deg;
	DegSort(map<int,double>& i2d)
		:idx2deg(i2d){}
	inline bool operator()(int a, int b){
		return idx2deg[a] < idx2deg[b];
	}
};

void DegreeSort::sort(ProbAdjPerm& pap)
{
	AdjMat::AdjMatT& adj    = *pap.getAdjMat();
	int n = adj.size1();
	boost::shared_ptr<PermMat::PermMatT> perm(new PermMat::PermMatT(n,n));
	*perm = boost::numeric::ublas::zero_matrix<double>(n,n);

	std::vector<double> idxs(adj.size2());
	for(unsigned int i=0;i<idxs.size();i++)
		idxs[i] = (double) i;

	// map index to degree
	std::map<int,double> idx2deg;
	for(unsigned int i=0;i<adj.size2();i++){
		ublas::matrix_column<AdjMat::AdjMatT> col(ublas::column(adj,i));
		idx2deg[i] = std::accumulate(col.begin(),col.end(),0);
	}

	// sort the ids by degree
	std::sort(idxs.begin(),idxs.end(),DegSort(idx2deg));
#ifndef NDEBUG
	for(int i=1;i<n;i++) { I(idx2deg[idxs[i-1]] <= idx2deg[idxs[i]]); }
#endif

	// record ordering in perm mat
	for(unsigned int i=0;i<idxs.size();i++)
		(*perm)(i,idxs[i]) = 1; // projects from idxs[i] to i, that is, to original pos
	
	boost::shared_ptr<AdjMat::AdjMatT> adj_new(new AdjMat::AdjMatT(n,n));

	noalias(*adj_new) = prod(*perm, adj);
	*adj_new = prod(*adj_new, trans(*perm));

#ifndef NDEBUG
	ublas::vector<int> idxs_vec(n);
	for(unsigned int i=0;i<idxs.size();i++)
		idxs_vec(idxs[i]) = i;
	ublas::vector<int> v(n);
	for(int i=0;i<n;i++) v(i) = i;
	ublas::vector<int> p = prod(trans(*perm),v);
	ublas::vector<int> tmp = p - idxs_vec;
	for(int i=0;i<n;i++) { I(tmp(i)==0); }
#endif

	// save perm mat in problem
	if(pap.getPermMat() == NULL)
		pap.setPermMat(perm);
	else{
		PermMat::PermMatT tmp(n,n);
		noalias(tmp) = prod(*perm, *pap.getPermMat());
		*pap.getPermMat() = prod(tmp, trans(*perm));
	}
	// save adj mat in problem
	pap.setAdjMat(adj_new);
}

