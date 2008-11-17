
#include	<algorithm>
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
		return idx2deg[a] > idx2deg[b];
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
		double sum=0;
		for(unsigned int j=0;j<adj.size1();j++)
			sum += adj(i,j);
		idx2deg[i] = sum;
	}

	// sort the ids by degree
	std::sort(idxs.begin(),idxs.end(),DegSort(idx2deg));

	// record ordering in perm mat
	for(unsigned int i=0;i<idxs.size();i++)
		(*perm)(i,idxs[i]) = 1; // projects from idxs[i] to i, that is, to original pos
	
	boost::shared_ptr<AdjMat::AdjMatT> adj_new(new AdjMat::AdjMatT(n,n));
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			(*adj_new)(i,j) = adj(idxs[i],idxs[j]);

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
	pap.setPermMat(perm);
	// save adj mat in problem
	pap.setAdjMat(adj_new);
}

