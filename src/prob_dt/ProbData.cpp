#include <boost/numeric/ublas/matrix.hpp>
#include "ProbData.hpp"


namespace ublas = boost::numeric::ublas;
using namespace std;

// 2008-11-15 Hannes Schulz <mail at hannes-schulz dot de>

ProbAdjPerm::ProbAdjPerm()
{
}

ProbAdjPerm::~ProbAdjPerm()
{
  // cleanup
}


void ProbAdjLapPerm::calculateLaplacian()
{
	AdjMat::AdjMatT& adj = *getAdjMat();
	int n = adj.size1();
	ublas::vector<double> deg_vec(n);
	for(int i=0; i<n; i++){
		int sum = 0;
		for (int j=0;j<n;j++)
			sum += adj(i,j)>0 ? 1 : 0;
		deg_vec(i) = sum;
	}
	//diagonal_matrix<double> deg(n, deg_vec.data());
	boost::shared_ptr<Laplacian::LaplacianT> L_ptr( new Laplacian::LaplacianT(n,n) );
	Laplacian::LaplacianT& L = *L_ptr;
	for(int u=0;u<n;u++)
		for(int v=0;v<n;v++)
		{
			if(u==v && deg_vec(v) > 0.0000001)
				L(u,v) = 1;
			else if(adj(u,v)>0)
				L(u,v) = 1 / sqrt(deg_vec(u) * deg_vec(v));
			else
				L(u,v) = 0;
		}
	setLaplacian(L_ptr);
}
ProbAdjLapPerm::ProbAdjLapPerm(const ProbAdjPerm& p)
{
	setAdjMat(p.getAdjMat());
	setPermMat(p.getPermMat());
}

