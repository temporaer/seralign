#include <stdexcept>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <Serialization.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <factory/factory.h>
#include <configuration.hpp>
#include "heatkernel_embedder.hpp"
#include <nana.h>
#undef A

using namespace std;
namespace lapack = boost::numeric::bindings::lapack;
namespace ublas = boost::numeric::ublas;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>

HeatkernelGraphEmbedder::HeatkernelGraphEmbedder()
{
}
void 
HeatkernelGraphEmbedder::configure()
{
	mT = gCfg().getFloat("HeatkernelGraphEmbedder.t");
}

HeatkernelGraphEmbedder::cloud_type 
HeatkernelGraphEmbedder::operator()(const ProbAdjPerm& pap, int dim)
{
	AdjMat::AdjMatT& A = *pap.getAdjMat();
	int n = A.size1();

	// calculate Normalized Laplacian
	ublas::matrix<double,ublas::column_major> Eigv(n,n);
	ublas::vector<double> degs(n);
	for(int i=0;i<n;i++){
		double sum=0.0;
		for(int j=0;j<n;j++)
			sum += A(i,j);
		degs(i) = sum;
	}
	for(int i=0; i<n; i++){
		for(int j=i; j<n; j++){
			if(i==j)                 Eigv(i,j) = 1.0;
			else if(A(i,j)>0.00001)  {
				double x = degs(i)*degs(j);
				if(x>0)
					Eigv(i,j) = -1.0/sqrt(x);
				else{
					cerr <<"Warning: x<=0: "<<x<<endl;
					Eigv(i,j) = 0;
				}
			}
			else                     Eigv(i,j) = 0;
			Eigv(j,i) = Eigv(i,j);
		}
	}


	// eigen-analysis of normalized laplacian
	ublas::vector<double> lambda(n);
	lapack::syev( 'V', 'L', Eigv, lambda, lapack::minimal_workspace() );


	// calculate heat kernel
	ublas::matrix<double,ublas::column_major> H(n,n);
	for(int i=0; i<n; i++){
		for(int j=0; j<min(dim,n-1); j++){
			double sum=0.0;
			for(int k=0; k<n; k++){
				double e = exp(-lambda(k) * mT);
				double f = Eigv(k,i) * Eigv(k,j);
				sum +=  e*f ;
			}
			if(fabs(sum)<1e-10) sum=0;
			H(j,i) = H(i,j) = sum;
			if(sum!=sum)
				cerr << "Warning: sum is NaN"<<endl;
		}
	}

	cloud_type cloud(n, point_type(dim,0));
	// create point cloud
	for(int i=0;i<n;i++){
		for(int d=0;d<min(dim, n-1);d++)
			cloud[i][d] = H(i,d);
	}

	return cloud;
}


HeatkernelGraphEmbedder::~HeatkernelGraphEmbedder()
{
  // cleanup
}

namespace{ registerInFactory<GraphEmbedder, HeatkernelGraphEmbedder> registerBase("HeatkernelGraphEmbedder"); }
