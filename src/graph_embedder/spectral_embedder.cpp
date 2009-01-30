#include <stdexcept>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <Serialization.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <factory/factory.h>
#include <configuration.hpp>
#include "spectral_embedder.hpp"

using namespace std;
namespace lapack = boost::numeric::bindings::lapack;
namespace ublas = boost::numeric::ublas;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>

SpectralGraphEmbedder::SpectralGraphEmbedder()
{
}
void 
SpectralGraphEmbedder::configure()
{
}

SpectralGraphEmbedder::cloud_type 
SpectralGraphEmbedder::operator()(const ProbAdjPerm& pap, int dim)
{
	ProbAdjLapPerm pap2(pap);
	pap2.calculateLaplacian();
	Laplacian::LaplacianT& L = *pap2.getLaplacian(); 
	int n = L.size1();
	Serialization ret(n);

	ublas::matrix<double,ublas::column_major> Eigv(L);
	ublas::vector<double> lambda(n);
	ublas::vector<double>::iterator best_lambda;

	lapack::syev( 'V', 'L', Eigv, lambda, lapack::minimal_workspace() );
	int start_lambda = 0;
	while(start_lambda<lambda.size()-1)
	{
		if(fabs(lambda(start_lambda)) > 1E-6)
			break;
		start_lambda++;
	}
	ublas::vector<double> v1 = ublas::column(Eigv,min(start_lambda + 0, n-1));
	ublas::vector<double> v2 = ublas::column(Eigv,min(start_lambda + 1, n-1));
	ublas::vector<double> v3 = ublas::column(Eigv,min(start_lambda + 2, n-1));
	ublas::vector<double> v4 = ublas::column(Eigv,min(start_lambda + 3, n-1));
	ublas::vector<double> v5 = ublas::column(Eigv,min(start_lambda + 4, n-1));
	//ublas::vector<int> ranks(n);
	//normalize_direction(v1,ranks,n);
	//normalize_direction(v2,ranks,n);
	//normalize_direction(v3,ranks,n);
	

	cloud_type cloud(n, point_type(dim,0));
	// create point cloud
	for(int i=0;i<n;i++){
		for(int d=0;d<min(dim,n-1);d++)
			cloud[i][d] = Eigv(i,min(start_lambda+d,n-1));
	}

	return cloud;
}


SpectralGraphEmbedder::~SpectralGraphEmbedder()
{
  // cleanup
}

namespace{ registerInFactory<GraphEmbedder, SpectralGraphEmbedder> registerBase("SpectralGraphEmbedder"); }
