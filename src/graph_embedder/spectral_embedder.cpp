#include <stdexcept>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <Serialization.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <factory/factory.h>
#include <configuration.hpp>
#include "spectral_embedder.hpp"
#include <nana.h>

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
SpectralGraphEmbedder::operator()(const ProbAdjPerm& pap, unsigned int dim)
{
	ProbAdjLapPerm pap2(pap);
	pap2.calculateLaplacian();
	Laplacian::LaplacianT& L = *pap2.getLaplacian(); 
	unsigned int n = L.size1();
	Serialization ret(n);

	ublas::matrix<double,ublas::column_major> Eigv(L);
	ublas::vector<double> lambda(n);
	ublas::vector<double>::iterator best_lambda;

	lapack::syev( 'V', 'L', Eigv, lambda, lapack::minimal_workspace() );
	unsigned int start_lambda = 0;
#ifndef NDEBUG
	for(unsigned int i=1; i<n; ++i){
		I(lambda(i) >= lambda(i-1));
	}
#endif
#if 0
	while(start_lambda<lambda.size()-1)
	{
		if(fabs(lambda(start_lambda)) > 1E-6)
			break;
		start_lambda++;
	}
#else
	//use last dim vectors -- does not work, it seems
	//start_lambda = max(0,(int)n-(int)dim-1);
#endif


#define NORMALIZE_EV 1
#if NORMALIZE_EV
	//Eigv
#endif

	cloud_type cloud(n, point_type(dim,0));
	// create point cloud
	for(unsigned int i=0;i<n;i++){
		for(unsigned int d=0;d<min(dim,n-1);d++)
			cloud[i][d] = Eigv(i,min(start_lambda+d,n-1));
			//cloud[i][d] = fabs(Eigv(i,min(start_lambda+d,n-1)));
	}

	return cloud;
}


SpectralGraphEmbedder::~SpectralGraphEmbedder()
{
  // cleanup
}

namespace{ registerInFactory<GraphEmbedder, SpectralGraphEmbedder> registerBase("SpectralGraphEmbedder"); }
