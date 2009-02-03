#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/bindings/lapack/lapack.hpp>
//#include <Serialization.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <factory/factory.h>
#include <configuration.hpp>
#include "feature_embedder.hpp"
#include <nana.h>

using namespace std;
//namespace lapack = boost::numeric::bindings::lapack;
namespace ublas = boost::numeric::ublas;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>

FeatureGraphEmbedder::FeatureGraphEmbedder()
{
}
void 
FeatureGraphEmbedder::configure()
{
}

void    FeatureGraphEmbedder::setAdjMatGen(boost::shared_ptr<AdjMatGen> g, unsigned int dim)
{
	mGen=g;
	mMats.clear();
	unsigned int n = g->getNumFeatures();
	while(n > dim){
		unsigned int m = max((int)dim, (int)(log(n)+0.5));
		ublas::matrix<double> proj(m,n);
		for(unsigned int i=0;i<m;i++)
			for(unsigned int j=0;j<n;j++)
				proj(i,j) = drand48() >= 0.5 ? 1 : -1;
		mMats.push_back(proj);
		n = m;
	}
	mOS.close();
	mOS.open("/tmp/cloud.dat");
}

FeatureGraphEmbedder::cloud_type 
FeatureGraphEmbedder::operator()(const ProbAdjPerm& pap, unsigned int dim)
{
	cloud_type cloud;
	ublas::vector<double> fw = mGen->getFeatureWeights();
	for(unsigned int i=0; i< pap.getAdjMat()->size1(); i++){
		AdjMatGen::feature_t f = element_prod(fw, mGen->getFeatures(i, pap.getBackground()));
		BOOST_FOREACH(const ublas::matrix<double>& p, mMats){
			f = ublas::prod(p,f);
		}
		cloud.push_back(f);
	}
	BOOST_FOREACH(const ublas::vector<double>& p, cloud){
		copy(p.begin(), p.end(), ostream_iterator<double>(mOS, "\t"));
		mOS << mGen->getClassID(pap.getBackground());
		mOS << endl;
	}
	mOS << endl;
	return cloud;
}


FeatureGraphEmbedder::~FeatureGraphEmbedder()
{
  // cleanup
}

namespace{ registerInFactory<GraphEmbedder, FeatureGraphEmbedder> registerBase("FeatureGraphEmbedder"); }
