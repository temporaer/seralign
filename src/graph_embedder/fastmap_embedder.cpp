#include <GraphFromAdj.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <factory/factory.h>
#include "fastmap_embedder.hpp"

using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;

// 2009-01-17 Hannes Schulz <mail at hannes-schulz dot de>
//
struct GDistProjectedDB_PositionSorter{
	GDistProjectedDB_PositionSorter(const Serialization::PosT& p)
		:mPos(p)
	{
	}
	inline bool operator()(int i, int j)const{
		return mPos(i) < mPos(j);
	}
	const Serialization::PosT&  mPos;
};

template<class T, class U>
bool normalize_direction(ublas::vector<T>& v, ublas::vector<U>& w, int n){
	float center =(v[0]+v[n-1])/(float)n;
	float lever  = center;
	for(int i=0;i<n;i++)
		lever += center-v[i];
	if(lever < center){
		reverse(w.begin(),w.end());
		reverse(v.begin(),v.end());
		// muss 1-v[i] sein, nich einfach revert!
		v = ublas::scalar_vector<T>(n,1) - v;
		return true;
	}
	return false;
}

FastmapGraphEmbedder::FastmapGraphEmbedder()
{
}


void FastmapGraphEmbedder::configure()
{
	GraphEmbedder::configure();
	mFastmapTries     = gCfg().getInt("GDistProjectedDB.fastmap_tries");
}

FastmapGraphEmbedder::cloud_type 
FastmapGraphEmbedder::operator()(const ProbAdjPerm& pap, int dim)
{
	bool doNormalize = false;

	GraphFromAdj g(pap);
	int n = pap.getAdjMat()->size1();
	tuple<float, unsigned int, unsigned int> bestPair(0,0,1);
	for(int i=0; i< 2*mFastmapTries; i++){
		tuple<float, unsigned int, unsigned int> p = getFastmapQuality(g, 0, n*drand48());
		if(p.get<0>() > bestPair.get<0>())
			bestPair = p;
	}
	g.setA(0,bestPair.get<1>());
	g.setB(0,bestPair.get<2>());

	// convert graph to point-cloud:

	// create graph, calculate distances starting at 1st node
	Serialization::PosT  pos(n), pos2(n);
	Serialization::RankT ranks(n), ranks2(n);
	for(int i=0;i<n;i++){
		pos(i)   = g.getProjection(0,i);
#ifndef NDEBUG
		if(pos(i)!=pos(i))
			cerr <<"Warning: Distance is NaN for index "<<i<< endl;
#endif
		ranks(i) = i;
	}

	// order the columns in the adjmat according to their distance to the
	// first column in the original serialization
	sort(ranks.begin(),ranks.end(),GDistProjectedDB_PositionSorter(pos));
	// order the positions ascendingly (-->matching the order of the ranks)
	sort(pos.begin(),pos.end());
	double div = *max_element(pos.begin(),pos.end());
	if(doNormalize)
		pos /= div; // copy, otherwise max changed in between


	normalize_direction(pos, ranks, n);

	bestPair = tuple<float, unsigned int, unsigned int>(0,0,1);
	for(int i=n/2 - mFastmapTries; i<=n/2+mFastmapTries; i++){
		if(i<0)  continue;
		if(i>=n) continue;
		//tuple<float, unsigned int, unsigned int> p = getFastmapQuality(g, 1, ranks[i]);
		tuple<float, unsigned int, unsigned int> p = getFastmapQuality(g, 1, n*drand48());
		if(p.get<0>() > bestPair.get<0>())
			bestPair = p;
	}
	g.setA(1,bestPair.get<1>());
	g.setB(1,bestPair.get<2>());

	// project everything on line between the two vertices
	for(int i=0;i<n;i++){
		pos2(i)   = g.getProjection(1,i);
#ifndef NDEBUG
		if(pos2(i)!=pos2(i)) cerr <<"Warning: Distance is NaN for index "<<i<< endl;
#endif
		ranks2(i) = i;
	}
	//double div2 = *max_element(pos2.begin(),pos2.end());
	if(doNormalize)
		pos2 /= div; // copy! otherwise max changed while dividing.
	
	normalize_direction(pos2, ranks2, n);

	// create point cloud
	cloud_type cloud(n, point_type(dim,0));
	for(int i=0;i<n;i++){
		int idx1 = ranks[i];
		cloud[idx1][0] = pos[i];
		int idx2 = ranks2[i];
		cloud[idx2][1] = pos2[i];
	}
	return cloud;
}

boost::tuple<float,unsigned int,unsigned int>  
FastmapGraphEmbedder::getFastmapQuality(GraphFromAdj& g, int k, int seed)
{
	unsigned int tmp1, tmp2;
	g.setA(k, seed);
	tmp1 = g.getFarthestFromA(k);
	g.setA(k, tmp1);
	tmp2 = g.getFarthestFromA(k);
	float dist = g.getProjection(k,tmp2);
	//cout << "getFastmapQuality: dist="<<dist<<" tmp1="<<tmp1<<" tmp2="<<tmp2<<endl; 
	return make_tuple(dist,tmp1,tmp2) ;
}

FastmapGraphEmbedder::~FastmapGraphEmbedder()
{
  // cleanup
}

namespace{ registerInFactory<GraphEmbedder, FastmapGraphEmbedder> registerBase("FastmapGraphEmbedder"); }
