#include <memory>
#include <configuration.hpp>
#include <GraphFromAdj.hpp> 
#include "gdist_projected_db.hpp"
#include <gdist_seriation_gen.hpp>
#include <factory/factory.h>

using namespace std;
using namespace boost;
using namespace boost::numeric;

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

void GDistProjectedDB::add(const string& pap_name, const ProbAdjPerm& pap)
{
	// serialize once using
	auto_ptr<SerGenAdj> seriation_gen    = genericFactory<SerGenAdj>::instance().create(mSeriationGenName);
	if(!seriation_gen.get())
		throw logic_error(string("Supplied SerGenAdj `") + mSeriationGenName + "' does not exist");
	seriation_gen->configure();

	Serialization ser = (*seriation_gen)(pap);

	// drop out if serialization empty
	int n = ser.getRanks().size();
	if(n==0) return;
	
	// convert graph to point-cloud:
	std::vector<ublas::vector<double> > cloud(n, ublas::scalar_vector<double>(3,0));

	// create graph, calculate distances starting at 1st node
	GraphFromAdj g(pap,ser.getRanks()[0], ser.getRanks()[n-1]);
	Serialization::PosT  pos(n), pos2(n);
	Serialization::RankT ranks(n), ranks2(n);
	for(int i=0;i<n;i++){
		pos(i)   = g.getDist(i);
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


	// try to find "canonical" direction: 
	float center =(pos[0]+pos[n-1])/pos.size();
	float lever  = center;
	for(int i=0;i<n;i++)
		lever += center-pos[i];
	if(lever < center){
		reverse(ranks.begin(),ranks.end());
		reverse(pos.begin(),pos.end());
	}


	tuple<float, unsigned int, unsigned int> bestPair(0,0,1);
	for(int i=n/2 - mFastmapTries; i<=n/2+mFastmapTries; i++){
		if(i<0)  continue;
		if(i>=n) continue;
		tuple<float, unsigned int, unsigned int> p = getFastmapQuality(g, ranks[i]);
		if(p.get<0>() > bestPair.get<0>())
			bestPair = p;
	}
	g.setIdx1(bestPair.get<1>());
	g.setIdx2(bestPair.get<2>());

	// project everything on line between the two vertices
	for(int i=0;i<n;i++){
		pos2(i)   = g.getDist(i);
#ifndef NDEBUG
		if(pos2(i)!=pos2(i)) cerr <<"Warning: Distance is NaN for index "<<i<< endl;
#endif
		ranks2(i) = i;
	}
	
	// try to find "canonical" direction: 
	float center2 =(pos2[0]+pos2[n-1])/pos2.size();
	float lever2  = center2;
	for(int i=0;i<n;i++)
		lever2 += center2-pos2[i];
	if(lever2 < center2){
		reverse(ranks2.begin(),ranks2.end());
		reverse(pos2.begin(),pos2.end());
	}

	// create point cloud
	for(int i=0;i<n;i++){
		int idx1 = ranks[i];
		cloud[idx1][0] = pos[i];
		int idx2 = ranks2[i];
		cloud[idx2][1] = pos2[i];
	}

	TICP icp(mICPMaxIters,mICPLambda);
	mDB.push_back(icp);
	mDB.back().registerModel(cloud.begin(), cloud.end());
	mIDs.push_back(pap_name);
}

void GDistProjectedDB::configure()
{
	mSeriationGenName = gCfg().getString("GDistProjectedDB.seriation_gen");
	mICPMaxIters      = gCfg().getInt("GDistProjectedDB.maxiter");
	mICPLambda        = gCfg().getFloat("GDistProjectedDB.lambda");
	mFastmapTries     = gCfg().getInt("GDistProjectedDB.fastmap_tries");
}

boost::tuple<float,unsigned int,unsigned int>  
GDistProjectedDB::getFastmapQuality(GraphFromAdj& g, int seed)
{
	unsigned int tmp1, tmp2;
	g.setIdx1(seed);
	tmp1 = g.getFarthestFromIdx1();
	g.setIdx1(tmp1);
	tmp2 = g.getFarthestFromIdx1();
	return make_tuple(g.getDistFromIdx1(tmp1),tmp1,tmp2) ;
}

