#include <memory>
#include <configuration.hpp>
#include <GraphFromAdj.hpp> 
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include "gdist_projected_db.hpp"
#include <gdist_seriation_gen.hpp>
#include <factory/factory.h>
#include <stats.hpp>
#include <nana.h>

using namespace std;
using namespace boost;
using namespace boost::numeric;
namespace lapack = boost::numeric::bindings::lapack;

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

GDistProjectedDB::TCloud
GDistProjectedDB::addSpectral(const string& pap_name, const ProbAdjPerm& pap)
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
	ublas::vector<double> v1 = ublas::column(Eigv,1);
	ublas::vector<double> v2 = ublas::column(Eigv,2);
	ublas::vector<double> v3 = ublas::column(Eigv,3);
	//ublas::vector<int> ranks(n);
	//normalize_direction(v1,ranks,n);
	//normalize_direction(v2,ranks,n);
	//normalize_direction(v3,ranks,n);
	

	TCloud cloud(n, point_type());
	// create point cloud
	for(int i=0;i<n;i++){
		cloud[i].id     = i;
		cloud[i].pos[0] = v1[i];
		cloud[i].pos[1] = v2[i];
		cloud[i].pos[2] = v3[i];
		cloud[i].pos[2] = 0;
	}

	TICP icp(mICPMaxIters,mICPLambda);
	mDB.push_back(icp);
	mDB.back().registerModel(cloud.begin(), cloud.end(), point_type::Translator());
	mIDs.push_back(pap_name);
	return cloud;
}
GDistProjectedDB::TCloud
GDistProjectedDB::add(const string& pap_name, const ProbAdjPerm& pap)
{
	// serialize once using
	//auto_ptr<SerGenAdj> seriation_gen    = genericFactory<SerGenAdj>::instance().create(mSeriationGenName);
	//if(!seriation_gen.get())
		//throw logic_error(string("Supplied SerGenAdj `") + mSeriationGenName + "' does not exist");
	//seriation_gen->configure();
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
	TCloud cloud(n, point_type());

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
	for(int i=0;i<n;i++){
		int idx1 = ranks[i];
		cloud[idx1].pos[0] = pos[i];
		cloud[idx1].id     = idx1;   // setting id once should be sufficient
		int idx2 = ranks2[i];
		cloud[idx2].pos[1] = pos2[i];
	}

	TICP icp(mICPMaxIters,mICPLambda);
	mDB.push_back(icp);
	mDB.back().registerModel(cloud.begin(), cloud.end(), point_type::Translator());
	mIDs.push_back(pap_name);
	return cloud;
}

void GDistProjectedDB::configure()
{
	mSeriationGenName = gCfg().getString("GDistProjectedDB.seriation_gen");
	mICPMaxIters      = gCfg().getInt("GDistProjectedDB.maxiter");
	mICPLambda        = gCfg().getFloat("GDistProjectedDB.lambda");
	mFastmapTries     = gCfg().getInt("GDistProjectedDB.fastmap_tries");
}

boost::tuple<float,unsigned int,unsigned int>  
GDistProjectedDB::getFastmapQuality(GraphFromAdj& g, int k, int seed)
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

