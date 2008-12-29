// 2008-12-19 Hannes Schulz <mail at hannes-schulz dot de>

#include "GraphFromAdj.hpp"
#include <limits.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <ProbData.hpp>
#include <nana.h>
#undef A
#define V(X) #X << "=" << (X) << " "

namespace bgl = boost;
using namespace std;


/**********************************************************
 *          GraphFromAdj Implementation                *
 **********************************************************/
struct GraphFromAdj::Impl{
	typedef double cost;
	typedef bgl::adjacency_matrix<bgl::undirectedS,bgl::no_property,bgl::property<bgl::edge_weight_t, cost> > Graph;
	typedef bgl::property_map<Graph, bgl::edge_weight_t>::type WeightMap;
	typedef Graph::vertex_descriptor vertex_descriptor;
	typedef Graph::edge_descriptor edge_descriptor;
	typedef Graph::vertex_iterator vertex_iterator;

	Impl(const ProbAdjPerm&, unsigned int idx1, unsigned int idx2);
	~Impl();
	double  getDist(int i);
	double  getDistFromIdx1(int i);
	double  getDistFromIdx2(int i);
	void    setIdx1(int i);
	void    setIdx2(int i);
	unsigned int     getFarthestFromIdx1();
	unsigned int     getFarthestFromIdx2();
	Graph               mG;
	WeightMap           mWeightMap;
	double              mTotalDistance;
	pair<vertex_descriptor,vertex_descriptor> mPivot; 
	std::vector<double> mDistanceMapA;
	std::vector<double> mDistanceMapB;
};

double
GraphFromAdj::Impl::getDist(int idx){
	vertex_descriptor v = bgl::vertex(idx, mG);
	if(v == mPivot.first)  return 0.0;
	if(v == mPivot.second) return mTotalDistance;
	
	double a = mDistanceMapA[v];
	double b = mDistanceMapB[v];
#ifndef NDEBUG
	if(a+b<mTotalDistance)
		cerr << "Warning: Triangle Inequality does not hold!"<<endl;
#endif
	double tmp = (mTotalDistance*mTotalDistance - a*a - b*b) / (-2*a*b);
	double gamma;
	if(0);
	else if(tmp<-1) gamma = M_PI;
	else if(tmp> 1) gamma = 0;
	else            gamma = acos( tmp );
	double s_gamma = sin(gamma);
	tmp = b/mTotalDistance*s_gamma;
	double x = (tmp>1)?0:(a * sqrt(1-tmp*tmp));
	return x;
}

unsigned int
GraphFromAdj::Impl::getFarthestFromIdx1(){
	double maxDist = -1E9;
	unsigned int    maxIdx  = -1;
	for(unsigned int i=0;i<mDistanceMapA.size();i++) {
		vertex_descriptor v = bgl::vertex(i, mG);
		if(mDistanceMapA[v]>maxDist){
			maxDist = mDistanceMapA[v];
			maxIdx  = i;
		}
	}
	return maxIdx;
}
unsigned int
GraphFromAdj::Impl::getFarthestFromIdx2(){
	double maxDist = -1E9;
	unsigned int    maxIdx  = -1;
	for(unsigned int i=0;i<mDistanceMapB.size();i++) {
		vertex_descriptor v = bgl::vertex(i, mG);
		if(mDistanceMapA[v]>maxDist){
			maxDist = mDistanceMapB[v];
			maxIdx  = i;
		}
	}
	return maxIdx;
}
double
GraphFromAdj::Impl::getDistFromIdx1(int i){ return mDistanceMapA[i]; }
double
GraphFromAdj::Impl::getDistFromIdx2(int i){ return mDistanceMapB[i]; }
void
GraphFromAdj::Impl::setIdx1(int i){ 
	vertex_descriptor s = vertex(i, mG);
	mPivot.first = s;
	bgl::dijkstra_shortest_paths(mG, s, bgl::distance_map(&mDistanceMapA[0]));
	mTotalDistance = mDistanceMapA[mPivot.second];
}
void
GraphFromAdj::Impl::setIdx2(int i){
	vertex_descriptor t = vertex(i, mG);
	mPivot.second = t;
	bgl::dijkstra_shortest_paths(mG, t, bgl::distance_map(&mDistanceMapB[0]));
	mTotalDistance = mDistanceMapB[mPivot.first];
}


// Impl constructor
GraphFromAdj::Impl::Impl(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2)
	:
		 mG(pap.getAdjMat()->size1())
		,mDistanceMapA(pap.getAdjMat()->size1())
		,mDistanceMapB(pap.getAdjMat()->size1())
{
	AdjMat::AdjMatT& A = *pap.getAdjMat(); 
	unsigned int n=A.size1();
	mWeightMap = get(bgl::edge_weight,mG);
	for(unsigned int i=0;i<n;i++)
		for(unsigned int j=0;j<n;j++){
			edge_descriptor e; bool inserted;
			boost::tie(e, inserted) = bgl::add_edge(i,j,mG);
			if(A(i,j)<0.0000001) mWeightMap[e] = INT_MAX;
			else                 mWeightMap[e] = A(i,j);
		}
	setIdx1(idx1);
	setIdx2(idx2);
}

// Impl destructor
GraphFromAdj::Impl::~Impl(){
}


/**********************************************************
 *          GraphFromAdj Interface                     *
 **********************************************************/
GraphFromAdj::GraphFromAdj(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2)
	:mImpl(new Impl(pap,idx1,idx2))
{
}

GraphFromAdj::~GraphFromAdj()
{
  // cleanup
}
double GraphFromAdj::getDist(int i)
{
	return mImpl->getDist(i);
}
double GraphFromAdj::getTotalDist(){
	return mImpl->mTotalDistance;
}
double GraphFromAdj::getDistFromIdx1(int i)
{ return mImpl->getDistFromIdx1(i); }
double GraphFromAdj::getDistFromIdx2(int i)
{ return mImpl->getDistFromIdx2(i); }
void   GraphFromAdj::setIdx1(int i)
{ mImpl->setIdx1(i);}
void   GraphFromAdj::setIdx2(int i)
{ mImpl->setIdx2(i);}
unsigned int    GraphFromAdj::getFarthestFromIdx1()
{ return mImpl->getFarthestFromIdx1(); }
unsigned int    GraphFromAdj::getFarthestFromIdx2()
{ return mImpl->getFarthestFromIdx2(); }


