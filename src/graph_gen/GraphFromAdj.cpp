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
	typedef std::vector<double> Map;
	typedef bgl::adjacency_matrix<bgl::undirectedS,bgl::no_property,bgl::property<bgl::edge_weight_t, cost> > Graph;
	typedef bgl::property_map<Graph, bgl::edge_weight_t>::type WeightMap;
	typedef Graph::vertex_descriptor vertex_descriptor;
	typedef Graph::edge_descriptor edge_descriptor;
	typedef Graph::vertex_iterator vertex_iterator;

	struct DistAndHeight{
		vertex_descriptor A;
		vertex_descriptor B;
		Map distA;
		Map distB;
		Map height;
		double totalDist;
		DistAndHeight(){}
		DistAndHeight(int n) :A(0),B(n),distA(n),distB(n),height(n),totalDist(1)
		{
		}
	};

	Impl(const ProbAdjPerm&, unsigned int idx1, unsigned int idx2, int dim=3);
	Impl(const ProbAdjPerm&, int dim=3);
	~Impl();
	double  getDist(int refid, int i);
	double  getDist(DistAndHeight&, int i);
	double  getDistFromA(int refidx, int i);
	double  getDistFromB(int refidx, int i);
	void    setA(int refidx, int i);
	void    setB(int refidx, int i);
	unsigned int     getFarthestFrom(int refidx, const Map&);
	unsigned int     getFarthestFromA(int refidx);
	unsigned int     getFarthestFromB(int refidx);
	Graph            mG;
	WeightMap        mWeightMap;
	std::vector<DistAndHeight> mDistAndHeights;
};

double  
GraphFromAdj::Impl::getDist(DistAndHeight& dh, int idx){
	vertex_descriptor v = bgl::vertex(idx, mG);
	if(v == dh.A)  return 0.0;
	if(v == dh.B)  return dh.distA[dh.B];
	I(dh.A != dh.B); 
	I(dh.distA.size() > 0);
	I(dh.distA.size() == dh.distB.size());
	I(dh.distA.size() == dh.height.size());
	
	double a = dh.distA[v];
	double b = dh.distB[v];
#ifndef NDEBUG
	if(a!=a) cerr << "Warning: NaN found in distA!"<<endl;
	if(b!=b) cerr << "Warning: NaN found in distB!"<<endl;
	if(a+b<dh.totalDist-1E10)
		cerr << "Warning: Triangle Inequality does not hold: "<< (dh.totalDist-a-b)<<endl;
#endif
	double tmp = (dh.totalDist*dh.totalDist - a*a - b*b) / (-2*a*b);
	double gamma;
	if(0);
	else if(tmp<-1) gamma = M_PI;
	else if(tmp> 1) gamma = 0;
	else            gamma = acos( tmp );
	double s_gamma  = sin(gamma);
#ifndef NDEBUG
	if(s_gamma != s_gamma) cerr << "Warning: s_gamma is NaN"<<endl;
#endif
	      
	dh.height[v] = a*b*s_gamma/dh.totalDist; // height

	tmp = b/dh.totalDist*s_gamma;
#ifndef NDEBUG
	if(dh.totalDist==0) cerr << "Warning: totalDist = 0"<<endl;
	if(tmp != tmp) cerr << "Warning: tmp is NaN"<<endl;
#endif
	return (tmp>1)?0:(a * sqrt(1-tmp*tmp));
}

double
GraphFromAdj::Impl::getDist(int depth, int idx){
	return getDist(mDistAndHeights[depth], idx);
}

unsigned int
GraphFromAdj::Impl::getFarthestFrom(int refid, const Map& m){
	return distance(m.begin(),max_element(m.begin(),m.end()));
}

unsigned int
GraphFromAdj::Impl::getFarthestFromA(int refid){
	I(mDistAndHeights.size() > refid);
	return getFarthestFrom(refid, mDistAndHeights[refid].distA);
}
unsigned int
GraphFromAdj::Impl::getFarthestFromB(int refid){
	I(mDistAndHeights.size() > refid);
	return getFarthestFrom(refid, mDistAndHeights[refid].distB);
}
double
GraphFromAdj::Impl::getDistFromA(int refid, int i){ 
	I(mDistAndHeights.size() > refid);
	return mDistAndHeights[refid].distA[i]; }
double
GraphFromAdj::Impl::getDistFromB(int refid, int i){ 
	I(mDistAndHeights.size() > refid);
	return mDistAndHeights[refid].distB[i]; }
void
GraphFromAdj::Impl::setA(int refid, int i){ 
	I(mDistAndHeights.size() > refid);
	vertex_descriptor s = vertex(i, mG);
	DistAndHeight& dh = mDistAndHeights[refid];
	I(dh.distA.size()>0);
	dh.A = s;
	bgl::dijkstra_shortest_paths(mG, s, bgl::distance_map(&dh.distA[0]));
	dh.totalDist = dh.distA[dh.B];
}
void
GraphFromAdj::Impl::setB(int refid, int i){
	I(mDistAndHeights.size() > refid);
	vertex_descriptor s = vertex(i, mG);
	DistAndHeight& dh = mDistAndHeights[refid];
	I(dh.distB.size()>0);
	dh.B = s;
	bgl::dijkstra_shortest_paths(mG, s, bgl::distance_map(&dh.distB[0]));
	dh.totalDist = dh.distB[dh.A];
}


// Impl constructor
GraphFromAdj::Impl::Impl(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2, int dim)
	:mG(pap.getAdjMat()->size1())
	,mDistAndHeights(dim, DistAndHeight(pap.getAdjMat()->size1()))
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
	setA(0,idx1);
	setB(0,idx2);
}

// Impl constructor
GraphFromAdj::Impl::Impl(const ProbAdjPerm& pap, int dim)
	: mG(pap.getAdjMat()->size1())
	 ,mDistAndHeights(dim,DistAndHeight(pap.getAdjMat()->size1()))
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
}

// Impl destructor
GraphFromAdj::Impl::~Impl(){
}


/**********************************************************
 *          GraphFromAdj Interface                     *
 **********************************************************/
GraphFromAdj::GraphFromAdj(const ProbAdjPerm& pap, unsigned int idx1, unsigned int idx2, int dim)
	:mImpl(new Impl(pap,idx1,idx2,dim))
{
}
GraphFromAdj::GraphFromAdj(const ProbAdjPerm& pap, int dim)
	:mImpl(new Impl(pap,dim))
{
}

GraphFromAdj::~GraphFromAdj()
{
  // cleanup
}
double GraphFromAdj::getDist(int refid,int i)
{
	return mImpl->getDist(refid, i);
}
double GraphFromAdj::getTotalDist(int refid){
	return mImpl->mDistAndHeights[refid].totalDist;
}
double GraphFromAdj::getDistFromA(int refid, int i)
{ return mImpl->getDistFromA(refid, i); }
double GraphFromAdj::getDistFromB(int refid, int i)
{ return mImpl->getDistFromB(refid, i); }
void   GraphFromAdj::setA(int refid, int i)
{ mImpl->setA(refid, i);}
void   GraphFromAdj::setB(int refid, int i)
{ mImpl->setB(refid, i);}
unsigned int    GraphFromAdj::getFarthestFromA(int refid)
{ return mImpl->getFarthestFromA(refid); }
unsigned int    GraphFromAdj::getFarthestFromB(int refid)
{ return mImpl->getFarthestFromB(refid); }


